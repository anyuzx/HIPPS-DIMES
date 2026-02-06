"""Optimization and dynamics models for HIPPS-DIMES."""

from .numerics import *  # noqa: F401,F403
from .numerics import _a2dmap_theory_gpu

class Optimize:
    """
    Optimizer for finding connectivity matrix A that matches target distance constraints.
    
    Uses Maximum Entropy principle to find the Gaussian polymer model that best
    matches the target mean squared distance map.
    
    Optimization Methods
    --------------------
    - 'IS' (Iterative Scaling): Default method. Multiplicative updates based on
      log-ratio of current vs target distances. Generally stable and reliable.
    - 'GD' (Gradient Descent): Additive updates with Nesterov momentum.
      Requires smaller learning rates (e.g., 1e-8).
    - 'DI' (Direct Inversion): One-shot pseudo-inverse of covariance matrix.
      Fast but may not satisfy all constraints; best for valid Euclidean distance matrices.
    
    Convergence Acceleration
    ------------------------
    **Momentum with Nesterov** (RECOMMENDED for fastest convergence):
    - Nesterov uses "look-ahead" correction to prevent overshooting
    - Enables higher momentum values (0.95) that would diverge with standard momentum
    - BEST setting: momentum=0.95, nesterov=True (~50% faster than momentum=0.9)
    
    **Standard Momentum** (fallback option):
    - Uses Polyak's heavy ball method
    - Safe with momentum=0.9, provides ~2-4x speedup over baseline
    
    Examples
    --------
    >>> # Standard optimization (no acceleration)
    >>> opt = Optimize(ddmap_target)
    >>> loss, entropy, dmap, A = opt.run(1000, learning_rate=10.0, method='IS')
    
    >>> # With Nesterov momentum (RECOMMENDED - fastest)
    >>> loss, entropy, dmap, A = opt.run(1000, learning_rate=10.0, method='IS', momentum=0.95, nesterov=True)
    
    >>> # With standard momentum (fallback)
    >>> loss, entropy, dmap, A = opt.run(1000, learning_rate=10.0, method='IS', momentum=0.9)
    """
    
    def __init__(self, ddmap_target, connectivity_matrix=None, use_gpu=False, gpu_float32=False):
        # ddmap_target is the targeted matrix we would like to match
        # note that ddmap_taret is the mean SQUARED distance matrix, not mean distance matrix
        self.ddmap_target = ddmap_target

        # get the size of system
        self.n = ddmap_target.shape[0]

        if connectivity_matrix is None:
            # initialize the connectivity matrix
            # here the connectivity matrix is initialized as a simple rouse chain whose spring constant is determined such\
            # that its radius of gyration is close to the target
            # we need to filter out both NaN and inf entries
            rg2 = .5 * np.nanmean(self.ddmap_target[~np.isinf(self.ddmap_target)])
            k = self.n / (4. * rg2)
            self.A = construct_connectivity_matrix_rouse(self.n, k)
        else:
            self.A = connectivity_matrix

        # initialize the loss
        self.loss = None

        # initialize the entropy
        self.entropy = None

        # Optional off-diagonal mask: True means keep/update A[i,j]; False means freeze at 0
        # Diagonal entries are always recomputed by a2a()
        self.edge_mask = None
        
        # GPU support
        self.use_gpu = use_gpu and is_gpu_available()
        if use_gpu and not is_gpu_available():
            console.print("[yellow]Warning: use_gpu=True but CuPy is not available. Falling back to CPU.[/yellow]")
        
        if self.use_gpu:
            # Optional float32 mode on GPU (faster on many GPUs, but can slightly change numerics)
            self.gpu_float32 = bool(gpu_float32)
            self._gpu_dtype = cp.float32 if self.gpu_float32 else cp.float64
            # Move data to GPU
            self._A_gpu = cp.asarray(self.A, dtype=self._gpu_dtype)
            self._ddmap_target_gpu = cp.asarray(self.ddmap_target, dtype=self._gpu_dtype)
            self._velocity_gpu = None
            # Cached GPU masks (to avoid per-iteration cp.asarray allocations)
            self._edge_mask_gpu = None
            self._freeze_mask_gpu = None
            # Theta for GD GPU path
            self._theta_gpu = None
            # Loss/entropy scalars on GPU (avoid per-iteration device->host sync)
            self._loss_gpu = None
            self._entropy_gpu = None
            # How often to sync loss/entropy for tqdm display (syncing every iter can dominate runtime)
            self.gpu_display_every = 10
            cp.cuda.Stream.null.synchronize()
        else:
            self.gpu_float32 = False
            self._gpu_dtype = None
            self._edge_mask_gpu = None
            self._freeze_mask_gpu = None
            self._theta_gpu = None
            self._loss_gpu = None
            self._entropy_gpu = None
            self.gpu_display_every = None

    def set_edge_mask(self, edge_mask):
        """Set/update the off-diagonal mask for optimization.

        Parameters
        ----------
        edge_mask : (n, n) array_like of bool or None
            True => keep/allow this off-diagonal entry A[i,j] to be optimized (i!=j).
            False => freeze A[i,j]=0 (i!=j). The diagonal is ignored and recomputed by a2a().
            If None, no masking is applied (all off-diagonals are free).
        """
        if edge_mask is None:
            self.edge_mask = None
            if self.use_gpu:
                self._edge_mask_gpu = None
                self._freeze_mask_gpu = None
            return

        m = np.array(edge_mask, dtype=bool, copy=True)
        if m.shape != (self.n, self.n):
            raise ValueError(f"edge_mask must have shape ({self.n},{self.n}), got {m.shape}")

        # Symmetrize and ignore diagonal
        m = np.logical_and(m, m.T)
        np.fill_diagonal(m, False)
        self.edge_mask = m
        # Cache GPU copies of masks (mask is typically constant across iterations)
        if self.use_gpu:
            self._edge_mask_gpu = cp.asarray(self.edge_mask)
            self._freeze_mask_gpu = ~self._edge_mask_gpu
            cp.fill_diagonal(self._freeze_mask_gpu, False)

    def _freeze_masked_edges(self):
        """Force masked (disallowed) off-diagonals to be exactly zero."""
        if self.edge_mask is None:
            return
        freeze = ~self.edge_mask
        np.fill_diagonal(freeze, False)
        self.A[freeze] = 0.0

    def run_masked(self, epoch, edge_mask, ddmap_target_masked=None, general_method='optimization', **kwargs):
        """Run optimization with a fixed off-diagonal mask, optimized for *sparse* constraints.

        Replaces the previous wrapper that called `run()` (which recomputed the full
        distance map via an eigendecomposition at every iteration).

        When the number of constrained pairs is small (typical in minimax / greedy growing),
        we compute only the needed mean-squared distances using grounded Laplacian solves:

            ddmap(i,j) = <||x_i-x_j||^2> = 3 * (e_i-e_j)^T (-A)^+ (e_i-e_j)

        Parameters
        ----------
        epoch : int
            Number of iterations.
        edge_mask : (n,n) bool
            Off-diagonal mask; True => allow/update A_ij for i!=j, False => freeze A_ij=0.
        ddmap_target_masked : (n,n) float or None
            If provided, overrides self.ddmap_target. Use NaN/inf on unconstrained pairs.
        general_method : str
            Only 'optimization' is supported here.
        kwargs :
            learning_rate : float
            lamd : float
            reg : {'L1','L2'}
            method : {'IS'}  (other methods fall back to the full `run()` implementation)
            enforce_nonnegative_connectivity_matrix : bool
            pin : int (default 0)  -- grounded node for Laplacian solve
            jitter : float (default 1e-10) -- diagonal jitter for numerical stability
        """
        if ddmap_target_masked is not None:
            if ddmap_target_masked.shape != (self.n, self.n):
                raise ValueError("ddmap_target_masked must match the current system size")
            self.ddmap_target = ddmap_target_masked

        if general_method != 'optimization':
            self.set_edge_mask(edge_mask)
            loss_array, entropy_array, dmap_maxent, A_final, _ = self.run(epoch, general_method=general_method, **kwargs)
            return loss_array, entropy_array, dmap_maxent, A_final

        method = kwargs.get('method', 'IS')
        if method != 'IS':
            self.set_edge_mask(edge_mask)
            loss_array, entropy_array, dmap_maxent, A_final, _ = self.run(epoch, general_method=general_method, **kwargs)
            return loss_array, entropy_array, dmap_maxent, A_final

        learning_rate = float(kwargs.get('learning_rate', 1.0))
        lamd = float(kwargs.get('lamd', 0.0))
        reg = str(kwargs.get('reg', 'L2')).upper()
        enforce_nonnegative = bool(kwargs.get('enforce_nonnegative_connectivity_matrix', False))
        pin = int(kwargs.get('pin', 0))
        jitter = float(kwargs.get('jitter', 1e-10))

        self.set_edge_mask(edge_mask)

        # -------- Determine constrained pairs (i<j) --------
        tgt = np.array(self.ddmap_target, dtype=float, copy=False)
        if self.edge_mask is None:
            loss_array, entropy_array, dmap_maxent, A_final, _ = self.run(epoch, general_method=general_method, **kwargs)
            return loss_array, entropy_array, dmap_maxent, A_final

        finite = np.isfinite(tgt) & (tgt > 0.0)
        allow = self.edge_mask & finite
        ii, jj = np.where(np.triu(allow, 1))
        m = ii.size
        if m == 0:
            self._freeze_masked_edges()
            self.A = a2a(0.5 * (self.A + self.A.T), fill_negative=enforce_nonnegative)
            dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)
            # Compute entropy
            K = -self.A
            eigvals_K = np.linalg.eigvalsh(K)
            entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)
            return [], [entropy], dmap_maxent, self.A

        dd_target_edges = tgt[ii, jj].astype(float)

        # -------- Precompute grounded incidence matrix B for constrained pairs --------
        if not (0 <= pin < self.n):
            raise ValueError(f"pin must be in [0,{self.n-1}]")

        keep = np.ones(self.n, dtype=bool)
        keep[pin] = False

        pos = -np.ones(self.n, dtype=int)
        pos[np.where(keep)[0]] = np.arange(self.n - 1)

        B = np.zeros((self.n - 1, m), dtype=float)
        cols = np.arange(m)

        gi = pos[ii]
        gj = pos[jj]
        valid_i = gi >= 0
        valid_j = gj >= 0
        if np.any(valid_i):
            B[gi[valid_i], cols[valid_i]] += 1.0
        if np.any(valid_j):
            B[gj[valid_j], cols[valid_j]] -= 1.0

        # -------- Main loop --------
        loss_array = []
        entropy_array = []
        with trange(epoch, desc="Performing masked optimization (sparse)", unit="iteration") as pbar:
            for t in pbar:
                # Grounded Laplacian Lg = (-A) with pinned node removed
                A_sym = 0.5 * (self.A + self.A.T)
                L = -A_sym
                Lg = L[np.ix_(keep, keep)].copy()
                if jitter > 0.0:
                    Lg.flat[::Lg.shape[0] + 1] += jitter

                # Solve Lg X = B for all constrained edges at once
                try:
                    c, lower = scipy.linalg.cho_factor(Lg, lower=True, check_finite=False, overwrite_a=False)
                    X = scipy.linalg.cho_solve((c, lower), B, check_finite=False)
                except np.linalg.LinAlgError:
                    # Fallback: eigen-based pseudo-inverse on the grounded system
                    w, V = np.linalg.eigh(Lg)
                    w_inv = np.zeros_like(w)
                    good = w > jitter
                    w_inv[good] = 1.0 / w[good]
                    X = (V * w_inv) @ (V.T @ B)

                # ddmap on constrained edges: dd = 3 * b^T Lg^{-1} b
                rho = np.sum(B * X, axis=0)
                dd_cur_edges = 3.0 * np.maximum(rho, 0.0)

                # Iterative-scaling gradient on constrained edges
                with np.errstate(divide='ignore', invalid='ignore'):
                    log_ratio = np.log(dd_cur_edges / dd_target_edges)
                log_ratio = np.nan_to_num(log_ratio, posinf=0.0, neginf=0.0)

                fhash = 0.5 * float(np.sum(dd_cur_edges))
                if fhash <= 0.0:
                    fhash = 1.0

                g = log_ratio

                # Regularization (proximal L1 handled after update)

                g = g / fhash

                # Update only allowed off-diagonals
                self.A[ii, jj] += learning_rate * g
                self.A[jj, ii] += learning_rate * g

                # Cleanup + enforce mask + rebuild diagonal
                if lamd > 0.0:
                    if reg == 'L2':
                        shrink = 1.0 / (1.0 + learning_rate * lamd)
                        self.A = self.A * shrink
                    elif reg == 'L1':
                        thresh = learning_rate * lamd
                        self.A = np.sign(self.A) * np.maximum(np.abs(self.A) - thresh, 0.0)

                self.A = np.nan_to_num(self.A)
                self.A = 0.5 * (self.A + self.A.T)
                self._freeze_masked_edges()
                self.A = a2a(self.A, fill_negative=enforce_nonnegative)

                # Loss on constrained edges only
                with np.errstate(divide='ignore', invalid='ignore'):
                    rel = (dd_cur_edges - dd_target_edges) / dd_target_edges
                    loss = float(np.nanmean(rel * rel) ** 0.5)
                self.loss = loss
                
                # Compute entropy using eigenvalues of -A (K = -A)
                K = -self.A
                eigvals_K = np.linalg.eigvalsh(K)
                entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)
                self.entropy = entropy
                
                pbar.set_postfix(loss=self.loss, entropy=self.entropy if self.entropy is not None else np.nan)
                loss_array.append(self.loss)
                entropy_array.append(self.entropy if self.entropy is not None else np.nan)

        # Compute full mean distance map once at the end (compat)
        dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)
        return loss_array, entropy_array, dmap_maxent, self.A
    def __compute_loss(self, ddmap_t):
        with np.errstate(divide='ignore', invalid='ignore'):
            loss = np.nanmean(
                np.power((ddmap_t - self.ddmap_target)/self.ddmap_target, 2.)) ** .5
        return loss

    def __update_parameter(self, t, learning_rate, lamd=0.0, reg='l2', method='IS', enforce_nonnegative_connectivity_matrix=False, momentum=0.0, nesterov=False):
        # updating using Iterative Scaling or Gradient Descent
        
        if self.use_gpu and method == 'IS':
            # GPU path for IS using CuPy
            return self.__update_parameter_gpu(t, learning_rate, lamd, reg, enforce_nonnegative_connectivity_matrix, momentum, nesterov)
        elif self.use_gpu and method == 'GD':
            # GPU path for GD using CuPy
            return self.__update_parameter_gpu_gd(t, learning_rate, lamd, reg, enforce_nonnegative_connectivity_matrix)

        # CPU path (original implementation)
        # compute the mean squared distance matrix at current iteration step
        # Also get eigenvalues of A to reuse for entropy computation (avoiding double eigendecomposition)
        dmap_t, eigvals_A = a2dmap_theory(self.A, force_positive_definite=True, return_eigenvalues=True)
        ddmap_t = ((3. * np.pi) / 8.) * np.power(dmap_t, 2.)
        # compute the ratio between the current value and the target
        with np.errstate(divide='ignore', invalid='ignore'):
            compare_ratio = ddmap_t / self.ddmap_target
        # compute the prefactor for iterative scaling
        fhash = np.nansum(ddmap_t) / 2.

        if method == 'IS':
            # compute the gradient
            gradient_t = np.nan_to_num(
                np.log(compare_ratio), posinf=0., neginf=0.) / fhash

            # enforce symmetry and apply optional off-diagonal mask
            gradient_t = 0.5 * (gradient_t + gradient_t.T)
            if self.edge_mask is not None:
                gradient_t *= self.edge_mask

            # Apply update with optional momentum (standard or Nesterov)
            if momentum > 0.0:
                # Initialize velocity on first iteration
                if t == 0:
                    self.velocity = np.zeros_like(self.A)
                # Update velocity
                self.velocity = momentum * self.velocity + gradient_t
                
                if nesterov:
                    # Nesterov Accelerated Gradient: use look-ahead correction
                    # Update = gradient + momentum * velocity (after velocity update)
                    # This is equivalent to computing gradient at A + momentum * velocity
                    self.A = self.A + learning_rate * (gradient_t + momentum * self.velocity)
                else:
                    # Standard momentum (Polyak's heavy ball)
                    self.A = self.A + learning_rate * self.velocity
            else:
                # Standard update without momentum
                self.A = self.A + learning_rate * gradient_t
            if lamd > 0.0:
                if reg == 'L2':
                    shrink = 1.0 / (1.0 + learning_rate * lamd)
                    self.A = self.A * shrink
                elif reg == 'L1':
                    thresh = learning_rate * lamd
                    self.A = np.sign(self.A) * np.maximum(np.abs(self.A) - thresh, 0.0)
        elif method == 'GD':
            if t == 0:
                self.theta = np.copy(self.A)

            # compute the gradient
            gradient_t = (ddmap_t - self.ddmap_target)

            # enforce symmetry and apply optional off-diagonal mask
            gradient_t = 0.5 * (gradient_t + gradient_t.T)
            if self.edge_mask is not None:
                gradient_t *= self.edge_mask

            # perform Nesterov update rule
            # gradient descent state
            theta_previous = np.copy(self.theta)

            #self.theta = self.A + np.maximum(np.minimum(learning_rate * gradient_t, step_cap),-step_cap)
            self.theta = self.A + learning_rate * gradient_t

            # if momentum_rate == None:
            #    momentum_rate = t/(t+3)

            # update the connectivity matrix
            self.A = self.theta + (t/(t+3)) * (self.theta - theta_previous)

            if lamd > 0.0:
                if reg == 'L2':
                    shrink = 1.0 / (1.0 + learning_rate * lamd)
                    self.A = self.A * shrink
                    self.theta = self.theta * shrink
                elif reg == 'L1':
                    thresh = learning_rate * lamd
                    self.A = np.sign(self.A) * np.maximum(np.abs(self.A) - thresh, 0.0)
                    self.theta = np.sign(self.theta) * np.maximum(np.abs(self.theta) - thresh, 0.0)

        # convert all nan to zero
        self.A = np.nan_to_num(self.A)

        # keep symmetry (numerical) and freeze masked edges before rebuilding diagonals
        self.A = 0.5 * (self.A + self.A.T)
        self._freeze_masked_edges()

        self.A = a2a(self.A, fill_negative=enforce_nonnegative_connectivity_matrix)
        # project to be negative semidefinite
        #self.A = nearestNSD(self.A, 0.0)

        # compute the loss
        self.loss = self.__compute_loss(ddmap_t)
        
        # Compute entropy reusing eigenvalues from a2dmap_theory
        # Note: eigenvalues of K = -A are simply -eigvals_A (no need for second eigendecomposition)
        eigvals_K = -eigvals_A
        self.entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)

    def __update_parameter_noisy(self, t, learning_rate, gaussian_noise_variance, method='IS', enforce_nonnegative_connectivity_matrix=False, momentum=0.0, nesterov=False):
        """
        Update parameters with Gaussian-noise regularization on constraints.

        Uses a proximal L2 shrink step after the standard IS/GD update:
        A <- (A + lr * grad) / (1 + lr * gaussian_noise_variance)
        """
        if self.use_gpu and method == 'IS':
            return self.__update_parameter_noisy_gpu(t, learning_rate, gaussian_noise_variance, enforce_nonnegative_connectivity_matrix, momentum, nesterov)
        elif self.use_gpu and method == 'GD':
            return self.__update_parameter_noisy_gpu_gd(t, learning_rate, gaussian_noise_variance, enforce_nonnegative_connectivity_matrix)

        dmap_t, eigvals_A = a2dmap_theory(self.A, force_positive_definite=True, return_eigenvalues=True)
        ddmap_t = ((3. * np.pi) / 8.) * np.power(dmap_t, 2.)
        with np.errstate(divide='ignore', invalid='ignore'):
            compare_ratio = ddmap_t / self.ddmap_target
        fhash = np.nansum(ddmap_t) / 2.

        if method == 'IS':
            gradient_t = np.nan_to_num(np.log(compare_ratio), posinf=0., neginf=0.) / fhash
            gradient_t = 0.5 * (gradient_t + gradient_t.T)
            if self.edge_mask is not None:
                gradient_t *= self.edge_mask

            if momentum > 0.0:
                if t == 0:
                    self.velocity = np.zeros_like(self.A)
                self.velocity = momentum * self.velocity + gradient_t
                if nesterov:
                    self.A = self.A + learning_rate * (gradient_t + momentum * self.velocity)
                else:
                    self.A = self.A + learning_rate * self.velocity
            else:
                self.A = self.A + learning_rate * gradient_t
        elif method == 'GD':
            if t == 0:
                self.theta = np.copy(self.A)

            gradient_t = (ddmap_t - self.ddmap_target)
            gradient_t = 0.5 * (gradient_t + gradient_t.T)
            if self.edge_mask is not None:
                gradient_t *= self.edge_mask

            theta_previous = np.copy(self.theta)
            self.theta = self.A + learning_rate * gradient_t
            self.A = self.theta + (t / (t + 3)) * (self.theta - theta_previous)

        if gaussian_noise_variance > 0.0:
            shrink = 1.0 / (1.0 + learning_rate * gaussian_noise_variance)
            self.A = self.A * shrink
            if method == 'GD' and self.theta is not None:
                self.theta = self.theta * shrink

        self.A = np.nan_to_num(self.A)
        self.A = 0.5 * (self.A + self.A.T)
        self._freeze_masked_edges()
        self.A = a2a(self.A, fill_negative=enforce_nonnegative_connectivity_matrix)

        self.loss = self.__compute_loss(ddmap_t)
        eigvals_K = -eigvals_A
        self.entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)

    def __update_parameter_gpu(self, t, learning_rate, lamd=0.0, reg='l2', enforce_nonnegative_connectivity_matrix=False, momentum=0.0, nesterov=False):
        """GPU-accelerated version of __update_parameter for IS method using CuPy."""
        # Compute mean squared distance matrix on GPU
        dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
        # Keep dtype stable (avoid implicit float64 upcast when using cp.pi)
        dd_const = cp.asarray((3.0 * np.pi) / 8.0, dtype=self._A_gpu.dtype)
        ddmap_t_gpu = dd_const * cp.square(dmap_t_gpu)
        
        # Compute ratio and gradient on GPU
        compare_ratio_gpu = ddmap_t_gpu / self._ddmap_target_gpu
        # Keep fhash on GPU to avoid an extra device->host sync each iteration
        fhash = cp.nansum(ddmap_t_gpu) / 2.
        
        # Compute gradient
        gradient_t_gpu = cp.nan_to_num(cp.log(compare_ratio_gpu), posinf=0., neginf=0.) / fhash
        
        # Enforce symmetry
        gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)
        
        # Apply optional edge mask (if set, need to transfer to GPU)
        if self.edge_mask is not None:
            # Use cached mask to avoid per-iteration host->device transfer
            gradient_t_gpu *= self._edge_mask_gpu
        
        # Apply update with optional momentum (standard or Nesterov)
        if momentum > 0.0:
            if t == 0:
                self._velocity_gpu = cp.zeros_like(self._A_gpu)
            self._velocity_gpu = momentum * self._velocity_gpu + gradient_t_gpu
            
            if nesterov:
                self._A_gpu = self._A_gpu + learning_rate * (gradient_t_gpu + momentum * self._velocity_gpu)
            else:
                self._A_gpu = self._A_gpu + learning_rate * self._velocity_gpu
        else:
            self._A_gpu = self._A_gpu + learning_rate * gradient_t_gpu
        
        if lamd > 0.0:
            if reg == 'L2':
                shrink = 1.0 / (1.0 + learning_rate * lamd)
                self._A_gpu = self._A_gpu * shrink
            elif reg == 'L1':
                thresh = learning_rate * lamd
                self._A_gpu = cp.sign(self._A_gpu) * cp.maximum(cp.abs(self._A_gpu) - thresh, 0.0)

        # Clean up NaN values
        self._A_gpu = cp.nan_to_num(self._A_gpu)
        
        # Keep symmetry
        self._A_gpu = 0.5 * (self._A_gpu + self._A_gpu.T)
        
        # Apply edge mask freeze (if set)
        if self.edge_mask is not None:
            # Use cached freeze mask
            self._A_gpu = cp.where(self._freeze_mask_gpu, 0.0, self._A_gpu)
        
        # Rebuild diagonal (a2a operation on GPU)
        # a2a sets diagonal to negative sum of off-diagonal elements in each row
        # and optionally clamps negative off-diagonal values
        if enforce_nonnegative_connectivity_matrix:
            # Clamp off-diagonal to be non-positive (spring constants >= 0)
            diag_vals = cp.diag(self._A_gpu).copy()
            self._A_gpu = cp.minimum(self._A_gpu, 0.0)
            cp.fill_diagonal(self._A_gpu, diag_vals)
        
        # Rebuild diagonal: A_ii = -sum(A_ij for j != i)
        off_diag_sum = cp.sum(self._A_gpu, axis=1) - cp.diag(self._A_gpu)
        cp.fill_diagonal(self._A_gpu, -off_diag_sum)
        
        # Compute loss on GPU (avoid CPU transfer for ddmap)
        loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
        self._loss_gpu = loss_gpu
        
        # Compute entropy on GPU using eigenvalues already on GPU
        # K = -A, so eigenvalues of K = -eigenvalues of A
        # entropy = -sum(log(λ_i)) for positive eigenvalues λ_i of K
        eigvals_K_gpu = -eigvals_A_gpu
        positive_mask = eigvals_K_gpu > 1e-12
        log_terms = cp.where(positive_mask, -cp.log(eigvals_K_gpu), 0.0)
        self._entropy_gpu = cp.sum(log_terms)
        
        # NOTE: We intentionally do NOT sync self.A back to CPU here.
        # Syncing an (n,n) matrix every iteration can dominate runtime for n~1000.
        # We sync only when needed (save_steps) and once at the end of run().

    def __update_parameter_noisy_gpu(self, t, learning_rate, gaussian_noise_variance, enforce_nonnegative_connectivity_matrix=False, momentum=0.0, nesterov=False):
        """GPU-accelerated version of noisy update for IS method using CuPy."""
        dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
        dd_const = cp.asarray((3.0 * np.pi) / 8.0, dtype=self._A_gpu.dtype)
        ddmap_t_gpu = dd_const * cp.square(dmap_t_gpu)

        compare_ratio_gpu = ddmap_t_gpu / self._ddmap_target_gpu
        fhash = cp.nansum(ddmap_t_gpu) / 2.

        gradient_t_gpu = cp.nan_to_num(cp.log(compare_ratio_gpu), posinf=0., neginf=0.) / fhash
        gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)
        if self.edge_mask is not None:
            gradient_t_gpu *= self._edge_mask_gpu

        if momentum > 0.0:
            if t == 0:
                self._velocity_gpu = cp.zeros_like(self._A_gpu)
            self._velocity_gpu = momentum * self._velocity_gpu + gradient_t_gpu
            if nesterov:
                self._A_gpu = self._A_gpu + learning_rate * (gradient_t_gpu + momentum * self._velocity_gpu)
            else:
                self._A_gpu = self._A_gpu + learning_rate * self._velocity_gpu
        else:
            self._A_gpu = self._A_gpu + learning_rate * gradient_t_gpu

        if gaussian_noise_variance > 0.0:
            shrink = 1.0 / (1.0 + learning_rate * gaussian_noise_variance)
            self._A_gpu = self._A_gpu * shrink

        self._A_gpu = cp.nan_to_num(self._A_gpu)
        self._A_gpu = 0.5 * (self._A_gpu + self._A_gpu.T)
        if self.edge_mask is not None:
            self._A_gpu = cp.where(self._freeze_mask_gpu, 0.0, self._A_gpu)

        if enforce_nonnegative_connectivity_matrix:
            diag_vals = cp.diag(self._A_gpu).copy()
            self._A_gpu = cp.minimum(self._A_gpu, 0.0)
            cp.fill_diagonal(self._A_gpu, diag_vals)

        off_diag_sum = cp.sum(self._A_gpu, axis=1) - cp.diag(self._A_gpu)
        cp.fill_diagonal(self._A_gpu, -off_diag_sum)

        loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
        self._loss_gpu = loss_gpu

        eigvals_K_gpu = -eigvals_A_gpu
        positive_mask = eigvals_K_gpu > 1e-12
        log_terms = cp.where(positive_mask, -cp.log(eigvals_K_gpu), 0.0)

        self._entropy_gpu = cp.sum(log_terms)

    def __update_parameter_gpu_gd(self, t, learning_rate, lamd=0.0, reg='l2', enforce_nonnegative_connectivity_matrix=False):
        """GPU-accelerated version of __update_parameter for GD method using CuPy."""
        # Initialize theta on first iteration
        if t == 0:
            self._theta_gpu = self._A_gpu.copy()
        
        # Compute mean squared distance matrix on GPU
        dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
        dd_const = cp.asarray((3.0 * np.pi) / 8.0, dtype=self._A_gpu.dtype)
        ddmap_t_gpu = dd_const * cp.square(dmap_t_gpu)
        
        # Compute gradient for GD
        gradient_t_gpu = ddmap_t_gpu - self._ddmap_target_gpu
        
        # Enforce symmetry
        gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)
        
        # Apply optional edge mask
        if self.edge_mask is not None:
            gradient_t_gpu *= self._edge_mask_gpu
        
        # Nesterov-like update (built into GD method)
        theta_previous_gpu = self._theta_gpu.copy()
        self._theta_gpu = self._A_gpu + learning_rate * gradient_t_gpu
        
        # Momentum update: A = theta + (t/(t+3)) * (theta - theta_previous)
        momentum_rate = t / (t + 3)
        self._A_gpu = self._theta_gpu + momentum_rate * (self._theta_gpu - theta_previous_gpu)
        
        if lamd > 0.0:
            if reg == 'L2':
                shrink = 1.0 / (1.0 + learning_rate * lamd)
                self._A_gpu = self._A_gpu * shrink
                self._theta_gpu = self._theta_gpu * shrink
            elif reg == 'L1':
                thresh = learning_rate * lamd
                self._A_gpu = cp.sign(self._A_gpu) * cp.maximum(cp.abs(self._A_gpu) - thresh, 0.0)
                self._theta_gpu = cp.sign(self._theta_gpu) * cp.maximum(cp.abs(self._theta_gpu) - thresh, 0.0)

        # Clean up NaN values
        self._A_gpu = cp.nan_to_num(self._A_gpu)
        
        # Keep symmetry
        self._A_gpu = 0.5 * (self._A_gpu + self._A_gpu.T)
        
        # Apply edge mask freeze (if set)
        if self.edge_mask is not None:
            self._A_gpu = cp.where(self._freeze_mask_gpu, 0.0, self._A_gpu)
        
        # Rebuild diagonal (a2a operation on GPU)
        if enforce_nonnegative_connectivity_matrix:
            diag_vals = cp.diag(self._A_gpu).copy()
            self._A_gpu = cp.minimum(self._A_gpu, 0.0)
            cp.fill_diagonal(self._A_gpu, diag_vals)
        
        # Rebuild diagonal: A_ii = -sum(A_ij for j != i)
        off_diag_sum = cp.sum(self._A_gpu, axis=1) - cp.diag(self._A_gpu)
        cp.fill_diagonal(self._A_gpu, -off_diag_sum)
        
        # Compute loss on GPU (avoid CPU transfer for ddmap)
        loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
        self._loss_gpu = loss_gpu
        
        # Compute entropy on GPU using eigenvalues already on GPU
        # K = -A, so eigenvalues of K = -eigenvalues of A
        # entropy = -sum(log(λ_i)) for positive eigenvalues λ_i of K
        eigvals_K_gpu = -eigvals_A_gpu
        positive_mask = eigvals_K_gpu > 1e-12
        log_terms = cp.where(positive_mask, -cp.log(eigvals_K_gpu), 0.0)
        self._entropy_gpu = cp.sum(log_terms)
        
        # NOTE: We intentionally do NOT sync self.A back to CPU here.
        # Sync only when needed (save_steps) and once at the end of run().

    def __update_parameter_noisy_gpu_gd(self, t, learning_rate, gaussian_noise_variance, enforce_nonnegative_connectivity_matrix=False):
        """GPU-accelerated version of noisy update for GD method using CuPy."""
        if t == 0:
            self._theta_gpu = self._A_gpu.copy()

        dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
        dd_const = cp.asarray((3.0 * np.pi) / 8.0, dtype=self._A_gpu.dtype)
        ddmap_t_gpu = dd_const * cp.square(dmap_t_gpu)

        gradient_t_gpu = ddmap_t_gpu - self._ddmap_target_gpu
        gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)
        if self.edge_mask is not None:
            gradient_t_gpu *= self._edge_mask_gpu

        theta_previous_gpu = self._theta_gpu.copy()
        self._theta_gpu = self._A_gpu + learning_rate * gradient_t_gpu
        momentum_rate = t / (t + 3)
        self._A_gpu = self._theta_gpu + momentum_rate * (self._theta_gpu - theta_previous_gpu)

        if gaussian_noise_variance > 0.0:
            shrink = 1.0 / (1.0 + learning_rate * gaussian_noise_variance)
            self._A_gpu = self._A_gpu * shrink
            self._theta_gpu = self._theta_gpu * shrink

        self._A_gpu = cp.nan_to_num(self._A_gpu)
        self._A_gpu = 0.5 * (self._A_gpu + self._A_gpu.T)
        if self.edge_mask is not None:
            self._A_gpu = cp.where(self._freeze_mask_gpu, 0.0, self._A_gpu)

        if enforce_nonnegative_connectivity_matrix:
            diag_vals = cp.diag(self._A_gpu).copy()
            self._A_gpu = cp.minimum(self._A_gpu, 0.0)
            cp.fill_diagonal(self._A_gpu, diag_vals)

        off_diag_sum = cp.sum(self._A_gpu, axis=1) - cp.diag(self._A_gpu)
        cp.fill_diagonal(self._A_gpu, -off_diag_sum)

        loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
        self._loss_gpu = loss_gpu

        eigvals_K_gpu = -eigvals_A_gpu
        positive_mask = eigvals_K_gpu > 1e-12
        log_terms = cp.where(positive_mask, -cp.log(eigvals_K_gpu), 0.0)
        self._entropy_gpu = cp.sum(log_terms)

    def run(self, epoch, general_method='optimization', save_steps=None, output_prefix=None, **kwargs):
        """
        Main function to run the optimization
        
        Parameters
        ----------
        epoch : int
            Number of iterations
        general_method : str
            'optimization' or 'direct'
        save_steps : list of int, optional
            Iteration steps at which to capture the connectivity matrix.
            Matrices are returned in the 5th return value (dict step -> matrix).
            If output_prefix is set, files are also saved as
            '{output_prefix}_connectivity_matrix_iter{step}.txt'
        output_prefix : str, optional
            Prefix for output files when saving connectivity matrix at save_steps
        **kwargs
            Additional arguments passed to __update_parameter:
            - learning_rate : float
                Learning rate for optimization
            - lamd : float
                Regularization weight
            - reg : str
                Regularization type ('L1' or 'L2')
            - method : str
                Optimization method ('IS' or 'GD')
            - enforce_nonnegative_connectivity_matrix : bool
                Enforce non-negative spring constants
            - momentum : float, optional
                Momentum coefficient for IS method (default: 0.0, i.e., no momentum).
                RECOMMENDED: Use 0.95 with nesterov=True for fastest convergence.
                Use 0.9 for conservative settings. Only applies when method='IS'.
            - nesterov : bool, optional
                If True and momentum > 0, use Nesterov Accelerated Gradient (NAG).
                NAG enables higher momentum (0.95) without divergence.
                RECOMMENDED: Use with momentum=0.95 for ~50% faster convergence.

        Returns
        -------
        loss_array : list of float
        entropy_array : list of float
        dmap_maxent : np.ndarray
        A : np.ndarray
            Final connectivity matrix.
        connectivity_at_steps : dict of int -> np.ndarray
            Connectivity matrix at each step in save_steps (empty dict if save_steps not set).
        """

        console = Console()

        # For GPU optimization runs we keep loss/entropy histories on device and copy once at the end,
        # avoiding per-iteration device->host synchronization.
        loss_array = []
        entropy_array = []
        use_gpu_hist = self.use_gpu and general_method == 'optimization'
        if use_gpu_hist:
            loss_hist_gpu = cp.empty(epoch, dtype=self._gpu_dtype)
            entropy_hist_gpu = cp.empty(epoch, dtype=self._gpu_dtype)
        
        # Convert save_steps to a set for fast lookup and validate
        save_steps_set = None
        connectivity_at_steps = {}  # step -> matrix (for library use and/or file save)
        if save_steps is not None:
            save_steps_list = list(save_steps)
            # Filter out invalid steps (must be positive and <= epoch)
            valid_steps = [s for s in save_steps_list if 1 <= s <= epoch]
            if len(valid_steps) < len(save_steps_list):
                invalid = [s for s in save_steps_list if s < 1 or s > epoch]
                console.print(f"[yellow]Warning: Some save_steps are out of range and will be ignored: {invalid}[/yellow]")
            save_steps_set = set(valid_steps)
            if len(save_steps_set) > 0:
                console.print(f"[green]Will save connectivity matrix at iterations: {sorted(save_steps_set)}[/green]")

        if general_method == 'optimization':
            with trange(epoch, desc="Performing optimization", unit="iteration") as pbar:
                for t in pbar:
                    self.__update_parameter(t, **kwargs)
                    if self.use_gpu:
                        # Record without syncing to host
                        loss_hist_gpu[t] = self._loss_gpu
                        entropy_hist_gpu[t] = self._entropy_gpu

                        # Only sync occasionally for display (tqdm set_postfix wants Python scalars)
                        if (t == 0) or ((t + 1) % self.gpu_display_every == 0) or (t + 1 == epoch):
                            self.loss = float(self._loss_gpu)
                            self.entropy = float(self._entropy_gpu)
                            pbar.set_postfix(loss=self.loss, entropy=self.entropy)
                    else:
                        # CPU path
                        pbar.set_postfix(loss=self.loss, entropy=self.entropy if self.entropy is not None else np.nan)
                        loss_array.append(self.loss)
                        entropy_array.append(self.entropy if self.entropy is not None else np.nan)
                    
                    # Save connectivity matrix at specified iteration steps
                    # Note: t is 0-indexed, so we check for t+1 to match iteration number
                    if save_steps_set is not None and (t + 1) in save_steps_set:
                        # For GPU mode, sync to CPU only when needed
                        if self.use_gpu:
                            self.A = cp.asnumpy(self._A_gpu)
                        step = t + 1
                        connectivity_at_steps[step] = np.copy(self.A)
                        if output_prefix is not None:
                            filename = '{}_connectivity_matrix_iter{}.txt'.format(output_prefix, step)
                            np.savetxt(filename, self.A)
                            console.print(f"[green]Saved connectivity matrix at iteration {step} to {filename}[/green]")
        elif general_method == 'direct':
            if not checkEMD(self.ddmap_target):
                raise ValueError(
                    'The distance matrix is a not valid Euclidean distance matrix. Direct inversion method is not applicable. Please use optimization method such as Iterative scaling or gradient descent')
            self.A = ddmap2a_direct(self.ddmap_target)
            # Compute loss for direct inversion
            ddmap_t = ((3. * np.pi) / 8.) * np.power(a2dmap_theory(self.A, force_positive_definite=True), 2.)
            loss_array.append(self.__compute_loss(ddmap_t))
            # Compute entropy for direct inversion
            eigvals_K = np.linalg.eigvalsh(-self.A)
            entropy_array.append(compute_entropy_from_A(self.A, eigvals=eigvals_K))

        # Finalize loss/entropy arrays for GPU runs
        if use_gpu_hist:
            loss_array = cp.asnumpy(loss_hist_gpu).tolist()
            entropy_array = cp.asnumpy(entropy_hist_gpu).tolist()
            # Ensure latest scalars are synced for callers that read model.loss/entropy
            self.loss = float(self._loss_gpu) if self._loss_gpu is not None else None
            self.entropy = float(self._entropy_gpu) if self._entropy_gpu is not None else None

        # Ensure self.A is up-to-date on CPU before returning / computing dmap
        if self.use_gpu and general_method == 'optimization':
            self.A = cp.asnumpy(self._A_gpu)
        dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)

        return loss_array, entropy_array, dmap_maxent, self.A, connectivity_at_steps

    def run_noisy(self, epoch, gaussian_noise_variance, general_method='optimization', save_steps=None, output_prefix=None, **kwargs):
        """
        Run optimization with independent Gaussian noise on constraints.

        Parameters
        ----------
        epoch : int
            Number of iterations
        gaussian_noise_variance : float
            Noise variance for constraints (L2 shrink strength).
        general_method : str
            Only 'optimization' is supported for the noisy solver.
        save_steps : list of int, optional
            Iteration steps at which to save the connectivity matrix.
        output_prefix : str, optional
            Prefix for output files (required if save_steps is provided)
        **kwargs
            Additional arguments passed to __update_parameter_noisy:
            - learning_rate : float
                Learning rate for optimization
            - method : str
                Optimization method ('IS' or 'GD')
            - enforce_nonnegative_connectivity_matrix : bool
                Enforce non-negative spring constants
            - momentum : float, optional
                Momentum coefficient for IS method (default: 0.0).
            - nesterov : bool, optional
                If True and momentum > 0, use Nesterov Accelerated Gradient (NAG).
        """
        if general_method != 'optimization':
            raise ValueError("run_noisy supports only general_method='optimization'")

        console = Console()

        if self.use_gpu:
            loss_hist_gpu = cp.empty(epoch, dtype=self._gpu_dtype)
            entropy_hist_gpu = cp.empty(epoch, dtype=self._gpu_dtype)
        else:
            loss_array = []
            entropy_array = []

        save_steps_set = None
        if save_steps is not None:
            save_steps_list = list(save_steps)
            valid_steps = [s for s in save_steps_list if 1 <= s <= epoch]
            if len(valid_steps) < len(save_steps_list):
                invalid = [s for s in save_steps_list if s < 1 or s > epoch]
                console.print(f"[yellow]Warning: Some save_steps are out of range and will be ignored: {invalid}[/yellow]")
            save_steps_set = set(valid_steps)
            if output_prefix is None:
                raise ValueError("output_prefix must be provided when save_steps is specified")
            if len(save_steps_set) > 0:
                console.print(f"[green]Will save connectivity matrix at iterations: {sorted(save_steps_set)}[/green]")

        with trange(epoch, desc="Performing noisy optimization", unit="iteration") as pbar:
            for t in pbar:
                self.__update_parameter_noisy(t, gaussian_noise_variance=gaussian_noise_variance, **kwargs)
                if self.use_gpu:
                    loss_hist_gpu[t] = self._loss_gpu
                    entropy_hist_gpu[t] = self._entropy_gpu

                    if (t == 0) or ((t + 1) % self.gpu_display_every == 0) or (t + 1 == epoch):
                        self.loss = float(self._loss_gpu)
                        self.entropy = float(self._entropy_gpu)
                        pbar.set_postfix(loss=self.loss, entropy=self.entropy)
                else:
                    pbar.set_postfix(loss=self.loss, entropy=self.entropy if self.entropy is not None else np.nan)
                    loss_array.append(self.loss)
                    entropy_array.append(self.entropy if self.entropy is not None else np.nan)

                if save_steps_set is not None and (t + 1) in save_steps_set:
                    if self.use_gpu:
                        self.A = cp.asnumpy(self._A_gpu)
                    filename = '{}_connectivity_matrix_iter{}.txt'.format(output_prefix, t + 1)
                    np.savetxt(filename, self.A)
                    console.print(f"[green]Saved connectivity matrix at iteration {t + 1} to {filename}[/green]")

        if self.use_gpu:
            loss_array = cp.asnumpy(loss_hist_gpu).tolist()
            entropy_array = cp.asnumpy(entropy_hist_gpu).tolist()
            self.loss = float(self._loss_gpu) if self._loss_gpu is not None else None
            self.entropy = float(self._entropy_gpu) if self._entropy_gpu is not None else None

        if self.use_gpu:
            self.A = cp.asnumpy(self._A_gpu)
        dmap_maxent = a2dmap_theory(self.A, force_positive_definite=True)

        return loss_array, entropy_array, dmap_maxent, self.A

class Dynamics:
    def __init__(self, input, M=None, k=None, model=None):
        if isinstance(input, int) and M is None and k is not None:
            if not isinstance(k, int) and not isinstance(k, float):
                sys.stdout.write('Spring constant should be a number')
                sys.exit(0)
            if model != 'rouse' and model is not None:
                sys.stdout.write("Please specify model to be 'rouse'")
                sys.exit(0)

            self.A = construct_connectivity_matrix_rouse(input, k)
            self.eigvalue, self.eigvector = np.linalg.eigh(self.A)
            self.N = input
        elif isinstance(input, int) and M is not None and k is not None:
            if not isinstance(k, int) and not isinstance(k, float):
                sys.stdout.write('Spring constant should be a number')
                sys.exit(0)
            if not isinstance(M, int):
                sys.stdout.write('Number of random cross links need to be an integer')
                sys.exit(0)
            if model != 'random':
                sys.stdout.write("Please specify model to be 'random'")
                sys.exit(0)
            self.A = construct_connectivity_matrix_random(input, M, k)
            self.eigvalue, self.eigvector = np.linalg.eigh(self.A)
            self.N = input
        elif isinstance(input, np.ndarray) and M is None and k is None:
            if len(input.shape) !=2 or input.shape[0] != input.shape[1]:
                sys.stdout.write('The connectivity matrix should be a square matrix')
                sys.exit(0)
            if not np.allclose(input, input.T):
                sys.stdout.write('The connectivity matrix should be a symmetrix real matrix')
                sys.exit(0)

            self.A = input
            self.eigvalue, self.eigvector = np.linalg.eigh(self.A)
            self.N = input.shape[0]

    def generateXYZ(self, force_positive_definite = False):
        self.xyz = a2xyz_sample(self.A, force_positive_definite = force_positive_definite)[0]
        self.modes = self.eigvector.T @ self.xyz

    def initialize(self, dt, zeta, beta):
        if not isinstance(dt, int) and not isinstance(dt, float):
            sys.stdout.write('Time step should be a number')
            sys.exit(0)
        if not isinstance(zeta, int) and not isinstance(zeta, float):
            sys.stdout.write('Friction coefficient step should be a number')
            sys.exit(0)
        if not isinstance(beta, int) and not isinstance(beta, float):
            sys.stdout.write('Temperature step should be a number')
            sys.exit(0)
        elif beta <= 0.0:
            sys.stdout.write('Temperature should be positive')
            sys.exit(0)

        self.zeta = zeta
        self.beta = beta
        self.dt = dt

    def updateModes(self, method='euler-maruyama', force_projection=None, update_zero_modes=True):
        try:
            self.zeta
            self.beta
            self.dt
        except AttributeError:
            sys.stdout.write('Please run initialize() first')
            sys.exit(0)

        # Pass force_projection directly to Ornstein_Uhlenbeck_update
        # The function will handle the conversion to long-term mean internally
        if force_projection is None:
            force_projection = 0.0
            
        self.modes = Ornstein_Uhlenbeck_update(self.modes, self.dt, - self.eigvalue, self.zeta, self.beta, 
                                             force_projection=force_projection, method=method, 
                                             update_zero_modes=update_zero_modes)

    def updateXYZ(self):
        self.xyz = self.eigvector @ self.modes

    def run(self, T, update=1, every=1, initial_conformation=None, method='euler-maruyama', update_zero_modes=True):
        """
        T: number of timesteps
        update: update x,y,z positions every this many timesteps
        every: save x,y,z positions to the trajectory every this many timesteps
        initial_conformation: initial conformation of the simulation
        update_zero_modes: if True, COM will diffuse; if False, COM is fixed at initial position
        """
        if not isinstance(T, int):
            sys.stdout.write('Number of steps should be an integer')
            sys.exit(0)

        if initial_conformation is None:
            try:
                self.xyz
                self.modes
            except AttributeError:
                self.generateXYZ()
        else:
            if initial_conformation.shape[0] != self.N:
                sys.stdout.write('Number of particles is not correct')
                sys.exit(0)
            if initial_conformation.shape[1] != 3:
                sys.stdout.write('The dimension should be three')
                sys.exit(0)
            self.xyz = initial_conformation
            self.modes = self.eigvector.T @ self.xyz

        self.traj = []
        for t in tqdm(range(T)):
            if t % update == 0:
                self.updateXYZ()
            if t % every == 0:
                self.updateXYZ()
                self.traj.append(self.xyz.copy())
            self.updateModes(method=method, update_zero_modes=update_zero_modes)

        self.traj = np.array(self.traj)

    def run_with_force(self, T, force_loci, force_amplitude, force_direction, force_duration=None, 
                      update=1, every=1, initial_conformation=None, method='euler-maruyama', update_zero_modes=True):
        """
        Run dynamics simulation with applied forces following the correct derivation.
        
        The equation of motion is: ξ dR/dt = A*R + B + f
        In eigenmodes: ξ dX/dt = Λ*X + V^T*B + f̃
        
        Parameters
        ----------
        T : int
            Total number of timesteps
        force_loci : list of int
            Indices of loci where force is applied
        force_amplitude : float
            Magnitude of the force
        force_direction : (3,) array_like
            Direction vector of the force (will be normalized)
        force_duration : int, optional
            Number of timesteps to apply force. If None, force is applied for entire simulation.
        update : int, optional
            Update x,y,z positions every this many timesteps
        every : int, optional
            Save x,y,z positions to the trajectory every this many timesteps
        initial_conformation : (N, 3) array_like, optional
            Initial conformation of the simulation
        method : str, optional
            Integration method ('euler-maruyama' or 'exact')
        update_zero_modes : bool, optional
            If True, COM will diffuse; if False, COM is fixed at initial position
        """
        if not isinstance(T, int):
            sys.stdout.write('Number of steps should be an integer')
            sys.exit(0)

        if initial_conformation is None:
            try:
                self.xyz
                self.modes
            except AttributeError:
                self.generateXYZ()
        else:
            if initial_conformation.shape[0] != self.N:
                sys.stdout.write('Number of particles is not correct')
                sys.exit(0)
            if initial_conformation.shape[1] != 3:
                sys.stdout.write('The dimension should be three')
                sys.exit(0)
            self.xyz = initial_conformation
            self.modes = self.eigvector.T @ self.xyz

        # Normalize force direction
        force_direction = np.array(force_direction)
        force_direction = force_direction / np.linalg.norm(force_direction)
        
        # Create force vector B (following correct derivation)
        B = np.zeros((self.N, 3))
        for locus in force_loci:
            B[locus] = force_amplitude * force_direction

        # Project force onto eigenmodes: V^T * B
        # V is self.eigvector, so V^T is self.eigvector.T
        # This transforms forces from real space to eigenmode space
        force_projection = self.eigvector.T @ B  # Shape: (n_modes, 3)

        # Allow COM to move naturally due to force
        # No COM constraint applied
        
        self.traj = []
        for t in tqdm(range(T)):
            if t % update == 0:
                self.updateXYZ()
            if t % every == 0:
                self.updateXYZ()
                self.traj.append(self.xyz.copy())
            
            # Apply force only during specified duration
            if force_duration is None or t < force_duration:
                self.updateModes(method=method, force_projection=force_projection, update_zero_modes=update_zero_modes)
            else:
                self.updateModes(method=method, update_zero_modes=update_zero_modes)

        self.traj = np.array(self.traj)

    def reset(self):
        self.generateXYZ()

    def run_breakable_bond(self, T, cutoff_distance, update=1, every=1, initial_conformation=None, 
                          method='euler-maruyama', update_zero_modes=True):
        """
        Run dynamics simulation with breakable bonds based on distance cutoff.
        
        IMPORTANT: Nearest-neighbor bonds (i, i+1) representing the polymer backbone are 
        NEVER broken, regardless of distance. Only non-local interactions (i, j where j > i+1) 
        are subject to the distance cutoff.
        
        At each time step:
        1. Check distances between all pairs of loci
        2. For nearest-neighbor bonds (i, i+1): always keep original connectivity value
        3. For non-local bonds (i, j where j > i+1):
           - If distance < cutoff: bond is present (use original connectivity value)
           - If distance >= cutoff: bond is broken (set connectivity to zero)
        4. Update connectivity matrix and recompute eigenvalues/eigenvectors
        5. Update normal modes using the new connectivity matrix
        
        Parameters
        ----------
        T : int
            Total number of timesteps
        cutoff_distance : float
            Distance threshold for breaking non-local bonds. Only applies to bonds 
            between non-adjacent loci (i, j where j > i+1).
        update : int, optional
            Update x,y,z positions every this many timesteps
        every : int, optional
            Save x,y,z positions to the trajectory every this many timesteps
        initial_conformation : (N, 3) array_like, optional
            Initial conformation of the simulation
        method : str, optional
            Integration method ('euler-maruyama' or 'exact'). Default is 'euler-maruyama'.
            NOTE: 'euler-maruyama' is STRONGLY RECOMMENDED for breakable bond simulations,
            as broken bonds can create positive eigenvalues (unstable modes). The euler-maruyama
            method handles these gracefully with controlled growth over small timesteps, while
            the 'exact' method may produce numerical issues.
        update_zero_modes : bool, optional
            If True, COM will diffuse; if False, COM is fixed at initial position
        
        Notes
        -----
        The input connectivity matrix is preserved and not modified. All calculations 
        use an internal copy (A_internal). Nearest-neighbor bonds are never broken 
        to maintain polymer chain connectivity.
        
        **Handling Positive Eigenvalues:**
        When bonds break, the connectivity matrix may develop positive eigenvalues (unstable modes).
        This causes exponential growth in those modes. However, the feedback mechanism (bonds
        break as distances increase) combined with small timesteps prevents catastrophic blowup.
        The system remains stable because: (1) growth is incremental with small dt, and (2) as
        distances grow, more bonds break, reducing the instability.
        """
        if not isinstance(T, int):
            sys.stdout.write('Number of steps should be an integer')
            sys.exit(0)
        if not isinstance(cutoff_distance, (int, float)) or cutoff_distance <= 0:
            sys.stdout.write('Cutoff distance should be a positive number')
            sys.exit(0)

        # Store original connectivity matrix for reference (keep input unchanged)
        self.A_original = np.copy(self.A)
        
        # Create internal connectivity matrix for breakable bond calculations
        self.A_internal = np.copy(self.A)
        
        if initial_conformation is None:
            try:
                self.xyz
                self.modes
            except AttributeError:
                self.generateXYZ()
        else:
            if initial_conformation.shape[0] != self.N:
                sys.stdout.write('Number of particles is not correct')
                sys.exit(0)
            if initial_conformation.shape[1] != 3:
                sys.stdout.write('The dimension should be three')
                sys.exit(0)
            self.xyz = initial_conformation
            self.modes = self.eigvector.T @ self.xyz

        self.traj = []
        for t in tqdm(range(T)):
            if t % update == 0:
                self.updateXYZ()
            if t % every == 0:
                self.updateXYZ()
                self.traj.append(self.xyz.copy())
            
            # Update internal connectivity matrix based on current distances
            self._update_connectivity_from_distances(cutoff_distance)
            
            # Update normal modes (eigvalue and eigvector already updated in _update_connectivity_from_distances)
            self.updateModes(method=method, update_zero_modes=update_zero_modes)

        self.traj = np.array(self.traj)

    def _update_connectivity_from_distances(self, cutoff_distance):
        """
        Update internal connectivity matrix based on current pairwise distances.
        
        Nearest-neighbor bonds (i, i+1) are never broken - they represent the polymer backbone.
        Only non-local interactions are subject to the distance cutoff.
        
        Parameters
        ----------
        cutoff_distance : float
            Distance threshold for bond breaking (only applies to non-nearest-neighbor bonds)
        """
        # Calculate pairwise distances
        distances = self._calculate_pairwise_distances()
        
        # Update internal connectivity matrix
        # Only loop over non-local bonds (j > i+1) since nearest-neighbor bonds never break
        for i in range(self.N):
            for j in range(i+2, self.N):
                # Non-local bonds: apply distance cutoff
                if distances[i, j] < cutoff_distance:
                    # Bond is present - use original connectivity value
                    self.A_internal[i, j] = self.A_original[i, j]
                    self.A_internal[j, i] = self.A_original[j, i]
                else:
                    # Bond is broken - set to zero
                    self.A_internal[i, j] = 0.0
                    self.A_internal[j, i] = 0.0
        
        # Update diagonal elements to maintain Laplacian property using a2a function
        self.A_internal = a2a(self.A_internal, fill_negative=False)
        
        # Recompute eigenvalues and eigenvectors using internal matrix
        self.eigvalue, self.eigvector = np.linalg.eigh(self.A_internal)
        
        # Update modes to maintain consistency
        self.modes = self.eigvector.T @ self.xyz

    def _calculate_pairwise_distances(self):
        """
        Calculate pairwise distances between all loci using scipy.
        
        Returns
        -------
        distances : (N, N) ndarray
            Matrix of pairwise distances. distances[i,j] is the distance between loci i and j.
        """
        # Use pdist+squareform for better performance on larger systems
        # pdist computes condensed distance matrix, squareform converts to full symmetric matrix
        distances = scipy.spatial.distance.squareform(
            scipy.spatial.distance.pdist(self.xyz, metric='euclidean')
        )
        return distances

    def restore_original_connectivity(self):
        """
        Restore the original connectivity matrix and recompute eigenvalues/eigenvectors.
        This is useful after running breakable bond simulations to return to the original state.
        """
        if hasattr(self, 'A_original'):
            self.A = np.copy(self.A_original)
            self.eigvalue, self.eigvector = np.linalg.eigh(self.A)
            # Update modes to maintain consistency
            self.modes = self.eigvector.T @ self.xyz
        else:
            sys.stdout.write('No original connectivity matrix found. Cannot restore.')
