# HIPPS-DIMES Acceleration and GPU Technical Note

This note documents the convergence and performance changes implemented in `HippsDimes.py`, including consolidated eigendecompositions, momentum/Nesterov acceleration for Iterative Scaling (IS), and GPU acceleration for both IS and Gradient Descent (GD).

## 1) Theory Recap (minimum needed context)

The optimization targets the maximum-entropy Gaussian model whose connectivity matrix `A` (Laplacian-like, negative semidefinite) reproduces a target mean-squared distance map:

```text
ddmap[i,j] = <||x_i - x_j||^2>
```

The model-predicted mean-squared distances are computed from the mean distance map:

```text
ddmap = (3*pi/8) * (dmap)^2
```

For IS, the update is driven by the log-ratio between current and target ddmaps:

```text
gradient = log(ddmap_current / ddmap_target) / fhash
```

For GD, the update is driven by the residual:

```text
gradient = ddmap_current - ddmap_target
```

## 2) Consolidated eigendecomposition

**Motivation:** `a2dmap_theory()` already performs an eigendecomposition of `A`. Entropy computation used to recompute an eigendecomposition of `K = -A`. We now reuse the eigenvalues from `a2dmap_theory()` to avoid an extra `O(n^3)` call.

**Implementation:**

```python
def a2dmap_theory(A, force_positive_definite=False, return_eigenvalues=False):
    # ...
    eigvalue, eigvector = scipy.linalg.eigh(A)
    # ...
    dmap = 2.0 * np.sqrt(2.0 / np.pi) * sigma
    
    if return_eigenvalues:
        return dmap, eigvalue
    return dmap
```

```python
# In __update_parameter (CPU path)
dmap_t, eigvals_A = a2dmap_theory(self.A, force_positive_definite=True, return_eigenvalues=True)
ddmap_t = ((3. * np.pi) / 8.) * np.power(dmap_t, 2.)

# Entropy uses eigenvalues of K = -A, so eigvals_K = -eigvals_A
eigvals_K = -eigvals_A
self.entropy = compute_entropy_from_A(self.A, eigvals=eigvals_K)
```

## 3) Momentum and Nesterov acceleration for IS

**Theory:** IS is a fixed-point iteration. Adding Polyak momentum accelerates convergence by mixing current gradient with past updates. Nesterov acceleration uses a look-ahead correction to stabilize higher momentum values.

**Implementation:**

```python
# In __update_parameter (IS path)
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
```

**Recommended settings:**

- `momentum=0.95, nesterov=True` (fastest stable setting)
- `momentum=0.9` (conservative)

## 4) GPU acceleration (CuPy)

GPU acceleration is optional and used only when `use_gpu=True` **and** CuPy is available. The GPU path avoids CPU-GPU transfers inside the iteration by keeping intermediate arrays on the GPU (including loss and entropy).

**CuPy availability detection:**

```python
try:
    import cupy as cp
    cp.cuda.runtime.getDeviceCount()
    _CUPY_AVAILABLE = True
    # ...
except ImportError:
    pass
```

**GPU version of `a2dmap_theory`:**

```python
def _a2dmap_theory_gpu(A_gpu, force_positive_definite=False, return_eigenvalues=False):
    eigvalue, eigvector = cp.linalg.eigh(A_gpu)
    temp = -1.0 / eigvalue
    temp = cp.where(cp.isinf(temp), 0.0, temp)
    # ...
    dmap = 2.0 * cp.sqrt(2.0 / cp.pi) * sigma
    if return_eigenvalues:
        return dmap, eigvalue
    return dmap
```

**GPU Iterative Scaling (IS):**

```python
dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
ddmap_t_gpu = ((3. * cp.pi) / 8.) * cp.power(dmap_t_gpu, 2.)
compare_ratio_gpu = ddmap_t_gpu / self._ddmap_target_gpu
fhash = float(cp.nansum(ddmap_t_gpu) / 2.)

gradient_t_gpu = cp.nan_to_num(cp.log(compare_ratio_gpu), posinf=0., neginf=0.) / fhash
gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)

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

# Loss/entropy computed on GPU
loss_gpu = cp.nanmean(cp.power((ddmap_t_gpu - self._ddmap_target_gpu) / self._ddmap_target_gpu, 2.)) ** 0.5
eigvals_K_gpu = -eigvals_A_gpu
log_terms = cp.where(eigvals_K_gpu > 1e-12, -cp.log(eigvals_K_gpu), 0.0)
self.entropy = float(cp.sum(log_terms))
```

**GPU Gradient Descent (GD):**

```python
dmap_t_gpu, eigvals_A_gpu = _a2dmap_theory_gpu(self._A_gpu, force_positive_definite=True, return_eigenvalues=True)
ddmap_t_gpu = ((3. * cp.pi) / 8.) * cp.power(dmap_t_gpu, 2.)
gradient_t_gpu = ddmap_t_gpu - self._ddmap_target_gpu
gradient_t_gpu = 0.5 * (gradient_t_gpu + gradient_t_gpu.T)

theta_previous_gpu = self._theta_gpu.copy()
self._theta_gpu = self._A_gpu + learning_rate * gradient_t_gpu
momentum_rate = t / (t + 3)
self._A_gpu = self._theta_gpu + momentum_rate * (self._theta_gpu - theta_previous_gpu)
```

## 5) Entropy computation consistency (CPU/GPU)

Entropy is computed from eigenvalues of `K = -A` as:

```text
H = -sum_i log(λ_i),   λ_i > 0
```

CPU and GPU now use consistent formulas:

```python
positive_mask = eigvals > zero_tol
log_terms = -np.log(eigvals[positive_mask])
```

## 6) CLI and programmatic usage

**CLI options added:**

```python
@click.option('--momentum', type=click.FloatRange(0, 1), default=0.0, show_default=True, help='Momentum coefficient for IS method...')
@click.option('--nesterov', is_flag=True, default=False, show_default=True, help='Use Nesterov Accelerated Gradient (NAG)...')
@click.option('--use-gpu', is_flag=True, default=False, show_default=True, help='Use GPU acceleration via CuPy...')
```

**Programmatic usage:**

```python
results = run_optimization(
    input_matrix=cmap,
    input_type='cmap',
    method='IS',
    iteration=5000,
    learning_rate=10.0,
    momentum=0.95,
    nesterov=True,
    use_gpu=True,
)
```

## 7) Practical recommendations

- **IS + Nesterov (momentum=0.95)** is the fastest stable convergence setting observed.
- **GPU acceleration** provides substantial speedup for large matrices (n >= 200), with the biggest gains when computation stays entirely on GPU.
- **Consolidated eigendecomposition** removes redundant `O(n^3)` work per iteration and is always beneficial.
