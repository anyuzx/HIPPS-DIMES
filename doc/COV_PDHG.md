# Variance-whitened PDHG and direct-Gram FISTA COV optimization

`fit_gaussian_noise_covariance_pdhg` is the sole PDHG implementation and
`fit_gaussian_noise_covariance_hybrid` is the public PDHG-to-FISTA workflow for
the Gaussian noise-aware COV objective,

$$
F(B)=-\frac32\log\det{}'B
+\frac12\sum_a\frac{[\mathcal D_a(B)-d_a]^2}{v_a},
$$

where `B` is the centered three-dimensional Gram matrix, `a=(i,j)` indexes an
observed locus pair, and

$$
\mathcal D_{ij}(B)=B_{ii}+B_{jj}-2B_{ij}.
$$

The translational mode remains zero, and every internal mode remains strictly
positive. The former scalar-step PDHG implementation has been replaced rather
than retained as a fallback.

## Observation graph requirement

Finite, positive off-diagonal squared-distance entries define an undirected
observation graph. COV requires this graph to connect all retained loci. If it
has multiple components, their relative translations are unconstrained and the
log-determinant term makes the stated objective unbounded below; a small
finite-iterate residual would not certify a finite optimum. The validator scans
the dense observed-pair mask in $O(N^2)$ time and uses only $O(N)$ additional
working memory.

Partially missing matrices are supported when
`ignore_missing_data=True` and their observation graph remains connected. With
the default `ignore_missing_data=False`, the high-level API interpolates every
remaining missing pair first, so the completed map becomes the target. A fully
missing locus must be repaired or removed explicitly before either pair policy
is applied. Repaired and interpolated pairs are genuine target constraints:
they enter both the data term and the variance constructed by the chosen COV
noise model.

## Variance whitening

Define

$$
\mathcal A=V^{-1/2}\mathcal D,\qquad b=V^{-1/2}d.
$$

The objective becomes

$$
F(B)=-\frac32\log\det{}'B+\frac12\|\mathcal A(B)-b\|_2^2.
$$

For the whitened dual variable `u`, PDHG applies

$$
u^{k+1}=\frac{u^k+\sigma[\mathcal A(\bar B^k)-b]}{1+\sigma},
$$

followed by the analytic proximal update

$$
B^{k+1}=\operatorname{prox}_{\tau(-3\log\det'/2)}
[B^k-\tau\mathcal A^*u^{k+1}].
$$

Whitening gives every observed-pair dual coordinate unit quadratic curvature.
For `relative_noise_std=c`, the standardized target is the constant `1/c`.
The low-level PDHG default initializes the dual with the whitened data
residual; it does not run a separate multi-candidate `auto` scoring pass.

For step-size selection, define the edge-space normal operator

$$
G=\mathcal A\mathcal A^*.
$$

It is symmetric, positive semidefinite, and entrywise nonnegative. Starting
from a strictly positive edge vector $x$, the Collatz bounds give

$$
\min_a\frac{(Gx)_a}{x_a}
\leq \lambda_{\max}(G)=\|\mathcal A\|^2
\leq \max_a\frac{(Gx)_a}{x_a}.
$$

The implementation repeatedly applies $G$ without materializing it and uses
the maximum ratio as a certified upper bound, even if the configured bound-gap
tolerance is not reached. A further safety factor determines the PDHG step
product `tau*sigma`; adaptation changes only the ratio `tau/sigma`. This
certificate avoids relying on a heuristic power-iteration start that could
miss the dominant eigenspace.

## Inverse-free runtime KKT residual

The primal proximal optimality equation gives

$$
\frac32(B^{k+1})^+
=\mathcal A^*u^{k+1}+\frac{B^{k+1}-B^k}{\tau}.
$$

Consequently, the eliminated stationarity residual can be evaluated during
PDHG without reconstructing the inverse Gram matrix:

$$
g_{\mathrm{elim}}^{k+1}
=-\frac{B^{k+1}-B^k}{\tau}
+\mathcal A^*[\mathcal A(B^{k+1})-b-u^{k+1}].
$$

Step adaptation compares the two terms in this common primal matrix space.
The inverse and physical connectivity are reconstructed only for initial dual
setup, requested checkpoints, the returned model, and the independent final
certificate.

Because PDHG objective values need not be monotone, `return_best=True` returns
the iterate with the smallest runtime eliminated KKT residual. The returned
physical Gram matrix is then certified independently by recomputing

$$
R(B)=\mathcal D^*\!\left[
\frac{\mathcal D(B)-D^{\mathrm{obs}}}{v}
\right]-\frac32B^+.
$$

Convergence requires

$$
\|R(B)\|_F\leq
\mathrm{atol}+\mathrm{rtol}\max\left(
\|\nabla F_{\mathrm{data}}\|_F,
\left\|\frac32B^+\right\|_F
\right).
$$

The default relative tolerance is `1e-5`. `max_iterations` is only a budget;
exhausting it does not imply convergence. Check `info["converged"]`,
`info["independent_kkt_converged"]`, and
`info["relative_eliminated_kkt_residual"]`.

`info["iterations"]` is the number of updates executed, whereas
`info["returned_iteration"]` identifies the selected model. The history column
`is_returned_iterate` marks that row. Reported top-level `objective`, `loss`,
and `entropy` are all evaluated at that same returned iterate.

## Direct-Gram monotone FISTA refinement

At the hybrid handoff, FISTA receives PDHG's centered physical Gram matrix
directly. The handoff does not invert through a connectivity matrix and does
not repeat the scalar initialization calibration. For the smooth data term

$$
h(B)=\frac12\|V^{-1/2}[\mathcal D(B)-d]\|_2^2,
$$

the gradient is evaluated with the same matrix-free weighted distance
operator and adjoint used by PDHG. The nonsmooth-side proximal map is the
closed-form positive-definite update for $-3\log\det'(B)/2$, applied only to
the centered internal modes. In a hybrid run, FISTA reuses PDHG's certified
weighted-operator bound and physical distance scale. It uses backtracking when
needed and restarts acceleration whenever an extrapolated step would violate
monotonic objective descent.

The direct FISTA function is an internal low-level implementation. It is not
exported from `HippsDimes` and is not a standalone CLI optimizer. This keeps
the supported user contract narrow: FISTA is entered only after a certified
PDHG handoff.

## Public use

The supported optimizer names are:

- `covariance_optimizer="hybrid"` (default), which runs variance-whitened
  PDHG to the independently certified handoff tolerance and then direct-Gram
  monotone FISTA;
- `covariance_optimizer="pdhg"`, which continues standalone PDHG to the final
  tolerance or update budget.

The default handoff tolerance is `1e-2`; the default final relative KKT
tolerance is `1e-5`. Both phases share the single `iteration` budget. A
standalone PDHG run can be selected for controlled solver comparisons:

```python
import hipps_dimes as HD

results = HD.run_optimization(
    input_matrix=target_ddmap,
    input_type="ddmap",
    method="COV",
    gaussian_noise_relative_std=0.10,
    covariance_optimizer="pdhg",
    covariance_relative_tolerance=1e-5,
    iteration=20000,
    use_gpu=True,
)

info = results["covariance_optimization"]
print(info["converged"])
print(info["relative_eliminated_kkt_residual"])
print(info["returned_iteration"])
```

The low-level canonical entry point is
`fit_gaussian_noise_covariance_pdhg(...)`. There are no separate `whitened`,
`preconditioned`, or legacy PDHG public functions.

## Scaling and memory

Both phases store dense centered matrices, so matrix storage is $O(N^2)$.
Observed pairs are represented by index and value vectors, and neither phase
materializes a dense Hessian or an $(N-1)^2$ linear-system workspace. FISTA
does reconstruct the inverse Gram matrix and perform an internal dense
eigendecomposition for every accepted update, so its arithmetic remains
cubic in $N$ even though its memory growth is quadratic. PDHG avoids that
inverse reconstruction during ordinary iterations and can therefore remain
the cheaper phase far from the solution.

CPU and GPU backends execute the same float64 algorithm. A requested GPU run
requires CuPy and an accessible CUDA device; it fails clearly instead of
silently falling back to CPU.

## Diagnostics and tuning

The PDHG history records the objective components, eliminated KKT residual,
common-space primal and dual components, step sizes, step-ratio adaptations,
Gram conditioning, and connectivity norm. Additional metadata include:

- `weighted_operator_norm` for the edge-space Collatz certificate, including
  lower and certified upper bounds and their relative gap;
- `inverse_reconstruction_count` and
  `inverse_reconstructed_each_iteration`;
- `pdhg["variance_whitened"]`;
- `pdhg["inverse_free_runtime_kkt"]`;
- `pdhg["weighted_residual_balancing"]`.

The hybrid history uses one global accepted-update axis and labels each row
with `phase` and `phase_iteration`. FISTA rows additionally report the
accepted step size, backtracking count, momentum coefficient, restart flag,
and objective decrease. Hybrid metadata include:

- `phase_iterations` and `phase_wall_seconds` for `pdhg` and `fista`;
- `initial_step_size`, `final_step_size`, `backtracking_reductions`, and
  `restart_count`;
- `handoff["physical_gram_used_directly"] = True` and
  `handoff["scalar_recalibration"] = None`.

The main low-level tuning parameters are `step_safety`, `initial_dual_step`,
`step_ratio`, residual-adaptation controls, and the operator-bound iteration and
gap-tolerance controls. The high-level API uses the validated defaults.

Progress callbacks distinguish `covariance_operator_norm` and
`covariance_optimization`, with `phase` set to `pdhg` or `fista` in hybrid
runs. Operator-norm events report their own relative bound gap and do not
contain an objective value.

When invoked through the CLI, a nonconverged COV run writes the requested
partial artifacts and then exits with status 1. Library calls return those same
artifacts, emit a `RuntimeWarning`, and expose the false convergence flag for
programmatic handling.

## Scaling evidence

These measurements are descriptive rather than runtime contracts:

- On the trusted real GM12878 `N=400` fixture, the `1e-2` hybrid reached the
  final independent KKT `9.972e-6` in a 27.046-second median wall time, versus
  34.364 seconds for continued PDHG. Its objective relative error against the
  stored trusted solution was `9.230e-10`; relative differences in squared
  distance, distance, Gram, and connectivity were respectively `9.818e-6`,
  `3.923e-6`, `2.382e-5`, and `1.609e-5`. Peak CuPy pool reservation was
  52.50 MiB for both paths.
- On the GM12878 whole-chromosome-3 `N=991` fixture, the public default hybrid
  reached independent KKT `9.99924e-6` in 8,784 total updates and 477.70
  seconds. PDHG reached the `1e-2` handoff after 1,046 updates and 54.79
  seconds; FISTA used 7,738 updates and 421.81 seconds. The minimum internal
  Gram eigenvalue was `3.186e-3`, maximum absolute Gram and connectivity row
  sums were `6.82e-13` and `1.06e-13`, and the retained CuPy pool reservation
  was 301.23 MiB.

Historical replacement evidence: the removed Newton implementation used 3,214
PDHG updates plus 20 refinement steps, took 474.7 seconds, returned KKT
`3.14e-4`, and reported globalization failure on the same `N=991` case. This
value is retained only to explain the replacement; Newton is no longer a
supported optimizer or code path.

The test suite protects the objective, independent KKT certificate, stored
solution agreement, SPD/centering invariants, and GPU workflow on the real
`N=400` fixture. The larger measurement remains benchmark evidence and is not
run in ordinary CI.
