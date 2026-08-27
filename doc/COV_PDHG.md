# PDHG and hybrid solvers for Gaussian noise-aware COV

`fit_gaussian_noise_covariance_pdhg` solves the same convex objective as the
noise-aware COV Newton solver,

\[
F(B)=-\frac32\log\det{}'B
+\frac12\sum_{a}\frac{(\mathcal D_a B-d_a)^2}{v_a},
\]

where `B` is the centered three-dimensional Gram matrix, `a=(i,j)` indexes the
observed locus pairs, and

\[
\mathcal D_{ij}(B)=B_{ii}+B_{jj}-2B_{ij}.
\]

The translational mode of `B` remains zero and every internal mode remains
strictly positive.

## Primal-dual splitting

Write

\[
g(B)=-\frac32\log\det{}'B,
\qquad
h(z)=\frac12\sum_a\frac{(z_a-d_a)^2}{v_a}.
\]

The conjugate of the pairwise Gaussian term is

\[
h^*(y)=d^Ty+\frac12\sum_a v_a y_a^2.
\]

The Chambolle--Pock/PDHG iteration is

\[
y_a^{k+1}=\frac{y_a^k+\sigma
[\mathcal D_a(\bar B^k)-d_a]}{1+\sigma v_a},
\]

\[
Y^k=B^k-\tau\mathcal D^*y^{k+1},
\]

followed by the analytic log-determinant proximal step. If an internal
eigendecomposition of `Y` is

\[
Y=U\operatorname{diag}(\lambda_p)U^T,
\]

then

\[
B^{k+1}=U\operatorname{diag}\left[
\frac{\lambda_p+\sqrt{\lambda_p^2+6\tau}}{2}
\right]U^T.
\]

Every updated internal eigenvalue is positive, even if `Y` is indefinite. The
extrapolated variable is

\[
\bar B^{k+1}=B^{k+1}+\theta(B^{k+1}-B^k).
\]

For the complete unique-pair distance operator on centered symmetric matrices,

\[
\|\mathcal D\|^2=2N.
\]

For any observed-pair subset the norm is no larger, so the implementation uses

\[
\tau\sigma=\frac{s^2}{2N},\qquad 0<s<1.
\]

Only the ratio `tau/sigma` is adapted; their product remains fixed.

## Optimality certificates

The solver monitors both KKT equations:

\[
\mathcal D^*y-\frac32B^+=0,
\]

and

\[
\mathcal D(B)-d-v\odot y=0.
\]

The second equation, combined with

\[
A=-3B^+,
\]

is the HIPPS-DIMES pair stationarity relation

\[
D_{ij}^{\rm fit}-D_{ij}^{\rm obs}
=\frac{v_{ij}}{2}A_{ij}.
\]

Because PDHG objective values need not be monotone, the implementation returns
the iterate with the smallest dual-eliminated relative KKT residual when
`return_best=True`. It then independently recomputes the pseudoinverse and
certificate from that returned Gram matrix. `info["converged"]` is true only
when the internal stopping conditions and this final certificate pass.

## Use

The high-level COV interface defaults to the hybrid optimizer. It runs PDHG to
an independently certified relative KKT residual of `1e-3`, then initializes
the existing centered Newton-CG solver from that feasible result. The phases
share one maximum-iteration budget.

```python
results = HippsDimes.run_optimization(
    input_matrix=target_ddmap,
    input_type="ddmap",
    method="COV",
    gaussian_noise_relative_std=0.10,
    covariance_optimizer="hybrid",
    covariance_relative_tolerance=1e-5,
    covariance_handoff_relative_tolerance=1e-3,
    iteration=10000,
    use_gpu=True,
)
```

The low-level hybrid and standalone PDHG solvers remain available directly:

```python
import HippsDimes

fitted_ddmap, gram, connectivity, info = (
    HippsDimes.fit_gaussian_noise_covariance_hybrid(
        target_ddmap,
        relative_noise_std=0.10,
        initialization="rouse",
        max_iterations=10000,
        handoff_relative_tolerance=1e-3,
        use_gpu=True,
    )
)

pdhg_only = (
    HippsDimes.fit_gaussian_noise_covariance_pdhg(
        target_ddmap,
        relative_noise_std=0.10,
        initialization="rouse",
        max_iterations=2000,
        use_gpu=True,
    )
)
```

For homoskedastic errors, pass `noise_variance=<positive scalar>` instead of
`relative_noise_std`.

The standalone solvers remain available through
`covariance_optimizer="pdhg"` and `covariance_optimizer="newton"`. Newton-CG
is not deleted or replaced; the hybrid path reuses it as its local phase.

## Main tuning parameters

- `step_safety`: controls the fixed product `tau*sigma`; default `0.99`.
- `step_ratio`: optional explicit `tau/sigma`. If omitted, a dimensionless,
  variance-aware value is selected automatically.
- `adaptive_steps`: balance the two KKT residuals while preserving the safe
  step-size product.
- `adaptation_interval`, `adaptation_threshold`, `adaptation_factor`: control
  residual balancing.
- `dual_initialization`: `auto`, `zero`, `residual`, or `connectivity`.
- `theta`: extrapolation coefficient in `[0,1]`; default `1`.
- `relative_tolerance`: production relative KKT tolerance; default `1e-5`.
  Explicitly pass `1e-8` for strict small-system validation.
- `handoff_relative_tolerance`: hybrid PDHG-to-Newton relative KKT threshold;
  default `1e-3`.

The history records objective components, primal, dual, and dual-eliminated
KKT residuals, step sizes, step-ratio adaptations, Gram conditioning, and
connectivity norm.

## Real-data regression

The test suite includes the processed experimental GM12878 Hi-C contact map
for `chr1:31,000,000-41,000,000` at 25 kb resolution (`N=400`) and a converged
reference Gram matrix. The ordinary CPU suite independently reconstructs its
target squared-distance observations and verifies the COV objective and KKT
certificate of the full-size reference solution.

When CuPy and a CUDA GPU are available, the `real_data` test also runs the
complete default hybrid contact-map workflow from Rouse initialization:

```bash
pytest -q -m real_data tests/test_covariance_pdhg.py
```

CPU-only environments skip this approximately 70-second end-to-end solve but
still execute the fast full-size objective/KKT regression.
