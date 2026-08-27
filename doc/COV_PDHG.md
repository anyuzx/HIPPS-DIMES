# PDHG solver for Gaussian noise-aware COV

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
the iterate with the smallest maximum relative KKT residual when
`return_best=True`.

## Use

```python
import HippsDimes

fitted_ddmap, gram, connectivity, info = (
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

The current Newton-CG solver is unchanged and remains the default COV solver.
PDHG is intentionally opt-in until it is benchmarked against the N=400 cases
and the converged nearest-EDM reference.

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

The history records objective components, both KKT residuals, step sizes,
step-ratio adaptations, Gram conditioning, and connectivity norm.
