# Globalized COV optimization

The `feature/cov-globalization` branch keeps the calibrated Gaussian COV
objective unchanged,

\[
F(B)=-\frac{3}{2}\log\det B+
\frac12\sum_{i<j}\frac{[D_{ij}(B)-D_{ij}^{\rm obs}]^2}{v_{ij}},
\qquad B\succ0,
\]

but changes the numerical path used to reach its unique minimizer.

## Default optimization path

The installed COV solver now uses

1. the existing observed-pair Rouse initialization;
2. the existing exact scalar minimization along the initial Gram-matrix ray;
3. variance continuation when the calibrated start is still far from the final
   optimum;
4. monotone proximal-gradient warm-up using the analytic proximal operator of
   `-3/2 logdet(B)`;
5. inexact Newton-CG with a gradient-dependent forcing tolerance;
6. Newton equations in whitened covariance coordinates;
7. an exact positive-definite scalar line search, with Armijo backtracking only
   as a fallback.

For easy fits, continuation is skipped automatically and the solver proceeds
directly to the final objective.

## Anchored covariance coordinates

The default `coordinate_parameterization="anchored"` fixes the final locus as
the coordinate origin and optimizes the `(N-1) x (N-1)` Gram matrix of
`r_i-r_anchor`. In this gauge,

\[
D_{i,a}=B_{ii},\qquad
D_{ij}=B_{ii}+B_{jj}-2B_{ij}.
\]

Consequently, distance evaluation, the data gradient, and every data-Hessian
action require only elementwise matrix operations and are `O(N^2)`, rather than
dense transformations through an orthonormal centered basis. The returned Gram
matrix is converted back to the usual centered representation, and the returned
connectivity matrix remains symmetric, row-sum zero, and negative semidefinite.

The former centered parameterization remains available through
`coordinate_parameterization="centered"` for numerical comparisons.

## Symmetric-coordinate preconditioner

The preconditioner is defined in the orthonormal symmetric `svec` basis. For an
off-diagonal covariance coordinate, the basis matrix is
`(E_ij+E_ji)/sqrt(2)`, so the data-Hessian diagonal is `2 w_ij`, not `w_ij`.
The entropy-Hessian diagonal is evaluated as

\[
H_{ii,ii}^{\rm ent}=\frac32(B^{-1}_{ii})^2,
\]

\[
H_{ij,ij}^{\rm ent}=\frac32
\left(B^{-1}_{ii}B^{-1}_{jj}+(B^{-1}_{ij})^2\right),
\qquad i<j.
\]

For anchored coordinates the exact data-Hessian diagonal is constructed in
`O(N^2)` time.

## New diagnostics

`results['covariance_optimization']['history']` now includes

- `phase` (`proximal_warmup` or `newton`);
- `continuation_factor`;
- `cg_forcing_tolerance` and `cg_converged`;
- `line_search_backtracks`;
- `full_step_accepted`;
- `line_search_method`;
- `feasible_step_bound`;
- `fallback_direction`.

The top-level COV diagnostics also report the continuation schedule,
parameterization, preconditioner type, and globalization settings.

## Programmatic controls

The public function accepts the following additional keyword arguments:

```python
fit_gaussian_noise_covariance(
    ...,
    continuation_factors=(0.1, 1.0),
    continuation_intermediate_tolerance=1e-3,
    continuation_intermediate_newton_iterations=5,
    proximal_warmup_iterations=8,
    proximal_switch_relative_gradient=0.25,
    continuation_activation_relative_gradient=1.0,
    cg_forcing_max=0.5,
    exact_line_search=True,
    line_search_max_step=1.0,
    use_whitened_newton=True,
    coordinate_parameterization="anchored",
)
```

These controls alter convergence behavior only. The final `factor=1` objective
and its stationary solution are unchanged.
