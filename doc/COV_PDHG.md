# Variance-whitened PDHG and hybrid COV optimization

`fit_gaussian_noise_covariance_pdhg` is the sole PDHG implementation for the
Gaussian noise-aware COV objective,

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

The implementation estimates

$$
\|\mathcal A\|^2
=\lambda_{\max}(\mathcal D^*V^{-1}\mathcal D)
$$

with matrix-free power iteration. The estimate includes a safety factor, and
PDHG preserves the resulting safe product `tau*sigma` while adapting only the
ratio `tau/sigma`.

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

## Public use

The established optimizer names are unchanged:

- `covariance_optimizer="pdhg"` selects standalone variance-whitened PDHG;
- `covariance_optimizer="hybrid"` runs the same PDHG to the handoff tolerance
  and then invokes Newton-CG;
- `covariance_optimizer="newton"` selects standalone Newton-CG.

`hybrid` remains the default because Newton refinement is effective for the
validated smaller systems. For large systems, select PDHG explicitly:

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

## Large-system Newton warning

For an optimized matrix with `N > 2000`, selecting `hybrid` or `newton` emits
a runtime warning before expensive work begins. The warning recommends
standalone PDHG but does not reject the request or change the selected
optimizer.

This threshold is operational guidance, not a mathematical restriction on
Newton. The existing Newton stage forms a dense `(N-1)^2` linear system
implicitly and uses a diagonal preconditioner. Its equations, CG solver, and
line search are unchanged.

## Diagnostics and tuning

The PDHG history records the objective components, eliminated KKT residual,
common-space primal and dual components, step sizes, step-ratio adaptations,
Gram conditioning, and connectivity norm. Additional metadata include:

- `weighted_operator_norm` for power-iteration setup diagnostics;
- `inverse_reconstruction_count` and
  `inverse_reconstructed_each_iteration`;
- `pdhg["variance_whitened"]`;
- `pdhg["inverse_free_runtime_kkt"]`;
- `pdhg["weighted_residual_balancing"]`.

The main low-level tuning parameters are `step_safety`, `initial_dual_step`,
`step_ratio`, residual-adaptation controls, and the operator-norm power-iteration
controls. The high-level API uses the validated defaults.

Progress callbacks distinguish `covariance_operator_norm`, the PDHG
`covariance_optimization` phase, `covariance_preconditioner`, and the Newton
`covariance_optimization` phase. Operator-norm events report their own relative
residual and do not contain an objective value.

## Scaling evidence

These measurements motivate the replacement and warning but are not runtime
contracts:

- On the real GM12878 `N=400` fixture, variance-whitened hybrid reached the
  `1e-3` handoff after 1,108 PDHG updates and converged after three Newton
  updates in 24.36 seconds. The former scalar-step hybrid required 3,881 PDHG
  and two Newton updates in 73.79 seconds. Objectives matched, and the relative
  Gram-matrix difference was `1.764e-9`.
- On the chromosome 3 `N=3801` target with 5,567,171 observed pairs, whitened
  PDHG reached relative KKT values `1.044e-1` at 1,000 updates, `4.861e-3` at
  5,000, and `1.053e-3` at 9,500. These values document progress, not final
  `1e-5` convergence.
- From the same `N=3801` PDHG checkpoint, one-update Newton diagnostics with
  CG caps from 50 through 500 ended with relative linear residuals from `4.58`
  to `9.77`, far above the requested `1e-4`. The accepted descent steps lowered
  the objective but did not improve the KKT certificate.

The test suite protects the objective and independent KKT certificate on the
real `N=400` fixture. The `N=3801` measurements remain benchmark evidence and
are not run in ordinary CI.
