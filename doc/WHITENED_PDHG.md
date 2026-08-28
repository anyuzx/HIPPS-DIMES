# Variance-whitened PDHG for large heteroskedastic COV problems

The implementation in `hipps_dimes/covariance_pdhg_whitened.py` solves the same noise-aware covariance objective as the existing PDHG and Newton-CG solvers:

\[
F(B)=-\frac32\log\det{}'B+\frac12\sum_a\frac{[D_a(B)-d_a]^2}{v_a}.
\]

It introduces three numerical changes intended for large systems with strongly heterogeneous variances.

## B. Variance whitening

Define

\[
\mathcal A=V^{-1/2}\mathcal D,\qquad b=V^{-1/2}d.
\]

Then

\[
F(B)=-\frac32\log\det{}'B+\frac12\|\mathcal A B-b\|_2^2.
\]

For the whitened dual variable `u`, the updates are

\[
u^{k+1}=\frac{u^k+\sigma[\mathcal A\bar B^k-b]}{1+\sigma},
\]

\[
B^{k+1}=\operatorname{prox}_{\tau(-3\log\det'/2)}
[B^k-\tau\mathcal A^*u^{k+1}].
\]

The dual quadratic curvature is now one for every pair. The code estimates

\[
\|\mathcal A\|^2=\lambda_{\max}(\mathcal D^*V^{-1}\mathcal D)
\]

by matrix-free power iteration and chooses a safe product `tau*sigma`.

For `relative_noise_std=c`, the standardized target is constant:

\[
\frac{d_{ij}}{\sqrt{v_{ij}}}=\frac1c.
\]

## C. Inverse-free runtime KKT

The primal proximal optimality equation gives

\[
\frac32(B^{k+1})^{-1}
=\mathcal A^*u^{k+1}+\frac{B^{k+1}-B^k}{\tau}.
\]

Therefore the eliminated gradient can be evaluated without constructing `B^{-1}`:

\[
g_{\rm elim}^{k+1}
=-\frac{B^{k+1}-B^k}{\tau}
+\mathcal A^*[\mathcal A B^{k+1}-b-u^{k+1}].
\]

The inverse and connectivity are reconstructed only for initial dual setup, requested checkpoints, the returned model, and the independent final certificate.

## D. Weighted residual balancing

Let

\[
r_p=-\frac{B^{k+1}-B^k}{\tau},
\]

\[
r_d=\mathcal A B^{k+1}-b-u^{k+1}.
\]

Step adaptation compares `||r_p||` with `||A* r_d||`, which live in the same primal matrix space and add directly to the eliminated gradient. The solver changes `tau/sigma` while preserving the safe product `tau*sigma` and resets extrapolation after each adaptation.

## Public functions

```python
fit_gaussian_noise_covariance_pdhg_whitened(...)
fit_gaussian_noise_covariance_preconditioned_pdhg(...)  # alias
fit_gaussian_noise_covariance_hybrid_whitened(...)
```

Example:

```python
fitted, gram, connectivity, info = (
    HippsDimes.fit_gaussian_noise_covariance_pdhg_whitened(
        target_ddmap,
        relative_noise_std=0.1,
        use_gpu=True,
        max_iterations=20000,
        relative_tolerance=1e-3,
    )
)
```

The original scalar-step PDHG and hybrid functions remain unchanged for controlled comparisons.

Useful diagnostics are available under `info["history"]`:

- `relative_eliminated_kkt_residual`
- `primal_component_norm`
- `dual_component_primal_norm`
- `component_balance_ratio`
- `primal_step`, `dual_step`, and `step_ratio`
- `step_adapted`
- `gram_condition_number`

The operator-norm setup diagnostics are stored in `info["weighted_operator_norm"]`, and `info["inverse_reconstruction_count"]` verifies that the inverse was not reconstructed on every PDHG iteration.
