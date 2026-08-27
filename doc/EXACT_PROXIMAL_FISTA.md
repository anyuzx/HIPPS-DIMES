# Exact proximal-gradient and monotone FISTA for noisy HIPPS-DIMES

This branch adds an independent parameter-space solver for the Gaussian
noise-aware maximum-entropy dual.

For unique locus pairs `a=(i,j)`, define the natural parameter

\[
\theta_{ij}=-A_{ij}/2,
\]

where `A` is the HIPPS-DIMES connectivity matrix.  With independent Gaussian
errors of variance \(v_{ij}\) in the observed mean-squared distances, the dual
objective is

\[
L(\theta)=\log Z(\theta)-\sum_{i<j}\theta_{ij}D_{ij}^{\rm obs}
+\frac12\sum_{i<j}v_{ij}\theta_{ij}^{2}.
\]

The exact smooth gradient is

\[
\nabla f(\theta)_{ij}=D_{ij}^{\rm fit}(\theta)-D_{ij}^{\rm obs},
\]

and the exact quadratic proximal map is

\[
\operatorname{prox}_{\eta g}(x)_{ij}
=\frac{x_{ij}}{1+\eta v_{ij}}.
\]

Therefore one exact proximal-gradient update is

\[
\theta_{ij}^{+}=
\frac{\theta_{ij}-\eta
\left(D_{ij}^{\rm fit}-D_{ij}^{\rm obs}\right)}
{1+\eta v_{ij}}.
\]

The implementation evaluates the gradient at a FISTA extrapolated point,
uses smooth-term backtracking, rejects non-normalizable connectivity matrices,
and restarts whenever the accelerated proposal would increase the full
objective.  Every accepted iterate is spectrally valid and the accepted
objective sequence is monotone.

The exact KKT condition is

\[
D_{ij}^{\rm fit}-D_{ij}^{\rm obs}+v_{ij}\theta_{ij}=0,
\]

or, in connectivity notation,

\[
D_{ij}^{\rm fit}-D_{ij}^{\rm obs}=\frac{v_{ij}}{2}A_{ij}.
\]

## Python API

```python
from hipps_dimes import fit_gaussian_noise_dual_fista

fitted_ddmap, gram, connectivity, info = fit_gaussian_noise_dual_fista(
    observed_ddmap,
    noise_variance=0.1,
    max_iterations=2000,
    relative_tolerance=1e-5,
)
```

A symmetric variance matrix is accepted for heteroskedastic errors.  A shared
relative standard deviation can be supplied instead:

```python
fitted_ddmap, gram, connectivity, info = fit_gaussian_noise_dual_fista(
    observed_ddmap,
    relative_noise_std=0.2,
)
```

Set `accelerated=False` to run exact proximal gradient without FISTA.

## Scope

This is an opt-in research solver.  It does not replace the existing `IS`,
`GD`, or `DI` paths on this branch.  The implementation is intentionally
separate so that its convergence and runtime can be benchmarked against PDHG
and covariance-space methods before changing the production default.
