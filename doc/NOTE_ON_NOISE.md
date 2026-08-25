# HIPPS-DIMES with noisy constraints (notes)

## Setup
HIPPS-DIMES (MaxEnt / exponential family dual):
\[
L(\theta) = \log Z(\theta) - \theta^\top a
\]
- Constraints: \(a\) is the vector of all pairwise targets \(a_{ij}=\langle r_{ij}^2\rangle\) (usually \(i<j\)).
- Parameters: \(\theta_{ij}\) are the Lagrange multipliers for each constraint.
  - With the HIPPS-DIMES Hamiltonian \(H=\frac12\sum_{i<j}K_{ij}r_{ij}^2\),
    \(\theta_{ij}=-K_{ij}/2\), with connectivity diagonals determined by
    row-sum zero.

---

## Why log-ratios appear in iterative scaling
- In classic iterative scaling (GIS/IIS) for MaxEnt models, updates often involve terms like:
  \[
  \Delta\theta_k \propto \log\frac{a_k}{\mathbb E_\theta[f_k]}
  \]
- This comes from the log-linear structure plus feature conditions (e.g., nonnegativity and certain boundedness/constant-sum assumptions) that enable monotone improvement.
- In practice, “IS-like” updates using \(\log(\text{current}/\text{target})\) are common even when those strict conditions don’t exactly hold; they often behave like a stable multiplicative correction or a preconditioned step.

**Key nuance:** log-ratio updates emphasize *relative* error:
\[
\log\frac{\mathbb E_\theta[f]}{\hat a} \approx \frac{\mathbb E_\theta[f]-\hat a}{\hat a}\quad(\text{small errors})
\]

---

## Incorporating uncertainty in targets \(a\) naturally

### Measurement model (independent Gaussian noise)
Assume observed targets:
\[
\hat a_\alpha = a_\alpha + \varepsilon_\alpha,\qquad \varepsilon_\alpha\sim\mathcal N(0,\sigma^2)
\]
with no correlation between constraints (pair-independent).

### Soft-constraint MaxEnt (primal)
Instead of enforcing \(\mathbb E_p[f]=\hat a\) exactly, use:
\[
\max_p\; H(p)\;-\;\frac{1}{2\sigma^2}\|\mathbb E_p[f]-\hat a\|_2^2
\]
Entropy remains the driving principle; uncertainty determines how “soft” constraints are.

### Equivalent regularized dual
This yields (up to an additive constant):
\[
\boxed{
L_{\text{noise}}(\theta)
= \log Z(\theta) - \theta^\top \hat a + \frac{\sigma^2}{2}\|\theta\|_2^2
}
\]

**Gradient:**
\[
\nabla_\theta L_{\text{noise}}(\theta)=\mathbb E_\theta[f]-\hat a + \sigma^2\theta
\]
**Stationary condition:**
\[
\mathbb E_\theta[f] = \hat a - \sigma^2\theta
\]
So the model does **not** match \(\hat a\) exactly; it trades off fit vs plausible parameter magnitude.

---

## What the L2 penalty looks like for pairwise MSD constraints
Index \(\alpha\) as unique pairs \((i<j)\). Then:
\[
\boxed{
\text{Penalty} = \frac{\sigma^2}{2}\sum_{i<j}\theta_{ij}^2
}
\]
If \(\sigma^2=0.005\):
\[
\boxed{
\text{Penalty} = \frac{0.005}{2}\sum_{i<j}\theta_{ij}^2
}
\]

### Relation to connectivity \(K_{ij}\)
For HIPPS-DIMES, \(K_{ij}=-2\theta_{ij}\) for \(i\neq j\), so:
\[
\frac{\sigma^2}{2}\sum_{i<j}\theta_{ij}^2
=
\frac{\sigma^2}{8}\sum_{i<j}K_{ij}^2
\]
(assuming you regularize only the independent off-diagonals).

**Avoid common pitfalls:**
- Don’t double count pairs:
  \[
  \sum_{i\neq j}K_{ij}^2 = 2\sum_{i<j}K_{ij}^2
  \]
- Don’t blindly use \(\sum_{i,j}K_{ij}^2\): it includes diagonals, which are not independent if you enforce row-sum zero.

---

## Iterative scaling + L2: what’s “correct”
- Adding \(\frac{\sigma^2}{2}\|\theta\|^2\) is fully principled from the Gaussian-noise-on-targets assumption.
- Whether you implement it as an extra “shrinkage” term inside an IS-like update depends on your update rule.

### Proximal shrink inside a generic gradient method
If your IS-like update proposes:
\[
\theta \leftarrow \theta + \eta\;\Delta(\theta)
\]
then apply exact L2 shrink after the step:
\[
\boxed{
\theta \leftarrow \frac{\theta + \eta\;\Delta(\theta)}{1+\eta\sigma^2}
}
\]
This is the proximal step for the L2 term relative to that proposed update.
However, composing it with an ad-hoc log-ratio direction does **not** make the
overall iteration an optimizer of the Gaussian dual: its fixed point generally
differs from \(\nabla L_{\mathrm{noise}}=0\).

HIPPS-DIMES therefore retains the historical noisy IS path only for backward
compatibility. The calibrated alternatives are COV, FISTA, CHOL, and CIS.
FISTA applies the analytic proximal map of the negative log-determinant in the
covariance cone with backtracking and monotone acceleration; it minimizes the
same objective as COV. CIS performs exact cyclic coordinate minimization of the
connectivity dual. All four calibrated solvers check the equivalent physical
stationarity condition
\[
D^{\mathrm{fit}}_{ij}-D^{\mathrm{obs}}_{ij}
=\frac{\sigma^2_{ij}}{2}K_{ij}.
\]

### Optional signed sparsity in addition to Gaussian noise

FISTA and CIS can retain the variance-derived quadratic term and add the
independent connectivity prior
\[
\lambda_1\sum_{i<j}|K_{ij}|.
\]
This permits both signs but favors exact zero couplings. The common dual is
\[
-\frac32\log\det P
+\frac12\sum_{i<j}K_{ij}\hat D_{ij}
+\frac18\sum_{i<j}\sigma_{ij}^2K_{ij}^2
+\lambda_1\sum_{i<j}|K_{ij}|.
\]
Its continuous KKT residual is
\[
\operatorname{soft}\!\left(
D^{\mathrm{fit}}_{ij}-D^{\mathrm{obs}}_{ij},\,2\lambda_1
\right)
-\frac{\sigma_{ij}^2}{2}K_{ij}.
\]
CIS minimizes the negative, zero, and positive branches of every coordinate
exactly. FISTA uses the equivalent covariance-space squared dead-zone loss,
whose gradient is the soft-thresholded residual divided by the variance. Both
reduce exactly to their original Gaussian implementations at \(\lambda_1=0\).
For positive \(\lambda_1\), FISTA uses the common physical KKT residual above
as its convergence certificate. This avoids continuing after the
covariance-gradient diagnostic reaches its floating-point floor. The
zero-coefficient path retains its original full covariance-gradient stopping
rule exactly.
In the CLI and `run_optimization`, select this model with `--method FISTA` or
`--method CIS` together with `--lamd LAMBDA --reg L1`.

---

## Practical interpretation
- Using \(\log(\text{current}/\text{target})\) is not “wrong”; it’s common in MaxEnt iterative scaling–type methods.
- Adding L2 regularization is also natural and corresponds directly to uncertainty in the constraints.
- Additive Gaussian noise on MSD targets suggests linear-residual least-squares is statistically matched; log-ratio emphasizes relative errors (often desirable when MSD spans orders of magnitude).

---
