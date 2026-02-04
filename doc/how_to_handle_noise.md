# HIPPS-DIMES with noisy constraints (notes)

## Setup
HIPPS-DIMES (MaxEnt / exponential family dual):
\[
L(\theta) = \log Z(\theta) - \theta^\top a
\]
- Constraints: \(a\) is the vector of all pairwise targets \(a_{ij}=\langle r_{ij}^2\rangle\) (usually \(i<j\)).
- Parameters: \(\theta_{ij}\) are the Lagrange multipliers for each constraint.
  - In a “spring/Laplacian” parameterization, off-diagonal connectivity relates as \(K_{ij} = -\theta_{ij}\) for \(i\neq j\), with diagonals determined by row-sum zero.

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
If your free parameters are off-diagonal \(K_{ij}\) with \(K_{ij}=-\theta_{ij}\) for \(i\neq j\), then:
\[
\frac{\sigma^2}{2}\sum_{i<j}\theta_{ij}^2
=
\frac{\sigma^2}{2}\sum_{i<j}K_{ij}^2
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

### Proximal (clean) way to add L2 shrink
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
This isolates the L2 effect from any ad-hoc normalization used in \(\Delta(\theta)\).

---

## Practical interpretation
- Using \(\log(\text{current}/\text{target})\) is not “wrong”; it’s common in MaxEnt iterative scaling–type methods.
- Adding L2 regularization is also natural and corresponds directly to uncertainty in the constraints.
- Additive Gaussian noise on MSD targets suggests linear-residual least-squares is statistically matched; log-ratio emphasizes relative errors (often desirable when MSD spans orders of magnitude).

---

