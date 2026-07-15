# Pawula Positivity Diagnostic Plan

This plan defines a non-production toy diagnostic for the generalized NLO
Fokker-Planck normal form used in the current derivation notes.

## Scope

The diagnostic is intentionally one-dimensional and synthetic. It does not use
physical KLM coordinate kernels, does not modify production evolution code, and
does not train score or Hessian-score models.

The purpose is to expose a basic mathematical caveat: once finite third-order
derivative terms are present, positivity is not automatic in the ordinary
Markov semigroup sense.

## Background

LO JIMWLK has the structure of an ordinary second-order Fokker-Planck equation:

\[
\partial_Y W=L_A(\chi^{AB}L_BW).
\]

In divergence form, this can be written as a probability current using the
score \(s_A=L_A\log W\). The corresponding diffusion is second order.

The NLO generalized current normal form used in the dense diagnostic work is

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12 L_A L_B(K_2^{AB}W)
-
\frac16 L_A L_B L_C(K_3^{ABC}W).
\]

The \(K_3\) term contains third-order derivatives. Pawula's theorem warns that
finite Kramers-Moyal truncations with third- or higher-order terms need not
preserve positivity unless the hierarchy continues in a compatible way. Thus,
the algebraic current representation is not by itself a proof of positivity.

The score and Hessian-score representation also assumes \(W>0\), since

\[
s_A=L_A\log W,
\qquad
H_{AB}=L_As_B.
\]

If \(W\) reaches zero or becomes negative in a numerical evolution, those
objects are no longer well-defined without additional regularization or a
different representation.

## Toy Diagnostic

The diagnostic uses a one-dimensional periodic grid
\(\theta\in[0,2\pi)\) and the scalar generalized normal form

\[
\partial_Y W
=
-\partial_\theta(K_1W)
+
\frac12\partial_\theta^2(K_2W)
-
\frac16\partial_\theta^3(K_3W).
\]

Periodic finite-difference derivative matrices are used for
\(\partial_\theta\), \(\partial_\theta^2\), and \(\partial_\theta^3\). The
generator matrix is

\[
G
=
-D_1\,\mathrm{diag}(K_1)
+
\frac12D_2\,\mathrm{diag}(K_2)
-
\frac16D_3\,\mathrm{diag}(K_3).
\]

One short Euler step is then

\[
W_{\rm after}=W+\Delta Y\,GW.
\]

## Cases

The script evaluates four synthetic cases:

1. LO-like diffusion:
   \(K_3=0\), \(K_2>0\).
2. Pure third-order term:
   \(K_1=0\), \(K_2=0\), \(K_3\ne0\).
3. Mixed NLO-like term:
   \(K_2>0\) and small \(K_3\ne0\).
4. Variable-coefficient case:
   nonconstant synthetic \(K_1\), positive \(K_2\), and nonzero \(K_3\).

## Reported Diagnostics

For each case the script reports:

- normalization drift after one Euler step;
- \(\min W_{\rm after}\);
- negative mass after one Euler step,
  \[
  \int d\theta\,\max(-W_{\rm after},0);
  \]
- scaling of negative mass versus \(\Delta Y\);
- a positive maximum-principle diagnostic near grid points where \(W\) is near
  zero;
- a discrete generator matrix off-diagonal sign diagnostic.

For an ordinary continuous-time Markov generator on a finite state space,
off-diagonal matrix elements should be nonnegative. Negative off-diagonal
entries are therefore a useful finite-dimensional warning sign, though not a
proof about the physical continuum NLO Hamiltonian.

## Interpretation

This workflow is a toy diagnostic only. It can demonstrate positivity risk of
finite third-order generalized Fokker-Planck operators. It does not prove
physical NLO JIMWLK positivity or non-positivity, and it does not validate or
invalidate physical coordinate kernels.
