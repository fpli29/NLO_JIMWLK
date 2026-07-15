# NLO JIMWLK Generalized Current: Executive Summary

This is a documentation-only consolidation of the dense small-lattice,
non-production NLO current work. It does not implement physical KLM kernels,
production evolution, or score/Hessian-score model training.

## Main Result

The validated density-side normal form is

\[
\partial_Y W=-L_A J^A_{\rm NLO},
\]

with current

\[
J^A_{\rm NLO}
=
K_1^A W
-\frac12 L_B(K_2^{AB}W)
+\frac16 L_B L_C(K_3^{ABC}W).
\]

Equivalently, after dividing by \(W\), the current velocity contains:

- a direct one-generator drift \(K_1^A\);
- two-generator score terms from \(K_2^{AB}\);
- three-generator score and Hessian-score terms from \(K_3^{ABC}\);
- coefficient-derivative drift terms when the kernels depend on Wilson lines.

The dense skeleton stores the normal-form tensors as

\[
K_1:(D),\qquad K_2:(D,D),\qquad K_3:(D,D,D),
\qquad D=8N_{\rm site}.
\]

## Difference From LO

At LO, the JIMWLK Hamiltonian is already in divergence form:

\[
\partial_Y W=L_A(\chi^{AB}L_BW)
       =-L_A(v^AW),
\qquad
v^A=-\chi^{AB}s_B,
\]

where \(s_A=L_A\log W\). The derivative-of-\(\chi\) term is absorbed by the
divergence-form current, so LO needs only the score.

At NLO, the two-generator sectors still require score-level information, with
additional coefficient drift for Wilson-line-dependent coefficients. The cubic
sectors require

\[
H_{AB}=L_A s_B=L_A L_B\log W,
\]

so a complete NLO current needs Hessian-score information or an equivalent
contracted-Hessian estimator.

## Validated Sector Map

| sector | normal-form contribution | density information required | status |
|---|---|---|---|
| \(K_{JSJ}\) | \(K_2\) | score | Appendix A passed |
| \(K_{JSSJ}\) | \(K_2\) | score + coefficient drift | Appendix A passed |
| \(K_{q\bar q}\) | \(K_2\) | score + coefficient drift | Appendix A passed |
| \(K_{JJSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score | Appendix A passed |
| \(K_{JJSSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score | Appendix A passed |

The observable-side Appendix A convention is calibrated by \(K_{JSJ}\) as
\(H_{\rm sector}s\). For the one-\(f\) cubic sectors, the tested
Hermitian-generator direct-action convention obeys

\[
\text{TeX target}=(-i)\times\text{raw direct action}.
\]

## Validation Milestone

All five Appendix A dipole sectors are implemented as non-production targets
and passed dense synthetic tests:

| sector | max residual |
|---|---:|
| \(K_{JSJ}\) | `1.6172982426630268e-16` |
| \(K_{JSSJ}\) | `1.3877787807814457e-16` |
| \(K_{q\bar q}\) | `5.6325111571547335e-17` |
| \(K_{JJSJ}\) | `1.1801832636420706e-15` |
| \(K_{JJSSJ}\) | `1.9087382418229414e-15` |

This validates the full five-kernel NLO Hamiltonian action on the dipole
against Appendix A in the current dense, synthetic, observable-side test setup.

## Remaining Work

The next step is physical kernel integration in a separate workflow:

- implement barred/unbarred physical coordinate kernels with explicit regulator
  and subtraction policy;
- keep the first integration on tiny dense lattices with no physics claims;
- replace or augment the finite-difference coefficient-derivative diagnostic
  with scalable analytic, sparse, or automatic-differentiation rules;
- design a score/Hessian-score strategy, preferably allowing contracted
  quantities such as \(K_3^{ABC}H_{BC}\) without materializing a full Hessian;
- only after those pieces are validated, attempt a non-production NLO flow
  experiment for stability and sanity checks.

This work is not production-ready and is intentionally separate from production
evolution code.
