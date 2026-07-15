# Analytic Coefficient Derivative Failure Modes

## Wrong Adjoint-Index Orientation

The calibrated rule is \(L^hS_A^{ab}=f^{hac}S_A^{cb}\). Acting on the second
adjoint index or flipping the sign produces sector-level FD mismatches.

## Wrong Generator Sign

The local convention is \(L^aU=i t^aU\) and
\(L^aU^\dagger=-iU^\dagger t^a\). Reversing either sign breaks the qbarq trace
derivative.

## Derivative-Order Reversal

Ordered second derivatives are not symmetric. Reversing \(L_B L_C\) changes
same-site terms by commutators.

## Accidental Hessian Symmetrization

Density closure uses ordered Hessian-score \(H_{BC}=L_Bs_C\). The analytic
coefficient path must not symmetrize derivative contractions to compensate for
an ordered-Hessian mismatch.

## Duplicate or Missing Cubic Normalization

Raw physical cubic kernels carry the WORKNLO complex convention. The skeleton
uses KLM-normalized coefficients after one explicit \((-i)\) factor. Applying
that factor twice, or not at all, changes the complex character of the cubic
current.

## Omitted Commutator Corrections

Cubic normal-form tensors include lower-order corrections from same-site
commutators. Any analytic cubic backend must differentiate the final
normal-form tensors, not only the pre-canonical cubic words.

## Finite-Grid Policy Mismatch

FD and analytic comparisons must use the same physical kernels, finite
coordinate set, singularity policy, and `KJSJIntegrationPolicy`.

## Silent Backend Fallback

`backend="analytic"` must not silently call the FD oracle for pending sectors.
The current implementation supports KJJSJ analytically and still raises for
pending KJJSSJ unless `backend="hybrid_local_fd"` is explicitly selected.

## Zero-Tensor Relative Error

For structurally zero cubic derivative contractions in two-generator sector
filters, absolute residuals are meaningful and relative residuals can be
unstable.

## Complex-to-Real Cast

Complex intermediate values are expected in qbarq trace derivatives and raw
cubic diagnostics. They must remain complex until an explicit realness check.
