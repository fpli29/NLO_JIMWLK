# Physical Density Closure Failure Modes

This note lists failure modes for the tiny-lattice density-side closure
diagnostic. It is not a production-evolution or positivity document.

## Finite-Difference Cancellation

The direct density operator contains nested finite differences through third
order. Too small a step amplifies roundoff; too large a step leaves truncation
error. The diagnostic therefore records a step scan and looks for a stable
window rather than forcing monotonic convergence to roundoff.

## Inconsistent Score/Hessian Backends

The score \(s_A\) and ordered Hessian-score \(H_{BC}=L_Bs_C\) must come from the
same \(\log W\). Mixing backends or step sizes can mimic a current-assembly
error.

## Frozen-Velocity Outer Derivative

The divergence \(-L_A(v^A W)\) requires differentiating both \(v^A(U)\) and
\(W(U)\). Reusing a velocity computed only at the unperturbed configuration
misses coefficient and score/Hessian variations.

## Missing Coefficient Derivatives

The velocity contains \(L_BK_2^{AB}\), \(L_BL_CK_3^{ABC}\), and first-derivative
contractions of \(K_3\). Omitting these terms generically breaks closure.

## Missing Commutator Corrections

Coincident-site cubic words generate lower-order \(K_2\) and \(K_1\)
corrections. Disabling them changes the normal-form coefficients and can worsen
closure when those words are exercised.

## Raw Versus Normalized Cubic Convention

The physical adapter supplies KLM-normalized cubic coefficients:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

Using raw-complex cubic coefficients in the normal-form current produces the
wrong complex character.

## Singular Coordinate Entries

The physical kernel layer uses explicit non-production singularity and
quadrature policies. Changing the policy between the direct and current sides
invalidates a closure comparison.

## Dtype and Casting Errors

Complex diagnostic values must remain complex until an explicit realness check.
Silently casting complex arrays to real can hide a cubic-convention error.

