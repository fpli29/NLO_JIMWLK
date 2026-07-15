# Physical Density-Side Closure Derivation

## Scope

This note records the algebra used by the tiny-lattice density-side closure
diagnostic. It is a non-production check of the dense physical normal-form
implementation. It does not implement production evolution, train
score/Hessian-score models, or claim positivity or regulator independence.

## Normal-Form Density Operator

Start from the diagnostic generalized Fokker-Planck normal form

\[
\partial_Y W =
-L_A(K_1^A W)
+ \frac12 L_A L_B(K_2^{AB} W)
- \frac16 L_A L_B L_C(K_3^{ABC} W).
\]

Define the density current

\[
J^A =
K_1^A W
-\frac12 L_B(K_2^{AB} W)
+\frac16 L_B L_C(K_3^{ABC} W).
\]

Then, preserving the written Lie-derivative order,

\[
\partial_Y W = -L_A J^A.
\]

Where \(W>0\), define

\[
s_A = L_A \log W,
\qquad
H_{BC}=L_B s_C.
\]

The ordered Hessian-score \(H_{BC}\) is not assumed symmetric. In particular,
coincident-site Lie derivatives can carry commutator effects, so replacing
\(H_{BC}\) by a symmetrized matrix is not part of this validation.

Using

\[
L_B W = s_B W
\]

and

\[
L_B L_C W = (H_{BC}+s_Bs_C)W,
\]

the velocity \(v^A=J^A/W\) is

\[
v^A
=
K_1^A
-\frac12
\left[
L_BK_2^{AB}+K_2^{AB}s_B
\right]
+\frac16
\left[
L_BL_CK_3^{ABC}
+(L_CK_3^{ABC})s_B
+(L_BK_3^{ABC})s_C
+K_3^{ABC}(H_{BC}+s_Bs_C)
\right].
\]

The closure checked by this workflow is therefore

\[
\boxed{
-L_A(v^A W)
=
-L_A(K_1^A W)
+\frac12L_AL_B(K_2^{AB}W)
-\frac16L_AL_BL_C(K_3^{ABC}W)
}.
\]

## Caveats

1. \(H_{BC}=L_Bs_C\) is ordered and must not be assumed symmetric.
2. Coincident-site commutator corrections must already be included in
   \(K_1,K_2,K_3\).
3. The score/Hessian-score representation is local where \(W>0\).
4. Algebraic closure does not imply positivity preservation.
5. Physical finite-grid prescriptions remain diagnostic and non-production.

