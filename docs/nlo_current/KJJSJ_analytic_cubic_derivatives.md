# KJJSJ Analytic Cubic Coefficient Derivatives

## Scope

This note derives the analytic/local coefficient derivatives for the
\(K_{JJSJ}\) sector only. \(K_{JJSSJ}\) remains pending.

The physical coordinate kernel is independent of Wilson lines. Lie derivatives
therefore act only on Wilson-line coefficient structures and on the
right-to-left conversion factors used before normal-form canonicalization.

## Normal-Form Assembly Path

The code path is:

```text
raw physical KJJSJ
-> physical_cubic_conventions.klm_normalized_cubic_kernel
-> physical_kernel_adapter output KJJSJ
-> three_generator_terms KJJSJ blocks
-> cubic_commutator_terms right-to-left conversion
-> lie_word_algebra canonicalization
-> normal-form K3, K2_comm, K1_comm
-> analytic derivatives
```

Analytic derivatives act on the normalized coefficient already used by the
skeleton. The \((-i)\) normalization is not re-applied in derivative code.

## Distinct-Site Blocks

The LLR block is

\[
A_{LLR}^{dea}(x,y,w)
=
\int_z K_{JJSJ}(w;x,y;z) f^{bde} S_z^{ba}.
\]

The LRR block is

\[
B_{LRR}^{ade}(w,x,y)
=
-\int_z K_{JJSJ}(w;x,y;z) f^{bde} S_z^{ab}.
\]

The virtual block is

\[
V^{deb}(x,y,w)=
\frac13\int_z K_{JJSJ}(w;x,y;z) f^{bde},
\]

and enters as an `LLL - RRR` combination after the right-to-left conversion of
the RRR term.

## First Derivatives

For the LLR coefficient,

\[
L_C A_{LLR}^{dea}
=
\int_z K_{JJSJ} f^{bde} L_C S_z^{ba}.
\]

The final left-basis LLR coefficient is

\[
C_{LLR}^{deh}=A_{LLR}^{dea}S_w^{ha},
\]

so its derivative is

\[
L_C C_{LLR}^{deh}
=(L_C A_{LLR}^{dea})S_w^{ha}
+A_{LLR}^{dea}(L_C S_w^{ha}).
\]

For LRR,

\[
L_C B_{LRR}^{ade}
=
-\int_z K_{JJSJ} f^{bde} L_C S_z^{ab},
\]

with left-basis coefficient

\[
C_{LRR}^{apq}=B_{LRR}^{ade}S_x^{pd}S_y^{qe}.
\]

The derivative uses the ordered product rule over \(B_{LRR}\), \(S_x\), and
\(S_y\).

The virtual LLL block is Wilson-line independent, so its first derivative is
zero. The virtual RRR block contains three adjoint conversion factors and is
differentiated by the same ordered product rule. After the current
canonicalization, the validated diagnostic has structurally zero virtual `K3`
derivative contractions; expected-zero virtual residuals are therefore judged
by absolute residuals, not by relative residual alone.

## Ordered Second Derivatives

The ordered second derivative for LLR is

\[
L_B L_C A_{LLR}^{dea}
=
\int_z K_{JJSJ} f^{bde} L_B L_C S_z^{ba}.
\]

The full left-basis block follows

\[
L_B L_C(A S)
=(L_B L_C A)S
+(L_C A)(L_B S)
+(L_B A)(L_C S)
+A(L_B L_C S).
\]

LRR and virtual RRR use the corresponding ordered product rule over three
factors. Same-site derivative ordering is preserved; \(L_B L_C\) is not
replaced by \(L_C L_B\). Distinct-site derivatives commute only when their
local Kronecker support makes the commutator vanish.

## Canonicalization and Commutator Corrections

The analytic implementation differentiates the same raw block coefficients and
then passes the derivative terms through the same canonicalization used by the
skeleton. This preserves the normal-form \(K_3\) convention and the lower-order
commutator corrections.

The \(K_{2,\rm comm}\) derivative is obtained from the quadratic terms produced
by canonicalization. The current \(K_{1,\rm comm}\) status is classified from
the canonical linear terms instead of assumed.
