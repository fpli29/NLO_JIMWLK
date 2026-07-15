# Two-Generator NLO Current Sector Summary

This document summarizes the current small-lattice status of the two-generator
NLO blocks. It does not include any full NLO flow implementation and does not
start the three-generator sector.

## \(K_{JSJ}\)

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
(S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}),
\]

\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
\]

Status:

\[
\text{LO-like symmetric score current.}
\]

The sign is fixed by the KLM convention
\(\partial_Y{\cal O}=-H{\cal O}\), \(\partial_YW=-H^\dagger W\), and by the
LO replacement \(K_{JSJ}\to -M/2\).

## \(K_{JSSJ}\)

\[
A_{JSSJ}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{JSSJ}
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}),
\]

\[
v_{JSSJ}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{JSSJ}^{ab}
+
A_{JSSJ}^{ab}s_x^a
\right].
\]

Status:

\[
\text{score + coefficient drift; order-one antisymmetric component measured.}
\]

The diagnostic report measured
\[
\|C-C^T\|/\|C\|=1.3444397747997991
\]
for an \(x\leftrightarrow y\)-symmetric synthetic kernel and
\[
\|C-C^T\|/\|C\|=1.2953094070090567
\]
with both \(x\leftrightarrow y\) and \(z\leftrightarrow z'\) symmetries.
Therefore \(K_{JSSJ}\) should be kept as a generic ordered-current block with
possible commutator drift.

## \(K_{q\bar q}\)

\[
A_{q\bar q}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{q\bar q}
\left[
2{\rm tr}(S^\dagger(z)t^aS(z')t^b)-S_A^{ab}(z)
\right],
\]

\[
v_{q\bar q}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{q\bar q}^{ab}
+
A_{q\bar q}^{ab}s_x^a
\right].
\]

Status:

\[
\text{score + coefficient drift; order-one antisymmetric component measured.}
\]

The \(z'=z\) subtraction identity
\[
2{\rm tr}(S^\dagger t^a S t^b)-S_A^{ab}=0
\]
passes numerically under the current adjoint convention. The diagnostic report
measured
\[
\|C-C^T\|/\|C\|=1.0928913228548376
\]
for an \(x\leftrightarrow y\)-symmetric synthetic kernel and
\[
\|C-C^T\|/\|C\|=0.7911469410486182
\]
with both \(x\leftrightarrow y\) and \(z\leftrightarrow z'\) symmetries.
Therefore \(K_{q\bar q}\) should also be kept as a generic ordered-current
block with possible commutator drift.

## Shared conclusion

All two-generator NLO pieces require at most the score \(s_A\). None of them
require Hessian-score.

Hessian-score first appears in the three-generator sector:

\[
K_{JJSJ},
\qquad
K_{JJSSJ}.
\]

The next step, after the two-generator representation is accepted, is to move
to \(K_{JJSJ}\). That is outside this workflow.

