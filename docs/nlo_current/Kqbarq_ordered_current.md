# K_qbarq Ordered Two-Generator Current

This note adds the quark-pair two-generator block to the small-lattice NLO
current checks. It does not implement the full NLO flow and does not touch the
three-generator sector.

## 1. Hamiltonian structure

The quark-pair contribution in the KLM Hamiltonian has the two-generator
structure

\[
H_{q\bar q}
=
\int_{x,y,z,z'}
K_{q\bar q}(x,y;z,z')
\left[
2J_L^a(x)
{\rm tr}\!\left(S^\dagger(z)t^aS(z')t^b\right)
J_R^b(y)
-
J_L^a(x)S_A^{ab}(z)J_R^b(y)
\right].
\]

For configuration-level current use the generalized/nonsinglet barred kernel,

\[
K_{q\bar q}\rightarrow \bar K_{q\bar q}.
\]

For singlet dipole validation, use the unbarred kernel.

## 2. Ordered \(J_L A J_R\) coefficient

Define

\[
A_{q\bar q}^{ab}(x,y;U)
=
\int_{z,z'}
\bar K_{q\bar q}(x,y;z,z')
\left[
2{\rm tr}\!\left(S^\dagger(z)t^aS(z')t^b\right)
-
S_A^{ab}(z)
\right].
\]

Then the Hamiltonian block is

\[
H_{q\bar q}
=
\int_{x,y}
A_{q\bar q}^{ab}(x,y;U)
J_L^a(x)J_R^b(y).
\]

The subtraction is consistent with the adjoint convention because when
\(z'=z\),

\[
2{\rm tr}\!\left(S^\dagger(z)t^aS(z)t^b\right)
-
S_A^{ab}(z)
=0,
\]

as verified by the coefficient tests.

## 3. Ordered \(J_L A J_R\) lemma

From the previous workflow, for

\[
H_{LR}[A]=A^{ab}J_L^aJ_R^b
\]

with

\[
\partial_Y{\cal O}=-H_{LR}{\cal O},
\]

and

\[
J_R^b(y)=S_y^{hb}L_y^h,
\]

the density-side current is

\[
J_{LR}^{(y,h)}
=
S_y^{hb}
L_x^a[A^{ab}W],
\]

assuming

\[
L_y^hS_y^{hb}=0
\]

under the current convention.

Therefore

\[
v_{LR}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA^{ab}
+
A^{ab}s_x^a
\right].
\]

For the quark term,

\[
\boxed{
v_{q\bar q}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{q\bar q}^{ab}
+
A_{q\bar q}^{ab}s_x^a
\right].
}
\]

## 4. Required score information

Because \(K_{q\bar q}\) is second order in charge generators,

\[
\boxed{
K_{q\bar q}\text{ needs the score }s_A,
\text{ but not the Hessian-score }L_As_B.
}
\]

As with \(K_{JSSJ}\), the coefficient should not be assumed symmetric. The
symmetry or asymmetry of the left-basis coefficient must be measured.

## 5. Left-basis coefficient

Convert to the left basis:

\[
C_{q\bar q}^{(x,a)(y,h)}
=
A_{q\bar q}^{ab}(x,y;U)S_y^{hb}.
\]

Then decompose

\[
C^{AB}=C^{(AB)}+C^{[AB]}.
\]

If \(C^{[AB]}\neq0\), keep the generic second-order current and the associated
commutator drift. The numerical report in
`reports/nlo_current/kqbarq_symmetry_report.md` is a diagnostic, not a proof.

