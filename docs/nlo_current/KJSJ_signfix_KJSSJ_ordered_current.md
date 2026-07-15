# KJSJ Sign Fix and Ordered KJSSJ Current

This note records the convention checks and the ordered two-generator current
needed before any full NLO JIMWLK flow implementation. It uses the local
left-Lie convention

\[
L_x^a F(U)=\frac{d}{d\epsilon}F(e^{i\epsilon t^a}U_x)\bigg|_{\epsilon=0},
\]

and the adjoint Wilson line

\[
S_A^{ab}(U)=2\,{\rm Re}\,{\rm tr}(t^a U t^b U^\dagger).
\]

With this convention the numerical tests verify
\[
J_L^a=S_A^{ab}J_R^b,\qquad J_R^a=S_A^{ba}J_L^b.
\]

## 1. KLM sign convention

Kovner--Lublinsky--Mulian use the observable-side convention

\[
\frac{d}{dY}{\cal O}=-H{\cal O},
\qquad
\partial_YW=-H^\dagger W.
\]

For a two-generator observable-side term
\[
H_2=C_H^{AB}[U]L_A L_B,
\]
the symmetric density-side diffusion-like tensor is
\[
D^{AB}=-2C_H^{(AB)}.
\]

For a symmetric LO-like divergence kernel,
\[
\chi^{AB}=-C_H^{AB}.
\]

## 2. Corrected \(K_{JSJ}\) result

The corrected generalized, configuration-level \(K_{JSJ}\) contribution is

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
(S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}).
\]

The corresponding score current is
\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B,
\qquad s_B=L_B\log W.
\]

With the LO replacement \(K_{JSJ}\to -M/2\), this gives the usual positive
LO square-difference sign.

## 3. Ordered \(J_L A J_R\) current lemma

Consider the ordered observable-side block

\[
H_{LR}[A]=\int_{x,y}A^{ab}(x,y;U)J_L^a(x)J_R^b(y),
\]

with observable evolution
\[
\partial_Y{\cal O}=-H_{LR}{\cal O}.
\]

Because the charge generators are ordered to the right of the explicit
coefficient, the density-side expression is

\[
(\partial_YW)_{LR}
=
-J_R^b(y)J_L^a(x)
\left[
A^{ab}(x,y;U)W
\right].
\]

Using
\[
J_R^b(y)=S_y^{hb}L_y^h,
\]
a canonical left-divergence form is
\[
(\partial_YW)_{LR}
=
-L_y^h
\left\{
S_y^{hb}
L_x^a
\left[
A^{ab}(x,y;U)W
\right]
\right\}
+
\Delta_{\rm div}[A,W].
\]

For the local convention above, the finite-difference test gives
\[
\sum_h L_y^hS_y^{hb}=0
\]
to numerical precision. Therefore \(\Delta_{\rm div}=0\) for this convention.
The ordered-current component can be written as

\[
J_{LR}^{(y,h)}
=
S_y^{hb}
L_x^a[A^{ab}W],
\]

and the corresponding velocity component is
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

This form keeps the derivative of the ordered coefficient \(A\). It is not the
same as assuming a symmetric diffusion kernel.

## 4. Application to \(K_{JSSJ}\)

The KLM Appendix B identity is used in the collision-free index form

\[
f^{abc}f^{def}S_z^{be}S_{z'}^{cf}
-
N_cS_z^{ad}
=
f^{apc}f^{def}S_z^{pe}
(S_{z'}^{cf}-S_z^{cf}).
\]

Equivalently, in operator notation,

\[
f^{abc}f^{def}J_L^a(x)S_z^{be}S_{z'}^{cf}J_R^d(y)
-
N_cJ_L^a(x)S_z^{ab}J_R^b(y)
\]

\[
=
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
J_L^a(x)J_R^b(y).
\]

Define
\[
A_{JSSJ}^{ab}(x,y;U)
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}).
\]

After converting the right generator to the left basis,
\[
C_{JSSJ}^{(x,a)(y,h)}
=
A_{JSSJ}^{ab}(x,y;U)S_y^{hb}.
\]

If the ordered-block identity is verified, the ordered current is
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

The generic second-order decomposition is
\[
C^{AB}=C^{(AB)}+C^{[AB]},
\]
\[
D^{AB}=-2C^{(AB)}.
\]

The antisymmetric part \(C^{[AB]}\) cannot be discarded unless it is shown to
vanish under the relevant kernel and configuration assumptions. It contributes
through the noncommuting Lie derivatives as a commutator drift.

\[
K_{JSSJ}\text{ needs }s_A\text{ but not }L_As_B.
\]

The small-lattice symmetry diagnostic in
`reports/nlo_current/kjssj_symmetry_report.md` decides whether the tested dense
coefficient is symmetric or whether the commutator drift must be kept.

