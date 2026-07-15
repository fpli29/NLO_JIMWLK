# K_JJSSJ Cubic Ordered Current

This note validates the distinct-site ordered-current structure for
\(K_{JJSSJ}\). It does not implement the full NLO flow and does not solve
coincident-site commutators.

## 1. Hamiltonian structure

Use the five-kernel KLM form. The \(K_{JJSSJ}\) term is

\[
H_{JJSSJ}
=
\int_{w,x,y,z,z'}
K_{JJSSJ}(w;x,y;z,z')f^{acb}
\left[
J_L^d(x)J_L^e(y)S_z^{dc}S_{z'}^{eb}J_R^a(w)
-
J_L^a(w)S_z^{cd}S_{z'}^{be}J_R^d(x)J_R^e(y)
+
\frac13
\left(
J_L^c(x)J_L^b(y)J_L^a(w)
-
J_R^c(x)J_R^b(y)J_R^a(w)
\right)
\right].
\]

The three pieces are:

1. LLR real-like piece \(J_LJ_LS(z)S(z')J_R\).
2. LRR real-like piece \(-J_LS(z)S(z')J_RJ_R\).
3. \(1/3\) virtual piece \((J_LJ_LJ_L-J_RJ_RJ_R)/3\).

This term is cubic in charge generators and therefore requires score and
Hessian-score information.

## 2. Reused distinct-site cubic current lemmas

This workflow reuses the distinct-site cubic lemmas validated for \(K_{JJSJ}\).

### LLR block

For
\[
H_{LLR}[A]=A^{DEA}J_L^DJ_L^EJ_R^A,
\]
with KLM observable convention
\[
\partial_Y{\cal O}=-H{\cal O},
\]
the tested distinct-site current structure is
\[
v_{LLR}^{(w,h)}
=
-S_w^{ha}
\frac{1}{W}L_y^eL_x^d[A^{dea}W].
\]

Expanded,
\[
v_{LLR}^{(w,h)}
=
-S_w^{ha}
\left[
L_y^eL_x^dA^{dea}
+
(L_x^dA^{dea})s_y^e
+
(L_y^eA^{dea})s_x^d
+
A^{dea}(H_{yx}^{ed}+s_y^es_x^d)
\right],
\]
where
\[
H_{yx}^{ed}=L_y^es_x^d.
\]

### LRR block

The distinct-site LRR current lemma reused from the \(K_{JJSJ}\) workflow is
\[
J_{LRR}^{(w,a)}
=
-J_R^e(y)J_R^d(x)[B^{ade}W].
\]

After right-to-left conversion this is tested as
\[
J_{LRR}^{(w,a)}
=
-S_y^{qe}S_x^{pd}L_y^qL_x^p[B^{ade}W]
\]
for distinct lattice sites.

The velocity expansion is
\[
v_{LRR}^{(w,a)}
=
-S_y^{qe}S_x^{pd}
\left[
L_y^qL_x^pB^{ade}
+
(L_x^pB^{ade})s_y^q
+
(L_y^qB^{ade})s_x^p
+
B^{ade}(H_{yx}^{qp}+s_y^qs_x^p)
\right].
\]

### Virtual LLL/RRR block

The virtual block reuses the tested \(1/3\) handling from the \(K_{JJSJ}\)
workflow. The LLL current is
\[
J_{LLL}^{(w,a)}=-L_y^bL_x^c[V^{cba}W].
\]

The RRR virtual sign is tested again with right derivatives and converted at
the outer \(w\) coordinate for distinct sites.

No new signs are assumed without finite-difference checks.

## 3. \(K_{JJSSJ}\) coefficient blocks

### LLR block

The LLR part is
\[
K_{JJSSJ}(w;x,y;z,z')f^{acb}
J_L^d(x)J_L^e(y)S_z^{dc}S_{z'}^{eb}J_R^a(w).
\]

Define
\[
\boxed{
A_{LLR}^{dea}(x,y,w;U)
=
\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}
S_z^{dc}S_{z'}^{eb}.
}
\]

### LRR block

The LRR part is
\[
-
K_{JJSSJ}(w;x,y;z,z')f^{acb}
J_L^a(w)S_z^{cd}S_{z'}^{be}J_R^d(x)J_R^e(y).
\]

Define
\[
\boxed{
B_{LRR}^{ade}(w,x,y;U)
=
-\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}
S_z^{cd}S_{z'}^{be}.
}
\]

### Virtual LLL/RRR block

The virtual piece is
\[
\frac13K_{JJSSJ}(w;x,y;z,z')f^{acb}
\left[
J_L^c(x)J_L^b(y)J_L^a(w)
-
J_R^c(x)J_R^b(y)J_R^a(w)
\right].
\]

Define
\[
\boxed{
V^{cba}(x,y,w)
=
\frac13\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')f^{acb}.
}
\]

The \(1/3\) factor is included in the coefficient builder and tested by a norm
ratio.

## 4. Kernel symmetry and antisymmetry

KLM states that \(K_{JJSSJ}\) is antisymmetric under simultaneous interchange
\[
x\leftrightarrow y,
\qquad
z\leftrightarrow z'.
\]

The primary synthetic diagnostic imposes
\[
K(w;x,y;z,z')=-K(w;y,x;z',z).
\]

An unconstrained synthetic kernel is also used as a stress comparison. Reports
state which convention is being tested.

## 5. Density-derivative requirements

Like \(K_{JJSJ}\), the \(K_{JJSSJ}\) term is cubic in charge generators. It
requires
\[
s_A=L_A\log W
\]
and
\[
H_{AB}=L_As_B=L_AL_B\log W.
\]

It may also require coefficient derivatives
\[
L_AA,\qquad L_AL_BA,
\]
which are derivatives of known Wilson-line coefficient functions, not learned
density derivatives.

\[
\boxed{
K_{JJSSJ}\text{ requires Hessian-score.}
}
\]

## 6. Scope limitation

This workflow validates distinct-site ordered-current structure only. It does
not claim coincident-site ordering is resolved. Coincident-site commutators must
be handled in a separate workflow:

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

