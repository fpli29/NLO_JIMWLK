# K_JJSJ Cubic Ordered Current

This note starts the three-generator sector with \(K_{JJSJ}\) only. It does not
implement the full NLO current and does not start \(K_{JJSSJ}\).

## 1. Hamiltonian structure

Use the five-kernel KLM form. The \(K_{JJSJ}\) term is

\[
H_{JJSJ}
=
\int_{w,x,y,z}
K_{JJSJ}(w;x,y;z) f^{bde}
\left[
J_L^d(x)J_L^e(y)S_z^{ba}J_R^a(w)
-
J_L^a(w)S_z^{ab}J_R^d(x)J_R^e(y)
+
\frac13
\left(
J_L^d(x)J_L^e(y)J_L^b(w)
-
J_R^d(x)J_R^e(y)J_R^b(w)
\right)
\right].
\]

The three pieces are:

1. LLR real-like piece \(J_LJ_LSJ_R\).
2. LRR real-like piece \(-J_LSJ_RJ_R\).
3. \(1/3\) virtual piece \((J_LJ_LJ_L-J_RJ_RJ_R)/3\).

This block is cubic in charge generators, so it introduces score and
Hessian-score dependence.

## 2. Ordered cubic current lemma: LLR block

Define

\[
H_{LLR}[A]
=
\int A^{dea}(x,y,w;U)
J_L^d(x)J_L^e(y)J_R^a(w).
\]

With observable evolution
\[
\partial_Y{\cal O}=-H_{LLR}{\cal O},
\]
the density-side operator is
\[
(\partial_YW)_{LLR}
=
+J_R^a(w)J_L^e(y)J_L^d(x)
\left[
A^{dea}W
\right].
\]

The sign follows from three integrations by parts and the outer KLM minus sign,
and is checked numerically in `tests/nlo_current/test_cubic_ordered_current.py`.

Using
\[
J_R^a(w)=S_w^{ha}L_w^h,
\]
and the verified identity
\[
L_w^hS_w^{ha}=0,
\]
the current component is
\[
\boxed{
J_{LLR}^{(w,h)}
=
- S_w^{ha}
L_y^eL_x^d[A^{dea}W].
}
\]

Thus
\[
(\partial_YW)_{LLR}
=
-L_w^hJ_{LLR}^{(w,h)},
\]
and
\[
\boxed{
v_{LLR}^{(w,h)}
=
- S_w^{ha}
\frac{1}{W}L_y^eL_x^d[A^{dea}W].
}
\]

Expanding,
\[
\frac{1}{W}L_y^eL_x^d[A^{dea}W]
=
L_y^eL_x^dA^{dea}
+
(L_x^dA^{dea})s_y^e
+
(L_y^eA^{dea})s_x^d
+
A^{dea}\left(L_y^es_x^d+s_y^es_x^d\right).
\]

Therefore
\[
\boxed{
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
A^{dea}\left(H_{yx}^{ed}+s_y^es_x^d\right)
\right],
}
\]
where
\[
H_{yx}^{ed}=L_y^es_x^d.
\]

This is the first place Hessian-score enters.

## 3. Ordered cubic current lemma: LRR block

Define
\[
H_{LRR}[B]
=
\int B^{ade}(w,x,y;U)
J_L^a(w)J_R^d(x)J_R^e(y).
\]

For the \(K_{JJSJ}\) real-like LRR term, \(B\) already includes the explicit
minus sign in the Hamiltonian coefficient.

The density-side operator is
\[
(\partial_YW)_{LRR}
=
+J_R^e(y)J_R^d(x)J_L^a(w)[B^{ade}W].
\]

For the distinct-coordinate smoke tests, the derivatives commute with
\(L_w^a\), giving the current component
\[
\boxed{
J_{LRR}^{(w,a)}
=
-J_R^e(y)J_R^d(x)[B^{ade}W].
}
\]

After converting the two right generators,
\[
J_R^d(x)=S_x^{pd}L_x^p,\qquad
J_R^e(y)=S_y^{qe}L_y^q,
\]
the tested distinct-coordinate expression is
\[
J_{LRR}^{(w,a)}
=
-S_y^{qe}S_x^{pd}L_y^qL_x^p[B^{ade}W].
\]

The corresponding velocity expansion is
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

The sign and ordering above are checked by finite differences for distinct
sites. Coincident-site commutator/conversion terms are intentionally not
promoted to a production formula in this workflow.

## 4. Ordered cubic current lemma: LLL and RRR virtual blocks

The virtual block is
\[
H_{\rm virt}
=
\frac13
\int V_0^{deb}(x,y,w;U)
\left[
J_L^d(x)J_L^e(y)J_L^b(w)
-
J_R^d(x)J_R^e(y)J_R^b(w)
\right],
\]
where the implemented coefficient includes the factor
\[
V^{deb}(x,y,w;U)
=
\frac13\int_zK_{JJSJ}(w;x,y;z)f^{bde}.
\]

For the LLL ordered block,
\[
H_{LLL}[V]=V^{deb}L_x^dL_y^eL_w^b,
\]
the density-side contribution is
\[
(\partial_YW)_{LLL}=+L_w^bL_y^eL_x^d[V^{deb}W],
\]
with current
\[
J_{LLL}^{(w,b)}=-L_y^eL_x^d[V^{deb}W].
\]

For the RRR virtual block, the observable-side sign is opposite:
\[
-V^{deb}R_x^dR_y^eR_w^b.
\]
Its density-side sign is therefore
\[
(\partial_YW)_{RRR}=-R_w^bR_y^eR_x^d[V^{deb}W].
\]
The workflow only checks the sign/shape mapping for this block; a full
coincident-site production conversion is left for a later full three-generator
implementation.

The symmetric canonical cubic normal form,
\[
\partial_YW\supset -\frac16L_AL_BL_C(K_3^{ABC}W),
\]
has current
\[
J^A_{\rm cubic}
=
\frac16L_BL_C(K_3^{ABC}W),
\]
and velocity
\[
v^A_{\rm cubic}
=
\frac16
\left[
L_BL_CK_3^{ABC}
+
(L_CK_3^{ABC})s_B
+
(L_BK_3^{ABC})s_C
+
K_3^{ABC}(H_{BC}+s_Bs_C)
\right].
\]
This note does not replace the ordered KLM blocks by a fully symmetrized
canonical tensor.

## 5. Application to \(K_{JJSJ}\)

The three observable-side coefficient blocks are:

### LLR block

\[
A_{LLR}^{dea}(x,y,w;U)
=
\int_z
K_{JJSJ}(w;x,y;z)
f^{bde}S_z^{ba}.
\]

### LRR block

The Hamiltonian contains
\[
-K_{JJSJ}(w;x,y;z)f^{bde}
J_L^a(w)S_z^{ab}J_R^d(x)J_R^e(y).
\]

Define
\[
B_{LRR}^{ade}(w,x,y;U)
=
-\int_z
K_{JJSJ}(w;x,y;z)
f^{bde}S_z^{ab}.
\]

### Virtual LLL/RRR block

\[
V^{deb}(x,y,w;U)
=
\frac13\int_z
K_{JJSJ}(w;x,y;z)f^{bde}.
\]

The virtual contribution is
\[
V^{deb}
\left[
J_L^d(x)J_L^e(y)J_L^b(w)
-
J_R^d(x)J_R^e(y)J_R^b(w)
\right].
\]

The \(1/3\) factor is kept in the coefficient builder and tested by norm ratio.

## 6. ML object requirements

The \(K_{JJSJ}\) term requires
\[
s_A=L_A\log W,
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
K_{JJSJ}\text{ is the first term that requires Hessian-score.}
}

