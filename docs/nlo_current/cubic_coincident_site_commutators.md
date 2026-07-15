# Cubic Coincident-Site Commutators

This note defines the symbolic commutator algebra needed before the cubic
ordered \(K_{JJSJ}+K_{JJSSJ}\) sector can be assembled into any canonical
current representation. It does not implement the full NLO flow.

## 1. Why distinct-site validation is insufficient

The previous cubic workflows validated ordered currents assuming distinct
lattice sites, where
\[
[L_x^a,L_y^b]=0
\quad\text{for}\quad x\neq y.
\]

Full lattice sums include coincident coordinates:
\[
x=y,\qquad x=w,\qquad y=w,\qquad x=y=w.
\]

On the same site,
\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

Therefore ordered products such as
\[
L_x^aL_x^bL_w^c,
\qquad
L_x^aL_w^cL_x^b,
\qquad
L_x^aL_x^bL_x^c
\]
cannot be freely symmetrized or reordered without generating lower-order terms.

## 2. Canonical ordering

Use combined derivative indices
\[
A=(x,a).
\]

The canonical order is lexicographic:
\[
(x,a)<(y,b)
\]
if either \(x<y\), or \(x=y\) and \(a<b\).

For neighboring same-site derivatives, the canonicalization rule is
\[
L_x^aL_x^b
=
L_x^bL_x^a
+
f^{abc}L_x^c
\]
when the pair must be swapped from \((a,b)\) to \((b,a)\).

For different sites,
\[
L_x^aL_y^b=L_y^bL_x^a,
\qquad x\neq y.
\]

This sign matches the implemented left perturbation convention
\[
L^aF(U)=dF(e^{i\epsilon t^a}U)/d\epsilon|_{\epsilon=0}.
\]

## 3. Operator-level normal ordering

For an ordered derivative word
\[
L_A L_B L_C
\]
acting on a coefficient-density product \(F[U]=K[U]W[U]\), canonicalization
returns
\[
L_A L_B L_C F
=
\sum_{\alpha} c_\alpha L_{\alpha_1}L_{\alpha_2}L_{\alpha_3}F
+
\sum_{\beta} d_\beta L_{\beta_1}L_{\beta_2}F
+
\sum_{\gamma} e_\gamma L_{\gamma_1}F.
\]

The cubic piece contributes to the Hessian-score current. The induced quadratic
and linear pieces are lower-order commutator corrections.

## 4. Classification of commutator-induced terms

The canonicalization output is classified into
\[
K_3^{ABC}
\quad\text{cubic},
\]
\[
K_{2,\rm comm}^{AB}
\quad\text{quadratic commutator correction},
\]
\[
K_{1,\rm comm}^{A}
\quad\text{linear commutator correction}.
\]

In canonical density-side normal form,
\[
J^A_{\rm cubic}
=
\frac16L_BL_C(K_3^{ABC}W),
\]
\[
J^A_{\rm quad,comm}
=
-\frac12L_B(K_{2,\rm comm}^{AB}W),
\]
\[
J^A_{\rm lin,comm}
=
K_{1,\rm comm}^{A}W.
\]

The current signs follow the schematic density normal form and must still be
matched to the full NLO coefficient assembly before production use.

## 5. Scope

\[
\boxed{
\text{This workflow validates algebraic commutator corrections only.}
}
\]

It does not fold the corrections into a production evolution code path and does
not add any score or Hessian-score model training.

