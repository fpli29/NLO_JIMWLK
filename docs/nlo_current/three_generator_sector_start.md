# Three-Generator Sector Start

This note records the first cubic workflow status. It only covers
\(K_{JJSJ}\). It does not implement the full NLO current and does not start
\(K_{JJSSJ}\).

## Two-generator sector status

\[
K_{JSJ}: \text{LO-like score current}.
\]

\[
K_{JSSJ},\ K_{q\bar q}: \text{generic ordered currents with score and coefficient drift}.
\]

The two-generator sector requires at most the score \(s_A=L_A\log W\). It does
not require Hessian-score.

## First three-generator term

\[
K_{JJSJ}
\]

is the first block that requires
\[
H_{AB}=L_As_B=L_AL_B\log W.
\]

The finite-difference diagnostics found nonzero Hessian-score terms for the
toy coupled density:

- symmetric synthetic \(K(w;x,y;z)\): Hessian probe `5.8478200094802535e-03`;
- antisymmetric synthetic \(K(w;x,y;z)\): Hessian probe `1.7848910623813330e-02`.

Both \(x,y\) kernel symmetry conventions were tested. The antisymmetric
synthetic kernel makes the LLR coefficient symmetric under the combined
exchange \((x,d)\leftrightarrow(y,e)\); the symmetric synthetic kernel gives
plain color antisymmetry in \(d,e\) at fixed \(x,y\).

## LLR working formula

For
\[
H_{LLR}[A]=A^{dea}J_L^d(x)J_L^e(y)J_R^a(w),
\]
the tested density-side current is
\[
J_{LLR}^{(w,h)}
=
-S_w^{ha}L_y^eL_x^d[A^{dea}W].
\]

The velocity is
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
\right].
\]

The distinct-site sign check passed with relative residuals below `4e-10`.

## LRR working formula

For
\[
H_{LRR}[B]=B^{ade}J_L^a(w)J_R^d(x)J_R^e(y),
\]
the tested distinct-site current is
\[
J_{LRR}^{(w,a)}
=
-J_R^e(y)J_R^d(x)[B^{ade}W].
\]

After converting right generators,
\[
J_{LRR}^{(w,a)}
=
-S_y^{qe}S_x^{pd}L_y^qL_x^p[B^{ade}W],
\]
for the distinct-coordinate validation setup.

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

The distinct-site sign check passed with zero residual at the reported
precision.

## Virtual \(1/3\) term

The implemented virtual coefficient is
\[
V^{deb}(x,y,w)
=
\frac13\int_zK_{JJSJ}(w;x,y;z)f^{bde}.
\]

The \(1/3\) factor is included in `kjjsj_V_virtual_from_kernel`. The coefficient
test verifies that removing the factor changes the virtual norm by exactly a
factor of three.

The LLL ordered virtual current is
\[
J_{LLL}^{(w,b)}=-L_y^eL_x^d[V^{deb}W].
\]

The RRR virtual sign was also checked in a distinct-site right-to-left smoke
test. The tested sign residuals were below `7e-10`.

## Remaining work

The tested formulas are small-lattice, distinct-site ordered-block checks.
Coincident-site commutator handling and a production-volume conversion are not
implemented here.

The next three-generator block is
\[
K_{JJSSJ}.
\]

That block remains untouched and should be handled in a separate workflow.

