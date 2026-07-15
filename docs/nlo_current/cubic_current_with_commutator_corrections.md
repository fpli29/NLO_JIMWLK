# Cubic Current With Commutator Corrections

After canonicalizing coincident-site Lie derivative words, cubic Hamiltonian
pieces should be represented schematically in density normal form as

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W).
\]

Canonicalization maps the raw cubic words to
\[
K_3\rightarrow K_{3,\rm canonical},
\]
and induces lower-order corrections
\[
K_2\rightarrow K_2+K_{2,\rm comm},
\]
\[
K_1\rightarrow K_1+K_{1,\rm comm}.
\]

The associated current is
\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

Therefore
\[
v^A
=
K_1^A
-
\frac12
\left[
L_BK_2^{AB}
+
K_2^{AB}s_B
\right]
+
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

Consequences:

- cubic canonical terms require Hessian-score \(H_{BC}=L_Bs_C\);
- quadratic commutator terms require the score \(s_B\);
- linear commutator terms require no density derivatives beyond \(W\);
- coefficient derivatives are derivatives of known Wilson-line coefficient
  functions, not learned density derivatives.

This is a current representation target. The present workflow validates the
symbolic lower-order terms and finite-difference algebra, but does not assemble
a production NLO current.

