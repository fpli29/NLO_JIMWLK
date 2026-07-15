# Non-Production NLO Current Skeleton Design

This design describes a dense small-lattice diagnostic skeleton. It is not a
production NLO evolution implementation.

## 1. Normal form

Use the density normal form
\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W).
\]

The corresponding current is
\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

The corresponding velocity is
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
\right],
\]
where
\[
s_A=L_A\log W,
\qquad
H_{BC}=L_Bs_C.
\]

## 2. Sector assignment

### \(K_{JSJ}\)

\[
K_{JSJ}: K_2.
\]

This is an LO-like score current.

### \(K_{JSSJ}\)

\[
K_{JSSJ}: K_2.
\]

This is a generic second-order current. It requires score and coefficient
derivatives, but not Hessian-score.

### \(K_{q\bar q}\)

\[
K_{q\bar q}: K_2.
\]

This is also a generic second-order current. It requires score and coefficient
derivatives, but not Hessian-score.

### \(K_{JJSJ}\)

\[
K_{JJSJ}:
K_3^{JJSJ}
+
K_{2,\rm comm}^{JJSJ}
+
K_{1,\rm comm}^{JJSJ}.
\]

The commutator diagnostic found
\[
K_{2,\rm comm}^{JJSJ}\neq0,
\]
and
\[
K_{1,\rm comm}^{JJSJ}=0
\]
in the synthetic diagnostic. The interface still allows nonzero \(K_1\).

### \(K_{JJSSJ}\)

\[
K_{JJSSJ}:
K_3^{JJSSJ}
+
K_{2,\rm comm}^{JJSSJ}
+
K_{1,\rm comm}^{JJSSJ}.
\]

The commutator diagnostic found both nonzero:
\[
K_{2,\rm comm}^{JJSSJ}\neq0,
\qquad
K_{1,\rm comm}^{JJSSJ}\neq0.
\]

## 3. Data model

The non-production dense structure is

```python
@dataclass
class NLOCurrentTerms:
    K1: np.ndarray
    K2: np.ndarray
    K3: np.ndarray
    metadata: dict
```

For \(N_{\rm site}=N\), \(N_c^2-1=8\), and combined dimension \(D=8N\):

```text
K1: shape (D,)
K2: shape (D,D)
K3: shape (D,D,D)
```

Metadata records sector contributions, norms, commutator toggles, warnings, and
diagnostic assumptions.

## 4. Velocity evaluator scope

The skeleton can evaluate
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

Coefficient derivative arrays are explicit diagnostic inputs:

\[
L_BK_2^{AB},
\quad
L_BL_CK_3^{ABC},
\quad
L_CK_3^{ABC},
\quad
L_BK_3^{ABC}.
\]

If they are omitted, the evaluator treats them as zero and records warnings.
No production automatic differentiation is attempted.

## 5. Non-production warning

\[
\boxed{
\text{This skeleton is not a production NLO evolution implementation.}
}
\]

It validates interfaces, signs, term assembly, shapes, and diagnostic velocity
evaluation only.

