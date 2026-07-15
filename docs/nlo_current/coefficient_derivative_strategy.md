# Coefficient-Derivative Diagnostic Strategy

This document describes the dense finite-difference backend used only for
small-lattice NLO current diagnostics.

## 1. Why coefficient derivatives are needed

From the NLO velocity formula,

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

the coefficient-derivative arrays are

\[
dK2^A=L_BK_2^{AB},
\]

\[
(LC\_K3)^{AB}=L_CK_3^{ABC},
\]

\[
(LB\_K3)^{AC}=L_BK_3^{ABC},
\]

\[
d2K3^A=L_BL_CK_3^{ABC}.
\]

These are derivatives of Wilson-line coefficient functions. They are separate
from learned or modeled density derivatives such as score and Hessian-score.

## 2. Diagnostic finite-difference backend

For tiny dense tests, coefficient callbacks have the form

```python
def K2_callback(U_fund, S_adj) -> np.ndarray:
    # returns K2 with shape (D,D)

def K3_callback(U_fund, S_adj) -> np.ndarray:
    # returns K3 with shape (D,D,D)
```

The backend uses central finite differences with the same left perturbation
convention as `src/nlo_current/su3_adjoint.py`:

\[
U_x \mapsto \exp(i\epsilon t^b)U_x.
\]

For combined index \(B=(x,b)\),

\[
L_BF(U)
\approx
\frac{F(e^{+i\epsilon t^b}U_x)-F(e^{-i\epsilon t^b}U_x)}{2\epsilon}.
\]

The backend produces

```text
dK2: shape (D,)
LC_K3_ABC: shape (D,D)
LB_K3_ABC: shape (D,D)
d2K3: shape (D,)
```

## 3. Definitions and contractions

For \(K_2^{AB}\),

\[
dK2^A=\sum_B L_BK_2^{AB}.
\]

For \(K_3^{ABC}\),

\[
(LC\_K3)^{AB}=\sum_C L_CK_3^{ABC},
\]

\[
(LB\_K3)^{AC}=\sum_B L_BK_3^{ABC},
\]

\[
d2K3^A=\sum_{B,C}L_BL_CK_3^{ABC}.
\]

The second derivative is ordered as written. Same-site commutator issues are
handled in the canonicalized \(K_3,K_2,K_1\) representation; this backend only
applies finite differences to already assembled coefficient callbacks.

## 4. Non-production warning

\[
\boxed{
\text{This finite-difference backend is not a production coefficient-derivative implementation.}
}
\]

It is only for dense diagnostics, product-rule tests, and velocity sensitivity
studies.
