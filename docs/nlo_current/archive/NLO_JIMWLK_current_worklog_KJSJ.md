# NLO Generalized Probability Current: Worklog and First-Term Derivation

This document records the current plan, formulas to derive, checkpoints requiring verification, and the detailed derivation of the first NLO JIMWLK term, \(K_{JSJ}\).

The purpose is to build a bridge from the existing LO deterministic probability-flow paper to a possible NLO generalized probability-current formulation.

---

## 1. Current Objective

We want to test whether the NLO JIMWLK Hamiltonian can be rewritten as a generalized probability-current equation,

\[
\partial_Y W[U] = - L_A J^A[U,Y],
\]

and then define

\[
v^A_{\rm gen}[U,Y] = \frac{J^A[U,Y]}{W[U,Y]}.
\]

At LO, the current is

\[
v^A_{\rm LO}=-\chi^{AB}s_B,
\qquad
s_B=L_B\log W.
\]

At NLO, terms with two charge generators \(J\) should still require only the score \(s_A\). Terms with three charge generators \(J\) should require the Hessian-score object

\[
L_A s_B = L_A L_B \log W.
\]

Thus the expected extension is

\[
\boxed{
\text{LO: score current}
\quad\longrightarrow\quad
\text{NLO: score + Hessian-score generalized current}.
}
\]

---

## 2. Required Derivations

We need to derive, term by term, how the NLO Hamiltonian contributes to the density-side current.

The NLO Hamiltonian contains the following structures:

| Term | Charge-generator order | Expected density derivative order | Current needs |
|---|---:|---:|---|
| \(K_{JSJ}\) | 2 | 2 | score \(s\) |
| \(K_{JSSJ}\) | 2 | 2 | score \(s\), possible kernel-derivative drift |
| \(K_{q\bar q}\) | 2 | 2 | score \(s\), possible kernel-derivative drift |
| \(K_{JJSJ}\) | 3 | 3 | score \(s\) and Hessian-score \(Ls\) |
| \(K_{JJSSJ}\) | 3 | 3 | score \(s\) and Hessian-score \(Ls\) |
| \(K_{JJJ}\) or \(1/3\) virtual pieces | 3 | 3 | score \(s\) and Hessian-score \(Ls\) |

The immediate plan is:

1. Convert each term from \(J_L,J_R\) to a single left basis \(L_A\).
2. Identify the coefficient tensor \(C^{AB}[U]\) or \(C^{ABC}[U]\).
3. Take the adjoint to obtain the density-side operator.
4. Rewrite it as a current:
   \[
   \partial_Y W=-L_AJ^A.
   \]
5. Divide by \(W\) to obtain the generalized velocity:
   \[
   v^A=J^A/W.
   \]
6. Identify whether the term needs only \(s_A\), or also \(L_As_B\).

---

## 3. Conventions to Check

These points must be checked against the exact NLO Hamiltonian convention before implementation.

### 3.1 Observable-side sign

The NLO JIMWLK paper writes

\[
\frac{d}{dY}{\cal O} = - H^{\rm JIMWLK}{\cal O}.
\]

For a density \(W\), this means

\[
\partial_Y W = -H^\dagger W.
\]

This global sign must be aligned with the LO convention used in the deterministic-flow paper.

### 3.2 Left/right generator relation

The paper gives

\[
J_L^a(x)=S_A^{ab}(x)J_R^b(x),
\]

\[
J_R^a(x)=S_A^{ba}(x)J_L^b(x).
\]

We use the left basis

\[
L_x^a\equiv J_L^a(x).
\]

Thus

\[
J_R^a(x)=S_x^{ba}L_x^b.
\]

### 3.3 Operator ordering

The NLO Hamiltonian states that all factors of \(J\) are ordered to the right of all factors of \(S\). Thus, in the Hamiltonian expression, the \(J\)'s do not act on the \(S\)'s appearing as coefficients.

This is crucial. When converting to the density side, derivatives will act on coefficient functions after integration by parts, and those derivative-of-coefficient terms become part of the current.

### 3.4 Singlet vs nonsinglet Hamiltonian

For dipole/NLO BK validation, the singlet Hamiltonian is sufficient.

For a configuration-level probability current acting on Wilson-line configurations, the nonsinglet/generalized Hamiltonian is safer. This means using barred kernels for the two-generator terms:

\[
K_{JSJ}\rightarrow \bar K_{JSJ},
\]

\[
K_{JSSJ}\rightarrow \bar K_{JSSJ},
\]

\[
K_{q\bar q}\rightarrow \bar K_{q\bar q}.
\]

The kernels \(K_{JJSJ}\) and \(K_{JJSSJ}\) remain unchanged.

---

## 4. Generic Current Formula for Up to Third Order

Assume the density-side operator can be written in normal form

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12 L_A L_B(K_2^{AB}W)
-
\frac16 L_A L_B L_C(K_3^{ABC}W).
\]

Then

\[
\partial_YW=-L_AJ^A,
\]

with

\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

Thus

\[
v^A=\frac{J^A}{W}
\]

is

\[
v^A
=
K_1^A
-
\frac12\frac{1}{W}L_B(K_2^{AB}W)
+
\frac16\frac{1}{W}L_BL_C(K_3^{ABC}W).
\]

Using

\[
s_A=L_A\log W,
\qquad
L_AW=Ws_A,
\]

the second-order contribution is

\[
\frac{1}{W}L_B(K_2^{AB}W)
=
L_BK_2^{AB}+K_2^{AB}s_B.
\]

The third-order contribution is

\[
\frac{1}{W}L_BL_C(K_3^{ABC}W)
=
L_BL_CK_3^{ABC}
+
(L_CK_3^{ABC})s_B
+
(L_BK_3^{ABC})s_C
+
K_3^{ABC}(L_Bs_C+s_Bs_C).
\]

Therefore the generalized velocity is

\[
\boxed{
v^A_{\rm gen}
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
K_3^{ABC}
\left(
L_Bs_C+s_Bs_C
\right)
\right].
}
\]

If \(K_3=0\), this reduces to the ordinary Itô/Fokker--Planck probability-flow formula. For LO JIMWLK,

\[
K_2^{AB}=2\chi^{AB},
\qquad
K_1^A=L_B\chi^{AB},
\]

and one obtains

\[
v^A_{\rm LO}=-\chi^{AB}s_B.
\]

---

# 5. First NLO Term: \(K_{JSJ}\)

## 5.1 Observable-side Hamiltonian term

The first NLO term has the same charge-generator structure as LO:

\[
H_{JSJ}
=
\int_{x,y,z}
K_{JSJ}(x,y;z)
\left[
J_L^a(x)J_L^a(y)
+
J_R^a(x)J_R^a(y)
-
2J_L^a(x)S_z^{ab}J_R^b(y)
\right].
\]

Here

\[
S_z^{ab}\equiv S_A^{ab}(z).
\]

For a configuration-level current one should replace

\[
K_{JSJ}\rightarrow \bar K_{JSJ}.
\]

The derivation below is identical with \(K\) or \(\bar K\).

---

## 5.2 Convert \(J_R\) to the left basis

Use

\[
J_R^a(x)=S_x^{ba}L_x^b.
\]

### First term

\[
J_L^a(x)J_L^a(y)
=
L_x^aL_y^a
=
\delta^{bc}L_x^bL_y^c.
\]

### Second term

\[
J_R^a(x)J_R^a(y)
=
S_x^{ba}L_x^bS_y^{ca}L_y^c.
\]

Using the Hamiltonian ordering, \(J\)'s do not act on the \(S\)'s in the coefficient, so

\[
J_R^a(x)J_R^a(y)
=
S_x^{bd}S_y^{cd}L_x^bL_y^c.
\]

### Third term

\[
-2J_L^a(x)S_z^{ab}J_R^b(y)
=
-2L_x^aS_z^{ab}S_y^{cb}L_y^c.
\]

Thus

\[
-2J_L^a(x)S_z^{ab}J_R^b(y)
=
-2S_z^{bd}S_y^{cd}L_x^bL_y^c,
\]

after relabeling dummy color indices.

---

## 5.3 Raw left-basis coefficient

The raw coefficient multiplying \(L_x^bL_y^c\) is

\[
M_{\rm raw}^{bc}(x,y;z)
=
\delta^{bc}
+
S_x^{bd}S_y^{cd}
-
2S_z^{bd}S_y^{cd}.
\]

Therefore

\[
H_{JSJ}
=
\int_{x,y,z}
K_{JSJ}(x,y;z)
M_{\rm raw}^{bc}(x,y;z)
L_x^bL_y^c.
\]

---

## 5.4 Symmetrization in \((x,b)\leftrightarrow(y,c)\)

Since \(K_{JSJ}(x,y;z)\) is symmetric under \(x\leftrightarrow y\), the second-order operator can be represented using the symmetrized coefficient

\[
M_{\rm sym}^{bc}(x,y;z)
=
\frac12
\left[
M_{\rm raw}^{bc}(x,y;z)
+
M_{\rm raw}^{cb}(y,x;z)
\right].
\]

This gives

\[
M_{\rm sym}^{bc}
=
\delta^{bc}
+
S_x^{bd}S_y^{cd}
-
S_z^{bd}S_y^{cd}
-
S_x^{bd}S_z^{cd}.
\]

Using orthogonality of adjoint Wilson lines,

\[
S_z^{bd}S_z^{cd}=\delta^{bc},
\]

we can write

\[
M_{\rm sym}^{bc}
=
S_z^{bd}S_z^{cd}
+
S_x^{bd}S_y^{cd}
-
S_z^{bd}S_y^{cd}
-
S_x^{bd}S_z^{cd}.
\]

Thus

\[
\boxed{
M_{\rm sym}^{bc}(x,y;z)
=
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
}
\]

This is the key result: the \(K_{JSJ}\) term has the same square-difference color structure as the LO kernel.

---

## 5.5 Define the density-side kernel contribution

Define the combined indices

\[
A=(x,b),
\qquad
B=(y,c).
\]

Then define

\[
\boxed{
\chi_{JSJ}^{(x,b)(y,c)}[U]
=
\int_z
K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
}
\]

For the nonsinglet/generalized Hamiltonian, use

\[
\boxed{
\chi_{JSJ}^{(x,b)(y,c)}[U]
=
\int_z
\bar K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
}
\]

There may be an overall sign depending on the convention relating observable-side evolution

\[
\frac{d}{dY}O=-HO
\]

to density-side evolution. This sign must be fixed by matching the LO limit to the existing LO divergence-form convention.

---

## 5.6 Current and velocity contribution

If the density-side contribution is aligned to the divergence-form convention,

\[
(\partial_YW)_{JSJ}
=
L_A(\chi_{JSJ}^{AB}L_BW),
\]

then

\[
L_BW=Ws_B,
\]

so

\[
(\partial_YW)_{JSJ}
=
L_A(\chi_{JSJ}^{AB}s_BW)
=
-L_A(v_{JSJ}^AW).
\]

Hence

\[
\boxed{
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
}
\]

Equivalently,

\[
\boxed{
J_{JSJ}^A
=
-\chi_{JSJ}^{AB}s_BW.
}
\]

---

## 5.7 Need for Hessian-score?

The \(K_{JSJ}\) term contains only two charge generators. Therefore it generates a second-order density operator and only requires the ordinary score,

\[
s_A=L_A\log W.
\]

It does not require

\[
L_As_B.
\]

Thus

\[
\boxed{
K_{JSJ}\text{ is an LO-like NLO correction to the score current.}
}
\]

---

## 5.8 Implementation form

The velocity action can be written as

\[
(v_{JSJ})_x^b
=
-
\sum_{y,c}
\chi_{JSJ}^{(x,b)(y,c)}s_y^c.
\]

Using the explicit kernel,

\[
(v_{JSJ})_x^b
=
-
\sum_{y,z,c,d}
K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right)
s_y^c.
\]

Equivalently,

\[
\boxed{
(v_{JSJ})_x^b
=
-
\sum_{z,d}
\left(S_x^{bd}-S_z^{bd}\right)
\sum_{y,c}
K_{JSJ}(x,y;z)
\left(S_y^{cd}-S_z^{cd}\right)
s_y^c.
}
\]

This is LO-like in structure, but \(K_{JSJ}(x,y;z)\) may not factorize as a simple square-root convolution. For \(12^2\), direct summation is feasible for a first test.

---

# 6. Checkpoints Before Coding

Before implementing the \(K_{JSJ}\) contribution, verify:

1. **Overall sign**  
   Replace \(K_{JSJ}\) by the LO kernel and confirm that the density-side evolution reproduces the existing LO convention
   \[
   \partial_YW=L_A(\chi^{AB}L_BW).
   \]

2. **Barred vs unbarred kernel**  
   Use \(\bar K_{JSJ}\) for configuration-level flow; use \(K_{JSJ}\) only for singlet observable validation.

3. **Left/right convention**  
   Confirm the code uses
   \[
   J_R^a(x)=S_x^{ba}J_L^b(x).
   \]

4. **Kernel symmetry**  
   Confirm
   \[
   K_{JSJ}(x,y;z)=K_{JSJ}(y,x;z),
   \]
   so that the symmetrized square-difference form is valid.

5. **Adjoint Wilson line convention**  
   Confirm whether the code uses \(S_A^{ab}\) or \(\mathrm{Ad}_U^{ab}\), and whether it matches the paper's \(J_L,J_R\) convention.

---

# 7. Next Term to Derive

The next term is

\[
K_{JSSJ}(x,y;z,z')
\left[
f^{abc}f^{def}
J_L^a(x)S_z^{be}S_{z'}^{cf}J_R^d(y)
-
N_cJ_L^a(x)S_z^{ab}J_R^b(y)
\right].
\]

It is also second order in charge generators. Therefore it should require only the score \(s_A\), not the Hessian-score \(L_As_B\). However, unlike \(K_{JSJ}\), it is not manifestly a positive square-root diffusion. We should derive its left-basis coefficient and then decide whether it is best represented as:

1. a divergence-form second-order current, or
2. a generic second-order Itô-form contribution with kernel-derivative drift.

