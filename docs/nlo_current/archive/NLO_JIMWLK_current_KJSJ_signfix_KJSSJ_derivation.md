# NLO JIMWLK Generalized Probability Current  
## Updated Worklog: \(K_{JSJ}\) Sign Fix and \(K_{JSSJ}\) Derivation

This note completes two updates:

1. Fix the overall sign in the first-term \(K_{JSJ}\) current using the KLM convention
   \[
   \frac{d}{dY}{\cal O}=-H^{\rm JIMWLK}{\cal O}.
   \]
2. Continue the second-term \(K_{JSSJ}\) derivation using the identity in Appendix B of Kovner--Lublinsky--Mulian.

The main result is:

\[
\boxed{
K_{JSJ}\text{ is LO-like and contributes }v^A=-\chi_{JSJ}^{AB}s_B.
}
\]

\[
\boxed{
K_{JSSJ}\text{ is still second order in charge generators, so it needs }s_A
\text{ but not }L_As_B.
}
\]

However, unlike \(K_{JSJ}\), the \(K_{JSSJ}\) term is not manifestly a positive square-difference diffusion. Its safe representation is a generic second-order current with a symmetric diffusion-like part and a possible antisymmetric commutator drift.

---

# 1. Reference Conventions from KLM

The NLO JIMWLK paper uses

\[
\frac{d}{dY}{\cal O}=-H^{\rm JIMWLK}{\cal O}.
\]

Therefore, for the probability density \(W[U]\),

\[
\partial_YW=-H^\dagger W.
\]

The left and right generators satisfy

\[
J_L^a(x)=S_A^{ab}(x)J_R^b(x),
\]

\[
J_R^a(x)=S_A^{ba}(x)J_L^b(x).
\]

We use the left-basis notation

\[
L_x^a\equiv J_L^a(x).
\]

Thus

\[
J_R^a(x)=S_x^{ba}L_x^b.
\]

The Hamiltonian ordering convention is also important: the Wilson-line factors \(S\) appearing explicitly in the Hamiltonian are coefficient functions, and the charge generators are ordered so as not to act on those explicit \(S\)'s inside the Hamiltonian expression. When adjointing to the density side, derivatives do act on the resulting coefficient functions, producing current/drift terms.

---

# 2. Checkpoint A: Overall Sign

Suppose a two-generator term in the observable-side Hamiltonian is schematically

\[
H_2=C_H^{AB}[U]L_A L_B.
\]

Because KLM uses

\[
\partial_Y{\cal O}=-H{\cal O},
\]

the density-side evolution is

\[
\partial_YW=-H_2^\dagger W.
\]

For a symmetric two-derivative coefficient, this means the diffusion-like density tensor is

\[
\boxed{D^{AB}=-2C_H^{AB}.}
\]

Equivalently, if we write a divergence-form kernel as

\[
\partial_YW=L_A(\chi^{AB}L_BW),
\]

then the observable-side second-order coefficient is

\[
C_H^{AB}=-\chi^{AB}.
\]

Therefore,

\[
\boxed{\chi^{AB}=-C_H^{AB}.}
\]

This is also consistent with the LO check in KLM Appendix A: replacing

\[
K_{JSJ}\rightarrow -\frac12M
\]

recovers the LO dipole evolution. Thus the first \(K_{JSJ}\) current must include an overall minus sign in the definition of \(\chi_{JSJ}\).

---

# 3. Corrected First-Term Result: \(K_{JSJ}\)

## 3.1 Observable-side structure

The \(K_{JSJ}\) term is

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

After converting to the left basis and symmetrizing in \((x,b)\leftrightarrow(y,c)\), the color structure is

\[
M_{JSJ}^{bc}(x,y;z)
=
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
\]

Thus the observable-side coefficient is

\[
C_{JSJ}^{(x,b)(y,c)}
=
\int_z
K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
\]

Because \(\chi=-C_H\), the density-side divergence kernel is

\[
\boxed{
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
}
\]

For configuration-level evolution, use the nonsinglet/generalized barred kernel:

\[
\boxed{
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
}
\]

## 3.2 LO check

In the LO singlet limit,

\[
K_{JSJ}\rightarrow -\frac12M(x,y;z).
\]

Then

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
\frac12
\int_z
M(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right),
\]

which has the expected positive LO square-difference form.

Thus the sign is fixed.

## 3.3 Current contribution

Since the term is LO-like,

\[
(\partial_YW)_{JSJ}
=
L_A(\chi_{JSJ}^{AB}L_BW).
\]

Using

\[
L_BW=Ws_B,
\]

we obtain

\[
(\partial_YW)_{JSJ}
=
L_A(\chi_{JSJ}^{AB}s_BW)
=
-L_A(v_{JSJ}^AW).
\]

Hence

\[
\boxed{v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.}
\]

Equivalently,

\[
\boxed{J_{JSJ}^A=-\chi_{JSJ}^{AB}s_BW.}
\]

This term only requires the ordinary score \(s_A\). It does not require \(L_As_B\).

---

# 4. Checkpoint C: Barred vs Unbarred Kernels

For color-singlet observables such as the dipole, the unbarred singlet Hamiltonian is sufficient.

For configuration-level probability current, the safer choice is the nonsinglet/generalized Hamiltonian. Therefore use

\[
K_{JSJ}\rightarrow \bar K_{JSJ},
\]

\[
K_{JSSJ}\rightarrow \bar K_{JSSJ},
\]

\[
K_{q\bar q}\rightarrow \bar K_{q\bar q}.
\]

The three-generator kernels are unchanged:

\[
K_{JJSJ}\rightarrow K_{JJSJ},
\qquad
K_{JJSSJ}\rightarrow K_{JJSSJ}.
\]

Thus, for a genuine Wilson-line configuration flow, the kernel set should be

\[
\boxed{
\bar K_{JSJ},
\quad
\bar K_{JSSJ},
\quad
\bar K_{q\bar q},
\quad
K_{JJSJ},
\quad
K_{JJSSJ}.
}
\]

---

# 5. Second Term: \(K_{JSSJ}\)

## 5.1 Starting point

The second NLO term is

\[
H_{JSSJ}
=
\int_{x,y,z,z'}
K_{JSSJ}(x,y;z,z')
\left[
f^{abc}f^{def}
J_L^a(x)S_z^{be}S_{z'}^{cf}J_R^d(y)
-
N_cJ_L^a(x)S_z^{ab}J_R^b(y)
\right].
\]

For configuration-level current, replace

\[
K_{JSSJ}\rightarrow \bar K_{JSSJ}.
\]

In the rest of this derivation we write \(K_{JSSJ}\), with the understanding that the barred version should be used for nonsinglet/configuration-level evolution.

---

## 5.2 Use the KLM Appendix B identity

KLM Appendix B gives the identity

\[
f^{abc}f^{def}J_L^a(x)S_z^{be}S_v^{cf}J_R^d(y)
-
N_cJ_L^a(x)S_z^{ab}J_R^b(y)
\]

\[
=
f^{adc}f^{bef}
\left[
S_z^{de}S_v^{cf}
-
S_z^{de}S_z^{cf}
\right]
J_L^a(x)J_R^b(y).
\]

Setting \(v=z'\), this becomes

\[
\boxed{
f^{abc}f^{def}J_L^a(x)S_z^{be}S_{z'}^{cf}J_R^d(y)
-
N_cJ_L^a(x)S_z^{ab}J_R^b(y)
=
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
J_L^a(x)J_R^b(y).
}
\]

Thus

\[
H_{JSSJ}
=
\int_{x,y,z,z'}
K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
J_L^a(x)J_R^b(y).
\]

This identity is important because it makes the subtraction explicit:

\[
z'=z
\quad\Longrightarrow\quad
S_{z'}^{cf}-S_z^{cf}=0.
\]

So the virtual subtraction cancels the coincident \(z'=z\) structure.

---

## 5.3 Convert \(J_R\) to the left basis

Use

\[
J_R^b(y)=S_y^{hb}L_y^h.
\]

Then

\[
J_L^a(x)J_R^b(y)
=
L_x^a\left(S_y^{hb}L_y^h\right).
\]

At the structural coefficient level, the two-derivative part is

\[
S_y^{hb}L_x^aL_y^h.
\]

Thus the observable-side two-derivative coefficient is

\[
\boxed{
C_{JSSJ}^{(x,a)(y,h)}[U]
=
\int_{z,z'}
K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
S_y^{hb}.
}
\]

Equivalently, using a barred kernel for configuration-level flow,

\[
\boxed{
C_{JSSJ}^{(x,a)(y,h)}[U]
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
S_y^{hb}.
}
\]

This is the central coefficient for the \(K_{JSSJ}\) current.

---

# 6. Symmetric and Antisymmetric Decomposition

Let

\[
A=(x,a),
\qquad
B=(y,h).
\]

The coefficient need not be manifestly symmetric:

\[
C_{JSSJ}^{AB}\neq C_{JSSJ}^{BA}.
\]

Therefore decompose

\[
C_{JSSJ}^{AB}
=
C_{JSSJ}^{(AB)}
+
C_{JSSJ}^{[AB]},
\]

where

\[
C_{JSSJ}^{(AB)}
=
\frac12
\left(
C_{JSSJ}^{AB}
+
C_{JSSJ}^{BA}
\right),
\]

\[
C_{JSSJ}^{[AB]}
=
\frac12
\left(
C_{JSSJ}^{AB}
-
C_{JSSJ}^{BA}
\right).
\]

The symmetric part gives a diffusion-like density tensor

\[
\boxed{D_{JSSJ}^{AB}=-2C_{JSSJ}^{(AB)}.}
\]

The antisymmetric part contributes through the Lie-algebra commutator. If

\[
[L_A,L_B]=f_{AB}^{\ \ C}L_C,
\]

then the antisymmetric part gives a first-order drift-like contribution

\[
\boxed{
\mu_{\rm asym}^{C}
=
\frac12 f_{AB}^{\ \ C}P_{JSSJ}^{[AB]},
}
\]

where \(P\) is the coefficient in the density-side ordering. If one writes the density-side second-derivative term as

\[
-L_A L_B(P^{AB}W),
\]

then \(P^{AB}=C^{BA}\). Equivalently, one may do the decomposition directly after placing the density-side operator into a fixed normal order.

This is a bookkeeping point: before coding, one should choose one canonical normal ordering and keep it throughout.

---

# 7. Density-Side Current for the \(K_{JSSJ}\) Term

In a generic second-order density normal form,

\[
(\partial_YW)_{JSSJ}
=
-L_A(\mu_{JSSJ}^AW)
+
\frac12L_AL_B(D_{JSSJ}^{AB}W),
\]

the current is

\[
J_{JSSJ}^A
=
\mu_{JSSJ}^AW
-
\frac12L_B(D_{JSSJ}^{AB}W).
\]

Dividing by \(W\),

\[
\boxed{
v_{JSSJ}^A
=
\mu_{JSSJ}^A
-
\frac12L_BD_{JSSJ}^{AB}
-
\frac12D_{JSSJ}^{AB}s_B.
}
\]

Using

\[
D_{JSSJ}^{AB}=-2C_{JSSJ}^{(AB)},
\]

the symmetric score part can be written as

\[
-\frac12D_{JSSJ}^{AB}s_B
=
C_{JSSJ}^{(AB)}s_B.
\]

Therefore

\[
\boxed{v_{JSSJ,{\rm score}}^A=C_{JSSJ}^{(AB)}s_B.}
\]

The remaining drift-like part is

\[
\boxed{
v_{JSSJ,{\rm drift}}^A
=
\mu_{JSSJ}^A
-
\frac12L_BD_{JSSJ}^{AB}.
}
\]

This includes derivative-of-kernel/coefficient terms and the commutator drift generated by the antisymmetric coefficient.

---

# 8. Does \(K_{JSSJ}\) Need Hessian-Score?

No.

The \(K_{JSSJ}\) term contains two charge generators. Therefore its density-side operator is at most second order. A second-order term requires only

\[
s_A=L_A\log W,
\]

not

\[
L_As_B=L_AL_B\log W.
\]

Thus

\[
\boxed{
K_{JSSJ}\text{ requires score }s_A\text{ but not Hessian-score }L_As_B.
}
\]

Hessian-score terms first appear in the three-generator pieces:

\[
K_{JJSJ},
\qquad
K_{JJSSJ},
\qquad
K_{JJJ}\text{ or the }1/3\text{ virtual replacements}.
\]

---

# 9. Exact Current Form for the Explicit Two-Derivative Piece

If we isolate the density-side two-derivative term as

\[
(\partial_YW)_2
=
-L_A L_B(P^{AB}W),
\]

then it is already a total divergence:

\[
(\partial_YW)_2
=
-L_AJ^A,
\]

with

\[
\boxed{J^A=L_B(P^{AB}W).}
\]

Therefore

\[
\boxed{
v^A
=
\frac{1}{W}L_B(P^{AB}W)
=
L_BP^{AB}+P^{AB}s_B.
}
\]

This expression is exact for the chosen density-side normal ordering.

However, for comparison with ordinary Fokker--Planck notation, it is often better to decompose \(P\) into symmetric and antisymmetric parts and write the result in the form

\[
v^A
=
\mu^A
-
\frac12L_BD^{AB}
-
\frac12D^{AB}s_B.
\]

Both forms are equivalent once the same normal ordering is used.

---

# 10. Current Status of Checkpoints A--C

## Checkpoint A: Overall sign

Resolved.

Because KLM uses

\[
d{\cal O}/dY=-H{\cal O},
\]

the density-side two-derivative kernel satisfies

\[
D=-2C_H,
\qquad
\chi=-C_H.
\]

Therefore the corrected \(K_{JSJ}\) kernel is

\[
\chi_{JSJ}=-\int K_{JSJ}(S_x-S_z)(S_y-S_z).
\]

## Checkpoint B: Symmetry of \(K_{JSSJ}\)

Partially resolved.

KLM states that \(K_{JSSJ}\) is symmetric under relevant exchanges such as \(x\leftrightarrow y\) and/or \(z\leftrightarrow z'\), depending on convention, and Appendix B gives the identity that rewrites the subtraction as

\[
S_{z'}-S_z.
\]

This confirms the subtraction structure and the \(z'=z\) cancellation.

However, the full coefficient

\[
C_{JSSJ}^{AB}
\]

is not manifestly symmetric after converting \(J_R\) to the left basis. Therefore one should keep

\[
C^{(AB)}
\quad\text{and}\quad
C^{[AB]}
\]

separate until a direct color/kernel symmetry check proves further simplification.

## Checkpoint C: Barred kernels

Resolved.

For configuration-level NLO probability current use

\[
\bar K_{JSJ},
\qquad
\bar K_{JSSJ},
\qquad
\bar K_{q\bar q},
\qquad
K_{JJSJ},
\qquad
K_{JJSSJ}.
\]

For dipole/NLO BK validation, unbarred singlet kernels are sufficient.

---

# 11. Practical Next Steps

The next checks should be:

1. **Numerical symmetry check for \(C_{JSSJ}^{AB}\)**  
   Construct \(C^{AB}\) on a small lattice and measure
   \[
   \frac{\|C^{AB}-C^{BA}\|}{\|C^{AB}\|}.
   \]

2. **Commutator drift check**  
   If \(C^{[AB]}\neq0\), compute
   \[
   \mu_{\rm asym}^{C}
   =
   \frac12 f_{AB}^{\ \ C}P^{[AB]}.
   \]

3. **Dipole validation**  
   Apply the \(K_{JSJ}\) and \(K_{JSSJ}\) pieces to a dipole and verify agreement with the Appendix A expressions.

4. **Configuration-level choice**  
   Use barred kernels for the actual probability-current flow, but keep unbarred kernels for comparison to the singlet dipole formulas.

---

# 12. Summary

The first term is now fixed:

\[
\boxed{
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right).
}
\]

\[
\boxed{v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.}
\]

The second term has central coefficient

\[
\boxed{
C_{JSSJ}^{(x,a)(y,h)}
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
S_y^{hb}.
}
\]

It contributes a generic second-order current,

\[
\boxed{
v_{JSSJ}^A
=
\mu_{JSSJ}^A
-
\frac12L_BD_{JSSJ}^{AB}
-
\frac12D_{JSSJ}^{AB}s_B,
}
\]

with

\[
\boxed{D_{JSSJ}^{AB}=-2C_{JSSJ}^{(AB)}.}
\]

It needs \(s_A\), but not \(L_As_B\).

