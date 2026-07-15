# Analytic Physical Coefficient Derivative Derivation

## Scope

This document records the currently implemented analytic/local coefficient
derivative backend. It is non-production and keeps the finite-difference backend
as the reference oracle.

The implemented analytic scope includes the two-generator \(dK2\)
contractions and the \(K_{JJSJ}\) cubic derivative contractions:

\[
dK2^A=L_BK_2^{AB}.
\]

\[
(LC\_K3)^{AB}=L_CK_3^{ABC},
\qquad
(LB\_K3)^{AC}=L_BK_3^{ABC},
\qquad
d2K3^A=L_BL_CK_3^{ABC}.
\]

\(K_{JJSSJ}\) remains pending for a true analytic backend.

## Local Generator Rules

The code convention is

\[
L_x^aF(U)=\left.\frac{d}{d\epsilon}
F(e^{i\epsilon t^a}U_x)\right|_{\epsilon=0}.
\]

Thus

\[
L_x^aU_y=i\delta_{xy}t^aU_y,
\qquad
L_x^aU_y^\dagger=-i\delta_{xy}U_y^\dagger t^a.
\]

For the adjoint Wilson line,

\[
S_A^{ab}=2\,\mathrm{Re\,Tr}(t^aUt^bU^\dagger),
\]

the finite-difference calibrated rule is

\[
L^hS_A^{ab}=f^{hac}S_A^{cb}.
\]

Ordered second derivatives are preserved:

\[
L^gL^hS_A^{ab}=f^{hac}f^{gcd}S_A^{db}.
\]

Same-site commutators use the project convention

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

## Implemented Two-Generator Sectors

### KJSJ

The skeleton coefficient is

\[
K_2^{(x,b)(y,c)}
=
-2\sum_z K_{JSJ}(x,y;z)
(S_x-S_z)^{bd}(S_y-S_z)^{cd}.
\]

The analytic derivative applies the product rule only to the two adjoint
factors because the physical coordinate kernel is Wilson-line independent.

### KJSSJ

The ordered coefficient is built as

\[
A^{ab}(x,y)=
\sum_{z,z'}K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}S_z^{de}(S_{z'}^{cf}-S_z^{cf}),
\]

then converted to the left basis:

\[
K_2^{(x,a)(y,h)}=A^{ab}(x,y)S_y^{hb}.
\]

The analytic derivative differentiates \(S_z\), \(S_{z'}-S_z\), and the final
right-to-left conversion factor \(S_y\), preserving the \(z,z'\) ordering used
by `two_generator_terms.py`.

### Kqbarq

The full coefficient is decomposed into trace and subtraction pieces:

\[
A^{ab}(x,y)=
\sum_{z,z'}K_{q\bar q}(x,y;z,z')
\left[
2\mathrm{Tr}(U_z^\dagger t^aU_{z'}t^b)-S_z^{ab}
\right].
\]

The implemented functions are:

- `analytic_dK2_Kqbarq_trace(...)`;
- `analytic_dK2_Kqbarq_subtraction(...)`;
- `analytic_dK2_Kqbarq(...)`.

The trace derivative uses the calibrated fundamental rules for \(U\) and
\(U^\dagger\). The subtraction derivative uses the calibrated adjoint rule.
The final right-to-left conversion \(S_y^{hb}\) is also differentiated.

## KJJSJ Cubic Derivatives

The physical adapter provides KLM-normalized cubic coefficients:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

The \(K_{JJSJ}\) analytic derivatives act on the normalized coefficient used
by the skeleton; the \((-i)\) factor is not re-applied in derivative code.
They differentiate the same right-to-left converted and canonicalized
normal-form tensors used by the dense skeleton, including lower-order
commutator corrections.

The implemented coefficient blocks are:

\[
A_{LLR}^{dea}(x,y,w)=
\int_zK_{JJSJ}(w;x,y;z)f^{bde}S_z^{ba},
\]

\[
B_{LRR}^{ade}(w,x,y)=
-\int_zK_{JJSJ}(w;x,y;z)f^{bde}S_z^{ab},
\]

\[
V^{deb}(x,y,w)=
\frac13\int_zK_{JJSJ}(w;x,y;z)f^{bde}.
\]

The first derivatives use the calibrated adjoint rule
\(L^hS_A^{ab}=f^{hac}S_A^{cb}\). The ordered second derivatives use
\(L^gL^hS_A^{ab}=f^{hac}f^{gcd}S_A^{db}\), preserving same-site ordering.

Validation against the finite-difference oracle is recorded in
`reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`. In the
nonzero synthetic scan, full first-derivative residuals are at the
\(10^{-11}\) level and full ordered \(d2K3\) residuals are at the \(10^{-8}\)
level for the stable second-difference window.

`backend="analytic"` now supports `KJJSJ` without global FD fallback and still
raises for pending `KJJSSJ`. The explicit `backend="hybrid_local_fd"` path
remains available and labeled.

## Source-Code Mapping

- Primitive rules: `src/nlo_current/analytic_lie_derivatives.py`
- Two-generator derivatives:
  `src/nlo_current/analytic_two_generator_derivatives.py`
- KJJSJ cubic derivatives:
  `src/nlo_current/analytic_cubic_derivatives.py`
- Structured backend:
  `src/nlo_current/physical_coefficient_derivatives.py`
- FD oracle: `src/nlo_current/coefficient_derivatives.py`
- Velocity wrapper integration: `src/nlo_current/physical_nlo_current.py`
