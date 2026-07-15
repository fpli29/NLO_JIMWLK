# NLO JIMWLK Generalized Probability Current: Derivation and Validation Summary

This document consolidates the dense small-lattice, non-production derivations
and validations completed under `docs/nlo_current/`, `reports/nlo_current/`,
and `src/nlo_current/`. It is a theory and implementation-orientation summary,
not a production evolution design.

No physical coordinate kernels are implemented here. No production evolution
code and no score or Hessian-score model training are included.

## 1. Scope And Conventions

Observable-side evolution is written as

\[
\frac{d}{dY}{\cal O}=-H{\cal O}.
\]

The density-side equation is

\[
\partial_Y W=-H^\dagger W.
\]

The local left/right generator convention used by the diagnostic code is

\[
L_x^aF(U)=\frac{d}{d\epsilon}F(e^{i\epsilon t^a}U_x)\bigg|_{\epsilon=0},
\qquad
R_x^aF(U)=\frac{d}{d\epsilon}F(U_xe^{i\epsilon t^a})\bigg|_{\epsilon=0}.
\]

The adjoint Wilson line connects the two bases:

\[
J_R^a(x)=S_A^{ba}(x)J_L^b(x),
\qquad
J_L^a(x)=S_A^{ab}(x)J_R^b(x).
\]

The dipole observable used for Appendix A validation is

\[
s(u,v)=\frac1{N_c}\mathrm{tr}[U^\dagger(u)U(v)].
\]

Notation mapping: \(U\) is the fundamental Wilson line used in the diagnostic
code. KLM often writes this same fundamental Wilson line as \(S\). In this
summary, \(S_A\) denotes the adjoint Wilson line; shorthand expressions such as
\(J_LSJ_R\) should be read with the representation implied by the surrounding
Hamiltonian block.

Kernel notation is also context-dependent. Configuration-level current formulas
may use barred kernels, while the Appendix A dipole validation uses the
unbarred singlet kernels appearing in the KLM dipole-action formulas.

Appendix A targets are calibrated as Hamiltonian actions:

\[
\text{TeX Appendix A target}=H_{\rm sector}s.
\]

The raw cubic direct-action functions use Hermitian generator finite
differences. For the tested one-\(f\) cubic sectors,

\[
\text{TeX target}=(-i)\times\text{raw direct action}.
\]

This convention is exposed explicitly by
`klm_normalized_cubic_direct_action(raw_direct) = (-i) * raw_direct`; the raw
direct-action functions are not redefined.

Primary sources used:

- `references/WORKNLO.tex`
- `docs/nlo_current/KLM_appendix_A_dipole_targets_notes.md`
- `src/nlo_current/su3_adjoint.py`
- `src/nlo_current/dipole_observable.py`
- `src/nlo_current/dipole_hamiltonian_action.py`
- `src/nlo_current/dipole_appendix_targets.py`

## 2. LO Reference Point

The LO density equation can be written directly in divergence form:

\[
\partial_Y W=L_A(\chi^{AB}L_BW).
\]

With score

\[
s_A=L_A\log W,
\]

one has

\[
L_BW=Ws_B.
\]

Therefore

\[
\partial_YW=L_A(\chi^{AB}s_BW)=-L_A(v^AW),
\]

with

\[
v^A=-\chi^{AB}s_B.
\]

Here \(\chi\) is in the project's divergence-form normalization. In a
\(\frac12D\) Fokker-Planck convention, \(D=2\chi\).

For a generic Itô-style second-order Fokker-Planck equation,

\[
\partial_YW=-L_A(b^AW)+\frac12L_AL_B(D^{AB}W),
\]

the current is

\[
J^A=b^AW-\frac12L_B(D^{AB}W),
\]

and the velocity contains both a coefficient-derivative drift and a score term:

\[
v^A=b^A-\frac12\left[L_BD^{AB}+D^{AB}s_B\right].
\]

LO JIMWLK is special because the operator is already in the divergence form
\(L_A(\chi^{AB}L_BW)\). In that form, the derivative-of-\(\chi\) drift is not a
separate velocity term; it cancels into the divergence representation, leaving
the score current \(v^A=-\chi^{AB}s_B\).

## 3. General Density Normal Form

The working dense normal form is

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

Dividing by \(W\), the velocity is

\[
v^A
=
K_1^A
-\frac12
\left[
L_BK_2^{AB}+K_2^{AB}s_B
\right]
+\frac16
\left[
L_BL_CK_3^{ABC}
+(L_CK_3^{ABC})s_B
+(L_BK_3^{ABC})s_C
+K_3^{ABC}(H_{BC}+s_Bs_C)
\right],
\]

where

\[
H_{BC}=L_Bs_C=L_BL_C\log W.
\]

\(H_{BC}\) is an ordered Hessian-score. It is not assumed symmetric,
especially when \(B\) and \(C\) contain coincident-site Lie derivatives.

The key rule is:

\[
\text{k-th order generator terms require derivatives of }\log W\text{ up to order }k-1.
\]

Thus:

- two-generator terms require score;
- three-generator terms require score and Hessian-score.

## 4. Two-Generator Sector

### 4.1 \(K_{JSJ}\)

The Hamiltonian block is

\[
H_{JSJ}
=
\int K_{JSJ}
\left[
J_LJ_L+J_RJ_R-2J_LSJ_R
\right].
\]

The left-basis kernel result is

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
(S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}).
\]

The velocity is the LO-like score current

\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
\]

The LO sign check is

\[
K_{JSJ}\rightarrow -M/2,
\]

which gives the positive LO diffusion kernel in the corrected left-basis
normalization.

Appendix A validation passed with residual
`1.6172982426630268e-16`.

### 4.2 Ordered \(J_L A J_R\) Lemma

For an ordered block

\[
H_{LR}[A]=\int A^{ab}(x,y;U)J_L^a(x)J_R^b(y),
\]

use

\[
J_R^b(y)=S_y^{hb}L_y^h.
\]

The tested identity is

\[
L_y^hS_y^{hb}=0.
\]

The finite-difference check found
`max_abs_sum_h_Lh_Shb = 1.4745149545802860e-10`.

The density-side current is

\[
J_{LR}^{(y,h)}
=
S_y^{hb}L_x^a[A^{ab}W],
\]

and the velocity is

\[
v_{LR}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA^{ab}
+
A^{ab}s_x^a
\right].
\]

The finite-difference residuals in
`reports/nlo_current/ordered_lr_current_fd_report.md` decrease from
`9.9575953683810059e-07` at `eps=1e-4` to
`9.9512954300495782e-09` at `eps=1e-6`.

### 4.3 \(K_{JSSJ}\)

The ordered coefficient is

\[
A_{JSSJ}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}).
\]

The current is

\[
v_{JSSJ}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{JSSJ}^{ab}
+
A_{JSSJ}^{ab}s_x^a
\right].
\]

The coefficient symmetry diagnostic found

\[
\frac{\|C-C^T\|}{\|C\|}\sim O(1),
\]

with `max_r_asym = 1.3444397747997991e+00`. Therefore \(K_{JSSJ}\) is a generic
ordered-current / coefficient-drift term, not a symmetric score-only current.

The full Appendix A target was taken from the combined equation, excluding the
\(\widetilde K\) part. Appendix A validation passed with residual
`1.3877787807814457e-16`.

### 4.4 \(K_{q\bar q}\)

The ordered coefficient is

\[
A_{q\bar q}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{q\bar q}(x,y;z,z')
\left[
2\mathrm{tr}(U^\dagger(z)t^aU(z')t^b)
-
S_A^{ab}(z)
\right].
\]

The current is

\[
v_{q\bar q}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{q\bar q}^{ab}
+
A_{q\bar q}^{ab}s_x^a
\right].
\]

The \(z'=z\) subtraction identity passed:

\[
2\mathrm{tr}(U^\dagger t^aUt^b)-S_A^{ab}=0.
\]

The coefficient asymmetry is order one:
`max_r_asym = 1.0928913228548376e+00`. This sector also remains a generic
ordered-current / coefficient-drift term, not a symmetric score-only current.

Appendix A validation passed with residual `5.6325111571547335e-17`.

## 5. Three-Generator Sector: Distinct-Site Cubic Current

### 5.1 LLR Block

For

\[
H_{LLR}[A]=A^{dea}J_L^dJ_L^eJ_R^a,
\]

the tested current is

\[
v_{LLR}^{(w,h)}
=
-S_w^{ha}
\frac1W L_y^eL_x^d[A^{dea}W].
\]

Expanding gives

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

This is the first point where Hessian-score appears in the NLO current.

### 5.2 LRR And Virtual Blocks

The LRR current was tested by converting right generators to left-basis
generators while preserving the written ordering. The virtual blocks are the
canonical LLL/RRR cubic words with explicit \(1/3\) normalization.

The \(1/3\) virtual factor was tested in both \(K_{JJSJ}\) and \(K_{JJSSJ}\):
changing it changes the result by the expected factor and fails the validation.

### 5.3 \(K_{JJSJ}\)

The coefficient blocks are

\[
A_{LLR}^{dea}(x,y,w)
=
\int_z
K_{JJSJ}(w;x,y;z)
f^{bde}S_z^{ba},
\]

\[
B_{LRR}^{ade}(w,x,y)
=
-\int_z
K_{JJSJ}(w;x,y;z)
f^{bde}S_z^{ab},
\]

\[
V^{deb}(x,y,w)
=
\frac13\int_z
K_{JJSJ}(w;x,y;z)f^{bde}.
\]

Validation:

- LLR sign passed.
- LRR sign passed.
- virtual \(1/3\) passed.
- Hessian-score term nonzero.
- Appendix A target passed using the cubic \((-i)\) convention.

From `reports/nlo_current/kjjsj_cubic_requirements_report.md`:

- `max_sign_residual = 6.7563524631768235e-10`
- `max_hessian_score_component = 1.1358480138403666e-03`

Appendix A full residual: `1.1801832636420706e-15`.

### 5.4 \(K_{JJSSJ}\)

The coefficient blocks are

\[
A_{LLR}^{dea}(x,y,w)
=
\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}
S_z^{dc}S_{z'}^{eb},
\]

\[
B_{LRR}^{ade}(w,x,y)
=
-\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}
S_z^{cd}S_{z'}^{be},
\]

\[
V^{cba}(x,y,w)
=
\frac13\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')f^{acb}.
\]

The KLM-like simultaneous antisymmetry is

\[
K(w;x,y;z,z')=-K(w;y,x;z',z).
\]

Validation:

- LLR sign passed.
- LRR sign passed.
- virtual \(1/3\) passed.
- Hessian-score nonzero.
- Appendix A passed including \(\widetilde K\), pure eight-kernel
  contribution, and virtual term.

From `reports/nlo_current/kjjssj_cubic_requirements_report.md`:

- `max_sign_residual = 4.5529412514427881e-09`
- `max_hessian_score_component = 5.1891223835296681e-05`

Appendix A full residual: `1.9087382418229414e-15`.

## 6. Coincident-Site Commutators

Distinct-site validation is not enough because

\[
[L_x^a,L_y^b]=0\quad x\neq y,
\]

but at coincident sites

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

The canonical word ordering uses the combined index \(A=(x,a)\). The swap rule
is

\[
L_x^aL_x^b
=
L_x^bL_x^a
+
f^{abc}L_x^c.
\]

Canonicalization has the structure

\[
L_AL_BL_CF
=
\sum c\,L_{A'}L_{B'}L_{C'}F
+
\sum d\,L_{A'}L_{B'}F
+
\sum e\,L_{A'}F.
\]

Thus cubic sectors induce lower-order terms:

\[
K_3
\rightarrow
K_{3,\rm canonical},
\]

\[
K_2
\rightarrow
K_2+K_{2,\rm comm},
\]

\[
K_1
\rightarrow
K_1+K_{1,\rm comm}.
\]

Diagnostic results:

- \(K_{JJSJ}\): nonzero \(K_{2,\rm comm}\), zero \(K_{1,\rm comm}\) in the
  specific diagnostic.
- \(K_{JJSSJ}\): nonzero \(K_{2,\rm comm}\), nonzero \(K_{1,\rm comm}\).

The zero \(K_{1,\rm comm}\) for \(K_{JJSJ}\) is a diagnostic result, not a
theorem.

From `reports/nlo_current/cubic_commutator_corrections_report.md`:

| sector | quadratic commutator norm | linear commutator norm |
|---|---:|---:|
| \(K_{JJSJ}\) | `6.1541551762714866e+00` | `0.0000000000000000e+00` |
| \(K_{JJSSJ}\) | `1.0455495097967299e+01` | `3.6795897670206137e+00` |

## 7. NLO Current Skeleton

The non-production dense skeleton stores normal-form tensors as

```python
@dataclass
class NLOCurrentTerms:
    K1: np.ndarray
    K2: np.ndarray
    K3: np.ndarray
    metadata: dict
```

with

\[
D=8N_{\rm site},
\]

\[
K_1:(D,),\quad K_2:(D,D),\quad K_3:(D,D,D).
\]

Sector map:

| sector | normal-form contribution | density derivatives needed |
|---|---|---|
| \(K_{JSJ}\) | \(K_2\) | score |
| \(K_{JSSJ}\) | \(K_2\) | score + coefficient drift |
| \(K_{q\bar q}\) | \(K_2\) | score + coefficient drift |
| \(K_{JJSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |
| \(K_{JJSSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |

The non-production skeleton assembled all five sectors. Commutator corrections
are included and toggleable. Diagnostic velocity evaluation works with
supplied score and Hessian-score arrays.

From `reports/nlo_current/nlo_current_skeleton_demo_report.md`, the
derivative-enabled tiny-lattice section found:

- `dK2_norm = 1.2191975926577348e+00`
- `LC_K3_norm = 6.6584804802537545e-01`
- `LB_K3_norm = 6.8126547385360958e-01`
- `d2K3_norm = 3.0245993342733984e-01`
- `velocity_difference_norm = 7.1629857844995859e-01`

## 8. Coefficient Derivatives

Coefficient derivatives are needed because the NLO velocity contains

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

The diagnostic finite-difference backend found

\[
\|v_{\rm with}-v_{\rm without}\| \neq 0.
\]

From `reports/nlo_current/coefficient_derivative_backend_report.md`:

- `dK2_norm = 8.3263427542975321e-01`
- `LC_K3_norm = 6.7375671547562643e-01`
- `LB_K3_norm = 5.8796833753329447e-01`
- `d2K3_norm = 2.5435422318879475e-01`
- `velocity_with_derivatives_norm = 4.5427403565727809e-01`
- `velocity_without_derivatives_norm = 6.5119367148233423e-02`
- `velocity_difference_norm = 4.4359374979437494e-01`

This backend is non-production only. It scales poorly:

- \(O(D)\) first derivatives of full arrays;
- \(O(D^2)\) nested second derivatives for \(K_3\).

## 9. Appendix A Dipole Validation

All five sectors passed:

\[
K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q},\quad K_{JJSJ},\quad K_{JJSSJ}.
\]

| sector | status | max residual | notes |
|---|---|---:|---|
| \(K_{JSJ}\) | passed | `1.6172982426630268e-16` | sign and LO check; symmetric zero-diagonal target kernel |
| \(K_{JSSJ}\) | passed | `1.3877787807814457e-16` | full combined-equation target excluding \(\widetilde K\) contamination |
| \(K_{q\bar q}\) | passed | `5.6325111571547335e-17` | subsection target partial; full target includes subtraction |
| \(K_{JJSJ}\) | passed | `1.1801832636420706e-15` | cubic \((-i)\) convention |
| \(K_{JJSSJ}\) | passed | `1.9087382418229414e-15` | \(\widetilde K\), pure eight-kernel, virtual \(1/3\) exercised |

Dedicated Appendix-target reports also found:

- \(K_{JSSJ}\) subsection-as-full residual was order one in generic tests, so
  the subsection formula is not the full target.
- \(K_{JJSJ}\) raw uncalibrated direct action fails, confirming the cubic
  \((-i)\) convention is meaningful.
- \(K_{JJSSJ}\) raw uncalibrated direct action fails, removing the virtual
  \(1/3\) changes the result, \(\widetilde K\) is nonzero, and both pure
  eight-kernel and tilde-K real terms are exercised.

\[
\boxed{
\text{The full five-kernel NLO Hamiltonian action on the dipole has been validated against Appendix A.}
}
\]

This validates the observable-side sector algebra, signs, normalizations, and
trace ordering in the dense diagnostic implementation. It does not validate
physical coordinate kernels or production evolution.

## 10. Positivity and Pawula caveat

The generalized current representation is algebraic on the region where
\(W>0\). The score and Hessian-score variables require this positivity:

\[
s_A=L_A\log W,
\qquad
H_{AB}=L_As_B.
\]

If an evolution step reaches \(W=0\) or \(W<0\), these score-based variables are
not defined without an additional regularization or a different representation.

The dense NLO normal form contains finite third-order derivative terms through
\(K_3\). Pawula's theorem warns that finite third-order Kramers-Moyal
truncations need not generate a positivity-preserving Markov semigroup.
Therefore ordinary Markov positivity is not guaranteed by the current
representation alone.

The non-production toy diagnostic in
`scripts/nlo_current/check_pawula_positivity_toy.py` checks a scalar periodic
model

\[
\partial_Y W
=-\partial_\theta(K_1W)
+\frac12\partial_\theta^2(K_2W)
-\frac16\partial_\theta^3(K_3W).
\]

Its report,
`reports/nlo_current/pawula_positivity_diagnostic_report.md`, shows that an
LO-like \(K_3=0\), \(K_2>0\) case has no short-step negative mass in the toy
setup, while a pure third-order case produces a positive maximum-principle
warning, negative off-diagonal generator entries, and negative mass for a
constructed near-zero density.

This is only a toy diagnostic. Positivity must be checked separately: first in
such controlled dense diagnostics and later with physical coordinate kernels.
The diagnostic does not prove physical NLO JIMWLK positivity or
non-positivity, and it does not claim production positivity.

## 11. Current Limitations And Next Steps

Remaining blockers before production:

1. Physical coordinate kernels:
   - unbarred versus barred kernels;
   - singularities and regulator/subtraction policy;
   - integration measure;
   - small-lattice diagnostic first.

2. Scalable coefficient derivatives:
   - analytic derivative rules;
   - automatic differentiation;
   - sparse/local structure.

3. Score/Hessian-score strategy:

\[
s_A=L_A\log W,
\qquad
H_{AB}=L_As_B.
\]

   A contracted-Hessian strategy may be preferable:

\[
K_3^{ABC}H_{BC}
\]

   rather than materializing the full Hessian.

4. Non-production NLO flow experiment:
   - tiny lattice only;
   - no physics claims;
   - stability and sanity checks.

5. Paper/demo writing.

This summary does not claim production readiness.

## Appendix: Source Provenance

K_JSJ sign fix and two-generator setup:

- `docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md`
- `reports/nlo_current/kjssj_symmetry_report.md`
- `reports/nlo_current/ordered_lr_current_fd_report.md`
- `reports/nlo_current/left_divergence_identity_report.md`

Kqbarq ordered current:

- `docs/nlo_current/Kqbarq_ordered_current.md`
- `reports/nlo_current/kqbarq_symmetry_report.md`
- `src/nlo_current/two_generator_terms.py`

Two-generator sector summary:

- `docs/nlo_current/two_generator_sector_summary.md`

KJJSJ cubic sector:

- `docs/nlo_current/KJJSJ_cubic_ordered_current.md`
- `reports/nlo_current/kjjsj_cubic_requirements_report.md`
- `reports/nlo_current/kjjsj_appendix_target_validation_report.md`
- `src/nlo_current/three_generator_terms.py`

KJJSSJ cubic sector:

- `docs/nlo_current/KJJSSJ_cubic_ordered_current.md`
- `reports/nlo_current/kjjssj_cubic_requirements_report.md`
- `reports/nlo_current/kjjssj_appendix_target_validation_report.md`
- `src/nlo_current/three_generator_terms.py`

Cubic commutators:

- `docs/nlo_current/cubic_coincident_site_commutators.md`
- `docs/nlo_current/cubic_current_with_commutator_corrections.md`
- `reports/nlo_current/cubic_commutator_corrections_report.md`
- `src/nlo_current/lie_word_algebra.py`
- `src/nlo_current/cubic_commutator_terms.py`

Non-production skeleton and velocity evaluator:

- `docs/nlo_current/nlo_current_skeleton_design.md`
- `docs/nlo_current/nlo_current_map_status.md`
- `reports/nlo_current/nlo_current_skeleton_demo_report.md`
- `src/nlo_current/nlo_current_skeleton.py`
- `src/nlo_current/nlo_velocity_evaluator.py`

Coefficient derivatives:

- `docs/nlo_current/coefficient_derivative_strategy.md`
- `docs/nlo_current/coefficient_derivative_backend_limitations.md`
- `reports/nlo_current/coefficient_derivative_backend_report.md`
- `src/nlo_current/coefficient_derivatives.py`

Appendix A formulas and validation:

- `docs/nlo_current/KLM_appendix_A_dipole_targets_notes.md`
- `docs/nlo_current/full_dipole_validation_plan.md`
- `docs/nlo_current/dipole_validation_status.md`
- `reports/nlo_current/full_dipole_validation_report.md`
- `reports/nlo_current/cubic_i_convention_calibration_report.md`
- `reports/nlo_current/kjssj_appendix_target_validation_report.md`
- `reports/nlo_current/kjjsj_appendix_target_validation_report.md`
- `reports/nlo_current/kjjssj_appendix_target_validation_report.md`
- `src/nlo_current/dipole_appendix_targets.py`
- `src/nlo_current/dipole_hamiltonian_action.py`

Pawula positivity diagnostic:

- `docs/nlo_current/pawula_positivity_diagnostic_plan.md`
- `scripts/nlo_current/check_pawula_positivity_toy.py`
- `tests/nlo_current/test_pawula_positivity_toy.py`
- `reports/nlo_current/pawula_positivity_diagnostic_report.md`

## Physical Density-Side Closure Validation

The dense non-production physical closure diagnostic checks

\[
{\cal G}_{\rm direct}[W]
=
{\cal G}_{\rm current}[W],
\]

where the direct side is the normal-form density operator and the current side
is \(-L_A(v^A W)\). The implemented derivation is recorded in
`docs/nlo_current/physical_density_closure_derivation.md`.

Test densities are positive by construction: a single-link trace density, a
dipole trace density, a multilink nonlinear density, and a constant-density
limit. Scores and ordered Hessian-scores are finite-difference derivatives of
the same \(\log W\); no learned model is used.

The physical projected scan uses an explicit diagnostic finite-grid policy,
active outer derivative index `0`, and the `KJSSJ`/`Kqbarq` physical sectors.
The best projected residual is `3.5453411040275995e-13` absolute. Sparse
synthetic cubic checks verify nonzero Hessian-score closure, raw-cubic
normalization failure, and commutator-correction failure modes.

This validates the dense diagnostic density-side algebra on the tested tiny
lattice. It does not validate production evolution, physical coordinate-kernel
regulator independence, or physical positivity.

Physical density closure sources:

- `docs/nlo_current/physical_density_closure_derivation.md`
- `src/nlo_current/test_densities.py`
- `src/nlo_current/physical_density_operator.py`
- `src/nlo_current/physical_current_divergence.py`
- `src/nlo_current/physical_density_closure.py`
- `scripts/nlo_current/check_physical_density_closure.py`
- `tests/nlo_current/test_physical_density_closure.py`
- `reports/nlo_current/physical_density_closure_report.md`
- `reports/nlo_current/physical_density_closure_failure_modes.md`

## Analytic Coefficient Derivative Status

The finite-difference coefficient derivative backend remains the reference
oracle. A local analytic backend is now implemented and FD-validated for the
two-generator \(dK2^A=L_BK_2^{AB}\) contractions:

- `KJSJ`;
- `KJSSJ`;
- `Kqbarq`, including trace, subtraction, and full contributions.

It is also implemented and FD-validated for the \(K_{JJSJ}\) cubic diagnostic
contractions:

- \(dK2_{\rm comm}^A\) from canonicalized quadratic commutator corrections;
- \((LC\_K3)^{AB}=L_CK_3^{ABC}\);
- \((LB\_K3)^{AC}=L_BK_3^{ABC}\);
- \(d2K3^A=L_BL_CK_3^{ABC}\).

The \(K_{JJSJ}\) validation uses nonzero synthetic dense tensors plus an
expected-zero physical two-site smoke check. Full first-derivative residuals
are at the \(10^{-11}\) level and ordered \(d2K3\) residuals are at the
\(10^{-8}\) level in the stable finite-difference window.

The local primitive rules are calibrated against finite differences in
`tests/nlo_current/test_analytic_lie_derivatives.py`. The physical
two-generator backend is compared against the FD oracle in
`tests/nlo_current/test_analytic_coefficient_derivatives.py` and
`scripts/nlo_current/check_analytic_coefficient_derivatives.py`.

The \(K_{JJSJ}\) cubic backend is compared against the FD oracle in
`tests/nlo_current/test_analytic_cubic_derivatives.py` and
`scripts/nlo_current/check_kjjsj_analytic_cubic_derivatives.py`.

The remaining pending cubic analytic sector is:

- `KJJSSJ`: `LC_K3`, `LB_K3`, and `d2K3`.

`backend="analytic"` supports completed sectors without FD fallback and raises
for pending `KJJSSJ`. The explicit `backend="hybrid_local_fd"` path is
available for mixed diagnostics but is not marked analytic-complete.

Analytic coefficient derivative sources:

- `docs/nlo_current/analytic_lie_derivative_conventions.md`
- `docs/nlo_current/analytic_physical_coefficient_derivative_derivation.md`
- `src/nlo_current/analytic_lie_derivatives.py`
- `src/nlo_current/analytic_two_generator_derivatives.py`
- `src/nlo_current/analytic_cubic_derivatives.py`
- `src/nlo_current/physical_coefficient_derivatives.py`
- `scripts/nlo_current/check_kjjsj_analytic_cubic_derivatives.py`
- `scripts/nlo_current/check_analytic_coefficient_derivatives.py`
- `scripts/nlo_current/benchmark_coefficient_derivative_backends.py`
- `reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`
- `reports/nlo_current/kjjsj_analytic_cubic_benchmark.md`
- `reports/nlo_current/kjjsj_analytic_cubic_failure_modes.md`
- `reports/nlo_current/analytic_coefficient_derivative_validation_report.md`
- `reports/nlo_current/analytic_coefficient_derivative_benchmark.md`
- `reports/nlo_current/analytic_coefficient_derivative_failure_modes.md`
