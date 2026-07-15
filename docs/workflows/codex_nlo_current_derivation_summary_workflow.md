# Codex Workflow: Consolidate NLO JIMWLK Generalized Current Derivations into a Single Markdown Summary

## Purpose

Before moving into physical KLM kernel integration, consolidate the completed derivations, convention checks, finite-difference tests, commutator algebra, skeleton interface, coefficient-derivative diagnostics, and Appendix A dipole validation into one coherent Markdown document.

This is a documentation / theory-summary workflow only.

Do **not** implement physical kernels.  
Do **not** modify production evolution code.  
Do **not** train score or Hessian-score models.  
Do **not** change existing validated formulas or tests unless a clear typo is found and reported.  
Do **not** invent missing derivations.

The output should be a self-contained but source-grounded summary suitable for later paper notes and for onboarding future Codex runs.

---

## Execution Mode

This workspace may not be a Git repository. If `.git` metadata is unavailable, continue in no-git mode.

If no-git mode is used:

1. Do not modify source code unless explicitly required for documentation links.
2. Prefer adding isolated documentation/report files under:
   - `docs/nlo_current/`
   - `reports/nlo_current/`
3. Maintain/update:
   ```text
   reports/nlo_current/file_manifest.md
   ```
   listing every created or modified file.

If this is a Git repository, create a branch:

```bash
git checkout -b nlo-current-derivation-summary
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect the existing files. Start with:

```text
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
docs/nlo_current/Kqbarq_ordered_current.md
docs/nlo_current/two_generator_sector_summary.md
docs/nlo_current/KJJSJ_cubic_ordered_current.md
docs/nlo_current/KJJSSJ_cubic_ordered_current.md
docs/nlo_current/three_generator_sector_summary.md
docs/nlo_current/cubic_coincident_site_commutators.md
docs/nlo_current/cubic_current_with_commutator_corrections.md
docs/nlo_current/nlo_current_skeleton_design.md
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/coefficient_derivative_strategy.md
docs/nlo_current/coefficient_derivative_backend_limitations.md
docs/nlo_current/full_dipole_validation_plan.md
docs/nlo_current/dipole_validation_status.md
docs/nlo_current/KLM_appendix_A_dipole_targets_notes.md

reports/nlo_current/kjssj_symmetry_report.md
reports/nlo_current/kqbarq_symmetry_report.md
reports/nlo_current/kjjsj_cubic_requirements_report.md
reports/nlo_current/kjjssj_cubic_requirements_report.md
reports/nlo_current/cubic_commutator_corrections_report.md
reports/nlo_current/nlo_current_skeleton_demo_report.md
reports/nlo_current/coefficient_derivative_backend_report.md
reports/nlo_current/full_dipole_validation_report.md
reports/nlo_current/cubic_i_convention_calibration_report.md
reports/nlo_current/kjssj_appendix_target_validation_report.md
reports/nlo_current/kjjsj_appendix_target_validation_report.md
reports/nlo_current/kjjssj_appendix_target_validation_report.md
reports/nlo_current/file_manifest.md
```

Also inspect source files only as needed for signatures and tested conventions:

```text
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/cubic_commutator_terms.py
src/nlo_current/nlo_current_skeleton.py
src/nlo_current/nlo_velocity_evaluator.py
src/nlo_current/coefficient_derivatives.py
src/nlo_current/dipole_hamiltonian_action.py
src/nlo_current/dipole_appendix_targets.py
```

Create:

```text
reports/nlo_current/derivation_summary_start_status.md
```

It should state:

- whether all required docs/reports exist;
- whether the current test suite passes;
- whether all five Appendix A targets are marked passed;
- whether this workflow is documentation-only.

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

before summarizing. Report the result, but do not rewrite code to fix failures in this workflow.

---

## Phase 1: Create the Main Summary Markdown

Create:

```text
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md
```

The document should be self-contained and organized as follows.

---

# NLO JIMWLK Generalized Probability Current: Derivation and Validation Summary

## 1. Scope and conventions

Include:

- Observable-side convention:
  \[
  \frac{d}{dY}{\cal O}=-H{\cal O}.
  \]
- Density-side convention:
  \[
  \partial_Y W=-H^\dagger W.
  \]
- Left/right generator relation:
  \[
  J_R^a(x)=S_A^{ba}(x)J_L^b(x),
  \qquad
  J_L^a(x)=S_A^{ab}(x)J_R^b(x).
  \]
- Dipole:
  \[
  s(u,v)=\frac1{N_c}\mathrm{tr}[U^\dagger(u)U(v)].
  \]
- Clarify the Appendix A target convention:
  - TeX Appendix A formulas are calibrated as \(H_{\rm sector}s\).
  - The current direct-action code uses the calibrated convention documented in the validation reports.
  - Cubic one-\(f\) sectors use the established relation:
    \[
    \text{TeX target}=(-i)\times\text{raw direct action}
    \]
    under the current Hermitian-generator direct-action implementation.

Mention source files used.

---

## 2. LO reference point

Summarize the LO divergence-form current:

\[
\partial_Y W=L_A(\chi^{AB}L_BW).
\]

With score

\[
s_A=L_A\log W,
\]

one has:

\[
L_BW=Ws_B,
\]

and therefore:

\[
\partial_YW=L_A(\chi^{AB}s_BW)=-L_A(v^AW),
\]

with

\[
v^A=-\chi^{AB}s_B.
\]

Also include the generic Itô-style current formula and explain why the derivative-of-\(\chi\) drift cancels in the divergence-form LO JIMWLK current.

---

## 3. General density normal form

State the working normal form:

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W).
\]

Current:

\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

Velocity:

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

where:

\[
H_{BC}=L_Bs_C=L_BL_C\log W.
\]

State the key rule:

\[
\text{k-th order generator terms require derivatives of }\log W\text{ up to order }k-1.
\]

Thus:

- two-generator terms require score;
- three-generator terms require score and Hessian-score.

---

## 4. Two-generator sector

### 4.1 \(K_{JSJ}\)

Include Hamiltonian block:

\[
H_{JSJ}
=
\int K_{JSJ}
\left[
J_LJ_L+J_RJ_R-2J_LSJ_R
\right].
\]

Include the left-basis kernel result:

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
(S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}).
\]

Velocity:

\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
\]

Mention LO sign check:

\[
K_{JSJ}\rightarrow -M/2
\]

gives positive LO diffusion kernel.

### 4.2 Ordered \(J_LAJ_R\) lemma

State the block:

\[
H_{LR}[A]=\int A^{ab}(x,y;U)J_L^a(x)J_R^b(y).
\]

Using:

\[
J_R^b(y)=S_y^{hb}L_y^h,
\]

and tested identity:

\[
L_y^hS_y^{hb}=0,
\]

the density-side current is:

\[
J_{LR}^{(y,h)}
=
S_y^{hb}L_x^a[A^{ab}W],
\]

and velocity:

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

Mention finite-difference residuals from the relevant report.

### 4.3 \(K_{JSSJ}\)

Include coefficient:

\[
A_{JSSJ}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}).
\]

Current:

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

State diagnostic result:

\[
\frac{\|C-C^T\|}{\|C\|}\sim O(1),
\]

so \(K_{JSSJ}\) must not be simplified to a symmetric pure \(v=-\chi s\) score current.

### 4.4 \(K_{q\bar q}\)

Include coefficient:

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

Current:

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

State that \(z'=z\) subtraction identity passed:

\[
2\mathrm{tr}(U^\dagger t^aUt^b)-S_A^{ab}=0.
\]

State diagnostic asymmetry ratio is order-one, so this also remains a generic ordered-current / commutator-drift term.

---

## 5. Three-generator sector: distinct-site cubic current

### 5.1 LLR block

State:

\[
H_{LLR}[A]=A^{dea}J_L^dJ_L^eJ_R^a.
\]

The tested current:

\[
v_{LLR}^{(w,h)}
=
-S_w^{ha}
\frac1W L_y^eL_x^d[A^{dea}W].
\]

Expand:

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

This is the first appearance of Hessian-score.

### 5.2 LRR and virtual blocks

Summarize the tested LRR current lemma and the canonical virtual LLL/RRR block.

State that the \(1/3\) virtual factor was tested and changing it changes the result by factor three or fails the validation.

### 5.3 \(K_{JJSJ}\)

Include coefficient blocks:

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

State validation:

- LLR sign passed.
- LRR sign passed.
- virtual \(1/3\) passed.
- Hessian-score term nonzero.
- Appendix A target passed using cubic \((-i)\) convention.

### 5.4 \(K_{JJSSJ}\)

Include coefficient blocks:

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

State KLM-like simultaneous antisymmetry:

\[
K(w;x,y;z,z')=-K(w;y,x;z',z).
\]

State validation:

- LLR sign passed.
- LRR sign passed.
- virtual \(1/3\) passed.
- Hessian-score nonzero.
- Appendix A passed including \(\widetilde K\), pure eight-kernel contribution, and virtual term.

---

## 6. Coincident-site commutators

Explain why distinct-site validation is not enough:

\[
[L_x^a,L_y^b]=0\quad x\neq y,
\]

but:

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

State canonical word ordering for combined index \(A=(x,a)\).

State the swap rule:

\[
L_x^aL_x^b
=
L_x^bL_x^a
+
f^{abc}L_x^c.
\]

State the canonicalized result structure:

\[
L_AL_BL_CF
=
\sum c\,L_{A'}L_{B'}L_{C'}F
+
\sum d\,L_{A'}L_{B'}F
+
\sum e\,L_{A'}F.
\]

Thus cubic sectors induce:

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

Include the diagnostic results:

- \(K_{JJSJ}\): nonzero \(K_{2,\rm comm}\), zero \(K_{1,\rm comm}\) in the specific diagnostic.
- \(K_{JJSSJ}\): nonzero \(K_{2,\rm comm}\), nonzero \(K_{1,\rm comm}\).

Do not state diagnostic zero as theorem.

---

## 7. NLO current skeleton

Describe:

```python
@dataclass
class NLOCurrentTerms:
    K1: np.ndarray
    K2: np.ndarray
    K3: np.ndarray
    metadata: dict
```

with shapes:

\[
D=8N_{\rm site},
\]

\[
K_1:(D,),\quad K_2:(D,D),\quad K_3:(D,D,D).
\]

State sector map:

| sector | normal-form contribution | density derivatives needed |
|---|---|---|
| \(K_{JSJ}\) | \(K_2\) | score |
| \(K_{JSSJ}\) | \(K_2\) | score + coefficient drift |
| \(K_{q\bar q}\) | \(K_2\) | score + coefficient drift |
| \(K_{JJSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |
| \(K_{JJSSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |

Mention:

- non-production skeleton assembled all five sectors;
- commutator corrections are included and toggleable;
- diagnostic velocity evaluation works with supplied score and Hessian-score.

---

## 8. Coefficient derivatives

State why coefficient derivatives are needed:

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

Include diagnostic finding:

\[
\|v_{\rm with}-v_{\rm without}\| \neq 0,
\]

and report the actual norms from:

```text
reports/nlo_current/coefficient_derivative_backend_report.md
```

Mention finite-difference backend is non-production only and scales poorly:

- \(O(D)\) first derivatives of full arrays;
- \(O(D^2)\) nested second derivatives for \(K_3\).

---

## 9. Appendix A dipole validation

Summarize final status:

All five sectors passed:

\[
K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q},\quad K_{JJSJ},\quad K_{JJSSJ}.
\]

Include a table:

| sector | status | max residual | notes |
|---|---|---:|---|
| \(K_{JSJ}\) | passed | ... | sign/LO check |
| \(K_{JSSJ}\) | passed | ... | full combined-equation target excluding \(\widetilde K\) contamination |
| \(K_{q\bar q}\) | passed | ... | subsection target partial; full target includes subtraction |
| \(K_{JJSJ}\) | passed | ... | cubic \((-i)\) convention |
| \(K_{JJSSJ}\) | passed | ... | \(\widetilde K\), pure eight-kernel, virtual \(1/3\) exercised |

Pull residuals from:

```text
reports/nlo_current/full_dipole_validation_report.md
reports/nlo_current/kjssj_appendix_target_validation_report.md
reports/nlo_current/kjjsj_appendix_target_validation_report.md
reports/nlo_current/kjjssj_appendix_target_validation_report.md
```

State clearly:

\[
\boxed{
\text{The full five-kernel NLO Hamiltonian action on the dipole has been validated against Appendix A.}
}
\]

---

## 10. Current limitations and next steps

List remaining blockers before production:

1. Physical coordinate kernels:
   - unbarred vs barred;
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
   Mention possible contracted-Hessian strategy:
   \[
   K_3^{ABC}H_{BC}
   \]
   rather than materializing full \(H\).

4. Non-production NLO flow experiment:
   - tiny lattice only;
   - no physics claims;
   - stability/sanity checks.

5. Paper/demo writing.

---

## Phase 2: Add Source Index / Provenance Appendix

At the end of the same document, add:

```markdown
## Appendix: Source provenance
```

Include a list mapping each major claim to source docs/reports.

Example:

```text
K_JSJ sign fix:
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md

Kqbarq target pass:
reports/nlo_current/full_dipole_validation_report.md
src/nlo_current/dipole_appendix_targets.py

Cubic commutators:
docs/nlo_current/cubic_coincident_site_commutators.md
reports/nlo_current/cubic_commutator_corrections_report.md
```

Do not fabricate line numbers if not available. File paths and section names are enough.

---

## Phase 3: Create a Short Executive Summary

Create:

```text
docs/nlo_current/NLO_JIMWLK_generalized_current_executive_summary.md
```

This should be 1–2 pages max and include:

1. Main result:
   \[
   \partial_YW=-L_AJ^A_{\rm NLO}
   \]
   with \(J^A\) involving \(K_1,K_2,K_3\).

2. Main difference from LO:
   - LO needs score only.
   - NLO cubic terms require Hessian-score.

3. Validation milestone:
   - all five sectors passed Appendix A dipole validation.

4. Remaining work:
   - physical kernels;
   - scalable coefficient derivatives;
   - score/Hessian-score estimator;
   - non-production flow experiment.

---

## Phase 4: Consistency Checks

Create:

```text
reports/nlo_current/derivation_summary_consistency_report.md
```

This report should check:

1. The summary document exists.
2. The executive summary exists.
3. The full test suite passes.
4. All five sectors are marked passed in `dipole_validation_status.md` or the latest validation report.
5. The summary does not claim production readiness.
6. The summary explicitly states:
   - synthetic/dense/non-production status;
   - coefficient derivatives are diagnostic finite-difference only;
   - physical kernels are next step;
   - Hessian-score estimator is not implemented.

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 5: Acceptance Criteria

Stop when all are true:

1. Main summary exists:
   ```text
   docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md
   ```

2. Executive summary exists:
   ```text
   docs/nlo_current/NLO_JIMWLK_generalized_current_executive_summary.md
   ```

3. Consistency report exists:
   ```text
   reports/nlo_current/derivation_summary_consistency_report.md
   ```

4. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

5. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

6. No production evolution code modified.
7. No physical kernels implemented in this workflow.
8. No score/Hessian-score model training implemented.
9. Summary does not claim production readiness.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether the main derivation summary was created.
4. Whether the executive summary was created.
5. Whether all five Appendix A sectors are listed as passed.
6. Whether production-readiness caveats are included.
7. Any source file/report that was missing or inconsistent.
8. Recommended next step:
   - physical kernel integration workflow;
   - then scalable coefficient derivatives / score-Hessian strategy.
