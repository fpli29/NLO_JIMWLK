# Codex Workflow: Analytic Cubic Coefficient Derivatives for \(K_{JJSJ}\)

## Purpose

Continue from the partial analytic/local physical coefficient-derivative backend.

Current validated state:

1. Analytic Lie-derivative primitives exist and are finite-difference validated:
   \[
   L^aU=i\,t^aU,
   \qquad
   L^aU^\dagger=-i\,U^\dagger t^a,
   \]
   \[
   L^hS_A^{ab}=f^{hac}S_A^{cb},
   \]
   together with ordered second derivatives and trace-word derivatives.

2. Analytic two-generator coefficient derivatives are complete:
   - \(K_{JSJ}\): analytic \(dK2\);
   - \(K_{JSSJ}\): analytic \(dK2\);
   - \(K_{q\bar q}\): analytic trace, subtraction, and full \(dK2\).

3. The cubic physical sectors remain pending:
   - \(K_{JJSJ}\);
   - \(K_{JJSSJ}\).

4. The current `analytic` backend raises for cubic sectors and does not silently fall back.

5. `hybrid_local_fd` is explicitly labeled and remains available as the reference path.

The goal of this workflow is to complete **only the \(K_{JJSJ}\) cubic analytic derivatives**:

\[
(LC\_K3)^{AB}=L_CK_3^{ABC},
\]

\[
(LB\_K3)^{AC}=L_BK_3^{ABC},
\]

\[
d2K3^A=L_BL_CK_3^{ABC}.
\]

Do not implement \(K_{JJSSJ}\) in this task.

---

## Hard Constraints

Do **not**:

- modify production evolution code;
- train score or Hessian-score models;
- optimize for large lattices;
- implement \(K_{JJSSJ}\);
- silently use finite differences inside `backend="analytic"`;
- apply the cubic \((-i)\) normalization more than once;
- assume ordered Lie derivatives commute at coincident sites;
- symmetrize \(LC_K3\), \(LB_K3\), or \(d2K3\);
- discard complex intermediate values;
- loosen tolerances merely to force agreement;
- mark \(K_{JJSJ}\) analytic-complete unless derivative oracle checks and density closure both pass.

Keep the existing finite-difference backend as the reference oracle.

---

## Execution Mode

If `.git` metadata is unavailable, continue in no-git mode.

If Git is available:

```bash
git checkout -b nlo-current-kjjsj-analytic-cubic-derivatives
```

If the worktree is dirty, report dirty files before editing.

Maintain:

```text
reports/nlo_current/file_manifest.md
```

with every created or modified file.

---

## Phase 0: Baseline Audit

Inspect:

```text
src/nlo_current/analytic_lie_derivatives.py
src/nlo_current/analytic_cubic_derivatives.py
src/nlo_current/physical_coefficient_derivatives.py
src/nlo_current/three_generator_terms.py
src/nlo_current/cubic_commutator_terms.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/physical_nlo_current.py
src/nlo_current/physical_cubic_conventions.py
src/nlo_current/nlo_current_skeleton.py
src/nlo_current/nlo_velocity_evaluator.py
src/nlo_current/physical_density_closure.py

tests/nlo_current/test_analytic_lie_derivatives.py
tests/nlo_current/test_analytic_coefficient_derivatives.py
tests/nlo_current/test_physical_density_closure.py

docs/nlo_current/analytic_lie_derivative_conventions.md
docs/nlo_current/analytic_physical_coefficient_derivative_derivation.md

reports/nlo_current/analytic_coefficient_derivative_validation_report.md
reports/nlo_current/analytic_coefficient_derivative_benchmark.md
reports/nlo_current/physical_density_closure_report.md
reports/nlo_current/file_manifest.md
```

Create:

```text
reports/nlo_current/kjjsj_analytic_cubic_start_status.md
```

Record:

- current test count;
- current analytic two-generator status;
- current pending cubic status;
- current finite-difference reference step range;
- current cubic normalization convention;
- current density-closure residuals;
- that this task covers \(K_{JJSJ}\) only.

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

Expected baseline is approximately:

```text
143 passed
```

If all tests pass but the count differs, record the actual result and continue.

---

## Phase 1: \(K_{JJSJ}\) Analytic Derivation Note

Create:

```text
docs/nlo_current/KJJSJ_analytic_cubic_derivatives.md
```

Document the exact coefficient blocks used by the normal-form assembly.

Use the existing implementation and validation notes as source of truth.

The distinct-site blocks are of the form:

\[
A_{LLR}^{dea}(x,y,w)
=
\int_z
K_{JJSJ}(w;x,y;z)\,
f^{bde}S_z^{ba},
\]

\[
B_{LRR}^{ade}(w,x,y)
=
-\int_z
K_{JJSJ}(w;x,y;z)\,
f^{bde}S_z^{ab},
\]

plus the virtual block and coincident-site commutator corrections.

The document must derive:

### First derivative of LLR coefficient

\[
L_C A_{LLR}^{dea}
=
\int_z
K_{JJSJ}\,
f^{bde}\,
L_CS_z^{ba}.
\]

### First derivative of LRR coefficient

\[
L_C B_{LRR}^{ade}
=
-\int_z
K_{JJSJ}\,
f^{bde}\,
L_CS_z^{ab}.
\]

### Ordered second derivative

\[
L_BL_C A_{LLR}^{dea}
=
\int_z
K_{JJSJ}\,
f^{bde}\,
L_BL_CS_z^{ba},
\]

and similarly for \(B_{LRR}\).

State explicitly:

- physical coordinate kernels are independent of Wilson lines;
- Lie derivatives act only on Wilson-line coefficient structures;
- derivative ordering is preserved;
- same-site derivatives obey the tested Lie algebra;
- distinct-site derivatives commute;
- raw physical cubic kernels are complex;
- analytic derivatives are applied to the already KLM-normalized coefficient used by the skeleton;
- the \((-i)\) normalization must not be re-applied inside derivative code.

Include code-path mapping for:

```text
raw physical KJJSJ
-> physical_cubic_conventions
-> normalized real coefficient
-> normal-form K3
-> analytic derivatives
```

---

## Phase 2: Audit the Existing \(K_{JJSJ}\) Assembly

Before writing derivative code, identify every \(K_{JJSJ}\) contribution to:

\[
K_3,\qquad K_{2,\rm comm},\qquad K_{1,\rm comm}.
\]

Create a structured inventory in:

```text
reports/nlo_current/kjjsj_normal_form_inventory.md
```

For each contribution, record:

```text
name:
source function:
tensor order:
index ordering:
Wilson-line dependence:
coordinate-kernel dependence:
raw/normalized dtype:
commutator origin:
expected derivative contribution:
```

Required categories:

1. distinct-site LLR \(K_3\);
2. distinct-site LRR \(K_3\);
3. virtual cubic \(K_3\);
4. quadratic commutator correction \(K_{2,\rm comm}\);
5. linear commutator correction \(K_{1,\rm comm}\), even if zero in the current diagnostic;
6. any metadata-only or structurally zero block.

Do not implement derivatives until this inventory is complete.

---

## Phase 3: Implement Analytic First Derivatives

Update:

```text
src/nlo_current/analytic_cubic_derivatives.py
```

Implement explicit \(K_{JJSJ}\) functions:

```python
def analytic_LC_K3_KJJSJ(
    U,
    physical_terms,
    *,
    sector_data=None,
):
    """
    Return (LC_K3)^{AB} = L_C K3^{ABC}
    for the K_JJSJ sector only.
    """

def analytic_LB_K3_KJJSJ(
    U,
    physical_terms,
    *,
    sector_data=None,
):
    """
    Return (LB_K3)^{AC} = L_B K3^{ABC}
    for the K_JJSJ sector only.
    """
```

Requirements:

- derive from the exact assembled \(K_3\), not from a simplified synthetic coefficient;
- preserve combined-index ordering;
- use local Kronecker support;
- use `left_derivative_adjoint(...)`;
- include LLR, LRR, and virtual structures correctly;
- include any \(K_3\)-level commutator-canonicalization convention already present in the skeleton;
- preserve dtype;
- no global finite-difference calls;
- no silent fallback.

Return optional per-block diagnostics:

```python
{
    "LLR": ...,
    "LRR": ...,
    "virtual": ...,
}
```

so oracle mismatches can be localized.

---

## Phase 4: Implement Analytic Ordered Second Derivative

Update:

```text
src/nlo_current/analytic_cubic_derivatives.py
```

Implement:

```python
def analytic_d2K3_KJJSJ(
    U,
    physical_terms,
    *,
    sector_data=None,
):
    """
    Return d2K3^A = L_B L_C K3^{ABC}
    for the K_JJSJ sector only.
    """
```

Requirements:

1. Use the ordered second derivative primitive:
   ```python
   second_left_derivative_adjoint(...)
   ```

2. Preserve:
   \[
   L_BL_C \neq L_CL_B
   \]
   at coincident sites.

3. Handle:
   - both derivatives acting on the same adjoint Wilson line;
   - distinct-site derivatives;
   - virtual cubic contributions;
   - all index contractions exactly as assembled.

4. Do not implement this by finite-differencing `analytic_LC_K3_KJJSJ`.

5. A temporary diagnostic helper using local FD may be created, but it must be named:

```text
hybrid_local_fd_d2K3_KJJSJ
```

and must not be used by `backend="analytic"`.

Return optional per-block diagnostics:

```python
{
    "LLR": ...,
    "LRR": ...,
    "virtual": ...,
}
```

---

## Phase 5: Analytic Derivatives of Commutator Corrections

Update or create:

```text
src/nlo_current/analytic_cubic_commutator_derivatives.py
```

Implement analytic derivative contributions associated with \(K_{JJSJ}\):

```python
def analytic_dK2_comm_KJJSJ(...)
def analytic_dK1_comm_KJJSJ(...)
```

Requirements:

- differentiate the actual canonicalized commutator-correction tensors;
- do not assume \(K_{1,\rm comm}=0\) as a theorem;
- if the current assembled \(K_{1,\rm comm}\) is structurally zero for \(K_{JJSJ}\), prove this from the current tensor construction and return an explicitly tagged zero;
- compare any nonzero \(dK2_{\rm comm}\) against FD;
- preserve same-site ordering and structure-constant signs.

Metadata should distinguish:

```text
structurally_zero
diagnostically_zero
nonzero
```

---

## Phase 6: Unified Backend Integration

Update:

```text
src/nlo_current/physical_coefficient_derivatives.py
```

Extend `backend="analytic"` so it supports:

- all completed two-generator sectors;
- \(K_{JJSJ}\) cubic sector;
- still raises for \(K_{JJSSJ}\).

Required behavior:

```python
compute_physical_coefficient_derivatives(
    ...,
    backend="analytic",
    sector_mask={"KJJSJ"},
)
```

must return:

- \(LC_K3\);
- \(LB_K3\);
- \(d2K3\);
- \(dK2_{\rm comm}\);
- \(dK1_{\rm comm}\), if relevant;
- by-block metadata;
- no global FD usage.

A request including \(K_{JJSSJ}\) with `backend="analytic"` must still raise a clear pending-sector error.

Do not silently switch to hybrid mode.

---

## Phase 7: FD Oracle Comparison Script

Create:

```text
scripts/nlo_current/check_kjjsj_analytic_cubic_derivatives.py
```

Use:

- deterministic SU(3) configurations;
- the smallest valid physical coordinate set;
- at least two random seeds;
- multiple FD reference steps.

Recommended FD steps:

```text
2e-3
1e-3
5e-4
2.5e-4
```

or the established stable range.

Compare:

\[
LC_K3^{\rm analytic}
\quad\text{vs}\quad
LC_K3^{\rm FD},
\]

\[
LB_K3^{\rm analytic}
\quad\text{vs}\quad
LB_K3^{\rm FD},
\]

\[
d2K3^{\rm analytic}
\quad\text{vs}\quad
d2K3^{\rm FD}.
\]

Also compare commutator-correction derivatives.

For each quantity report:

- maximum absolute residual;
- relative residual;
- real residual;
- imaginary residual;
- norm of analytic result;
- norm of FD result;
- convergence with FD step;
- expected-zero status;
- runtime;
- number of coefficient evaluations.

Do not rely on relative error for structurally zero tensors.

---

## Phase 8: Required Tests

Create or update:

```text
tests/nlo_current/test_analytic_cubic_derivatives.py
```

Required tests:

1. \(K_{JJSJ}\) LLR \(LC_K3\) agrees with FD.
2. \(K_{JJSJ}\) LRR \(LC_K3\) agrees with FD.
3. \(K_{JJSJ}\) virtual \(LC_K3\) agrees with FD.
4. Full \(LC_K3\) agrees with FD.
5. LLR \(LB_K3\) agrees with FD.
6. LRR \(LB_K3\) agrees with FD.
7. Virtual \(LB_K3\) agrees with FD.
8. Full \(LB_K3\) agrees with FD.
9. LLR ordered \(d2K3\) agrees with FD.
10. LRR ordered \(d2K3\) agrees with FD.
11. Virtual \(d2K3\) agrees with FD.
12. Full \(d2K3\) agrees with FD.
13. Same-site derivative order is not silently reversed.
14. Distinct-site derivatives commute where expected.
15. \(K_{2,\rm comm}\) derivative agrees with FD.
16. \(K_{1,\rm comm}\) status is correctly classified.
17. Cubic normalization is applied exactly once.
18. No global FD call occurs under `backend="analytic"`.
19. No silent fallback occurs.
20. No `ComplexWarning`.
21. No imaginary information is silently lost.
22. Sector mask for `KJJSJ` returns zero for all other sectors.
23. Requesting `KJJSSJ` with `backend="analytic"` still raises pending-sector error.

Use monkeypatching or call counters to prove no global FD backend is invoked.

---

## Phase 9: Revalidate Velocity

Update the physical NLO velocity path as needed to support analytic \(K_{JJSJ}\) derivatives.

Compare:

\[
v_{\rm analytic}^{K_{JJSJ}}
\quad\text{vs}\quad
v_{\rm FD}^{K_{JJSJ}}.
\]

Test with:

- nonzero score;
- nonzero ordered Hessian-score;
- zero score;
- zero Hessian-score;
- constant density;
- commutator corrections on/off.

Required diagnostics:

1. full \(K_{JJSJ}\) velocity agreement;
2. Hessian-score contribution agreement;
3. coefficient-derivative contribution agreement;
4. commutator-correction contribution agreement;
5. expected real character after KLM normalization.

---

## Phase 10: Revalidate Density-Side Closure

Use the existing closure framework.

Required checks:

1. Sparse cubic \(K_{JJSJ}\)-only closure using `backend="analytic"`.
2. Positive density with nonzero Hessian-score.
3. Constant-density limit.
4. Omit Hessian-score toggle still fails.
5. Omit coefficient derivatives toggle still fails.
6. Omit commutator corrections toggle still fails where exercised.
7. Remove cubic normalization toggle still fails.
8. Analytic closure agrees with FD-reference closure.

Report:

\[
R_{\rm abs},
\qquad
R_{\rm rel},
\]

and the difference between analytic- and FD-derived velocities.

Do not mark analytic \(K_{JJSJ}\) complete if derivative tensors disagree but closure passes by accidental cancellation.

---

## Phase 11: Benchmark

Update:

```text
scripts/nlo_current/benchmark_coefficient_derivative_backends.py
```

Add a \(K_{JJSJ}\)-only benchmark comparing:

- finite difference;
- analytic;
- hybrid local FD.

Create/update:

```text
reports/nlo_current/kjjsj_analytic_cubic_benchmark.md
```

Record:

- wall time;
- coefficient evaluation count;
- observed speedup;
- residual against FD;
- memory if readily available.

Do not claim asymptotic production scaling from tiny lattices.

---

## Phase 12: Documentation and Reports

Create:

```text
reports/nlo_current/kjjsj_analytic_cubic_validation_report.md
reports/nlo_current/kjjsj_analytic_cubic_failure_modes.md
```

The validation report must include:

- formulas implemented;
- code paths;
- derivative-order convention;
- cubic normalization convention;
- per-block residuals;
- full residuals;
- FD step scan;
- velocity comparison;
- density-closure comparison;
- benchmark;
- pending \(K_{JJSSJ}\) status.

The failure-mode report must include:

- wrong adjoint-index orientation;
- wrong structure-constant sign;
- LLR/LRR index swap;
- reversed derivative order;
- duplicate \((-i)\) normalization;
- omitted virtual block;
- omitted commutator correction;
- accidental complex-to-real cast;
- hidden FD fallback;
- structurally zero relative-error instability.

Update:

```text
docs/nlo_current/analytic_physical_coefficient_derivative_derivation.md
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md
docs/nlo_current/nlo_current_map_status.md
reports/nlo_current/analytic_coefficient_derivative_validation_report.md
reports/nlo_current/analytic_coefficient_derivative_benchmark.md
reports/nlo_current/physical_density_closure_report.md
reports/nlo_current/file_manifest.md
```

Do not state that the full cubic analytic backend is complete. Only \(K_{JJSJ}\) may be marked complete if all acceptance criteria pass.

---

## Phase 13: Acceptance Criteria

This workflow is complete only if:

1. \(K_{JJSJ}\) analytic \(LC_K3\) passes FD checks.
2. \(K_{JJSJ}\) analytic \(LB_K3\) passes FD checks.
3. \(K_{JJSJ}\) analytic ordered \(d2K3\) passes FD checks.
4. LLR, LRR, and virtual blocks are individually validated.
5. Commutator-correction derivative status is validated.
6. Same-site derivative order is preserved.
7. Cubic normalization is applied exactly once.
8. `backend="analytic"` performs no global finite-difference calls.
9. No silent backend fallback occurs.
10. No complex warning or imaginary loss occurs.
11. Analytic and FD \(K_{JJSJ}\) velocities agree.
12. \(K_{JJSJ}\)-only density closure passes with the analytic backend.
13. Hessian-omission and coefficient-derivative-omission toggles remain meaningful.
14. Runtime is benchmarked.
15. \(K_{JJSSJ}\) remains explicitly pending.
16. The full suite passes:

```bash
python3 -m pytest tests/nlo_current -q
```

17. No production evolution is implemented.
18. No score/Hessian-score training is implemented.
19. No physical positivity or regulator-independence claim is made.

If any block fails, stop and classify the mismatch by:

- LLR;
- LRR;
- virtual;
- commutator correction;
- derivative ordering;
- cubic normalization;
- dtype;
- FD convergence window.

Do not force a pass.

---

## Final Codex Response Required

At the end, report:

1. Files created and modified.
2. Test commands and results.
3. Exact \(K_{JJSJ}\) analytic formulas implemented.
4. LLR residuals.
5. LRR residuals.
6. Virtual residuals.
7. Full \(LC_K3\) residual.
8. Full \(LB_K3\) residual.
9. Full \(d2K3\) residual.
10. Commutator-correction derivative status.
11. Velocity agreement.
12. Density-closure agreement.
13. Runtime and observed speedup.
14. Whether any global FD fallback remains.
15. Whether any complex warning or imaginary loss remains.
16. Whether \(K_{JJSJ}\) analytic completion is justified.
17. Confirmation that \(K_{JJSSJ}\) remains pending.
18. Remaining blockers.

Recommended next stage only after this workflow passes:

```text
K_JJSSJ analytic cubic coefficient derivatives
```

followed by:

```text
fully analytic physical coefficient-derivative backend
```

and then:

```text
matrix-free K3:H contraction
```
