# Codex Workflow: Analytic Cubic Coefficient Derivatives for \(K_{JJSSJ}\)

## Purpose

Continue from the validated partial analytic/local coefficient-derivative backend.

Current validated state:

1. Analytic two-generator derivatives are complete:
   \[
   K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q}.
   \]

2. Analytic \(K_{JJSJ}\) cubic derivatives are complete and FD-validated:
   \[
   LC\_K3,\qquad LB\_K3,\qquad d2K3,
   \]
   including the relevant commutator-correction derivative path.

3. The remaining analytic gap is the \(K_{JJSSJ}\) sector:
   \[
   (LC\_K3)^{AB}=L_CK_3^{ABC},
   \]
   \[
   (LB\_K3)^{AC}=L_BK_3^{ABC},
   \]
   \[
   d2K3^A=L_BL_CK_3^{ABC}.
   \]

4. `backend="analytic"` must still raise whenever \(K_{JJSSJ}\) is requested until this workflow passes.

The goal of this workflow is to complete and validate the analytic/local \(K_{JJSSJ}\) coefficient derivatives, then revalidate velocity and density-side closure.

This task must finish the fully analytic physical coefficient-derivative backend for all five unbarred physical sectors.

---

## Hard Constraints

Do **not**:

- change the physical \(K_{JJSSJ}\) kernel formula;
- change the established raw-complex to KLM-normalized-real convention;
- silently use global finite differences under `backend="analytic"`;
- symmetrize ordered Lie derivatives;
- assume \(z\neq z'\);
- merge coincident-site and distinct-site product rules without explicit tests;
- omit virtual \(1/3\) structures;
- omit pure-eight-kernel or \(\widetilde K\)-related contributions that are present in the current normal-form assembly;
- discard complex intermediates;
- weaken tests merely to obtain a pass;
- train score/Hessian models;
- add production evolution;
- make positivity or regulator-independence claims.

Keep the existing finite-difference implementation as the reference oracle.

If an assembly convention is ambiguous, derive it from the currently passing physical dipole and density-closure tests. Do not guess.

---

## Execution Mode

If `.git` metadata is available, create a branch:

```bash
git checkout -b nlo-current-kjjssj-analytic-cubic-derivatives
```

If Git is unavailable or push access is broken, continue in no-git mode.

Before editing, run:

```bash
git status --short 2>/dev/null || true
```

Record the starting filesystem state.

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
src/nlo_current/analytic_cubic_commutator_derivatives.py
src/nlo_current/physical_coefficient_derivatives.py
src/nlo_current/three_generator_terms.py
src/nlo_current/cubic_commutator_terms.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/physical_cubic_conventions.py
src/nlo_current/physical_kernel_adapter.py
src/nlo_current/physical_kernels.py
src/nlo_current/physical_nlo_current.py
src/nlo_current/physical_density_closure.py

tests/nlo_current/test_analytic_cubic_derivatives.py
tests/nlo_current/test_kjjssj_coefficients.py
tests/nlo_current/test_kjjssj_cubic_current.py
tests/nlo_current/test_cubic_commutator_end_to_end.py
tests/nlo_current/test_physical_density_closure.py

docs/nlo_current/KJJSJ_analytic_cubic_derivatives.md
docs/nlo_current/analytic_physical_coefficient_derivative_derivation.md
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md

reports/nlo_current/kjjsj_analytic_cubic_validation_report.md
reports/nlo_current/analytic_coefficient_derivative_validation_report.md
reports/nlo_current/physical_density_closure_report.md
reports/nlo_current/file_manifest.md
```

Create:

```text
reports/nlo_current/kjjssj_analytic_cubic_start_status.md
```

Record:

- current test count;
- current \(K_{JJSJ}\) residuals;
- current closure residuals;
- current backend routing;
- current cubic normalization convention;
- current finite-grid policy;
- current known \(K_{JJSSJ}\) tensor blocks;
- that this task covers \(K_{JJSSJ}\) only.

Run the baseline suite:

```bash
python3 -m pytest tests/nlo_current -q
```

Expected count is approximately:

```text
151 passed
```

If all tests pass but the count differs, record the actual result and continue.

Because the current full suite is slow, also record the runtime.

---

## Phase 1: Exact \(K_{JJSSJ}\) Normal-Form Inventory

Before implementing derivatives, create:

```text
reports/nlo_current/kjjssj_normal_form_inventory.md
```

Identify every \(K_{JJSSJ}\)-derived contribution to:

\[
K_3,\qquad K_{2,\mathrm{comm}},\qquad K_{1,\mathrm{comm}}.
\]

The inventory must include all currently assembled categories:

1. distinct-site LLR cubic block;
2. distinct-site LRR cubic block;
3. virtual LLL block;
4. virtual RRR block;
5. virtual \(1/3\) factor;
6. pure-eight-kernel contribution;
7. \(\widetilde K\)-related contribution, if routed through the \(K_{JJSSJ}\) assembly path;
8. quadratic commutator correction;
9. linear commutator correction;
10. structurally zero or metadata-only blocks.

For each block record:

```text
block_name:
source_file:
source_function:
input kernel:
raw or normalized:
tensor order:
combined-index ordering:
site ordering:
color-index ordering:
Wilson-line factors:
coordinate weights:
coincident-site policy:
expected first derivative:
expected second derivative:
commutator origin:
dtype:
```

This inventory is a hard prerequisite.

Do not write derivative code before the inventory is complete.

---

## Phase 2: Derivation Document

Create:

```text
docs/nlo_current/KJJSSJ_analytic_cubic_derivatives.md
```

Document the exact coefficient structures used by the code.

The core distinct-site blocks are of the form:

\[
A_{LLR}^{dea}
=
\int_{z,z'}
K_{JJSSJ}\,
f^{acb}
S_z^{dc}
S_{z'}^{eb},
\]

\[
B_{LRR}^{ade}
=
-\int_{z,z'}
K_{JJSSJ}\,
f^{acb}
S_z^{cd}
S_{z'}^{be}.
\]

Use the exact index ordering from the code even if it differs from this schematic display.

Derive the first derivative product rule:

\[
L_C(S_zS_{z'})
=
(L_CS_z)S_{z'}
+
S_z(L_CS_{z'}).
\]

Derive the ordered second derivative:

\[
\begin{aligned}
L_BL_C(S_zS_{z'})
={}&
(L_BL_CS_z)S_{z'}
+(L_CS_z)(L_BS_{z'})
\\
&+
(L_BS_z)(L_CS_{z'})
+S_z(L_BL_CS_{z'}).
\end{aligned}
\]

State explicitly that the two cross terms are generally distinct because color and derivative ordering are retained.

Document separately:

### Distinct-site case

For \(z\neq z'\), derivatives on different sites commute.

### Coincident-site case

For \(z=z'\), ordered Lie derivatives must preserve:

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

Do not replace ordered second derivatives by a symmetrized Hessian.

Document:

- virtual terms;
- pure-eight-kernel terms;
- commutator corrections;
- cubic normalization;
- exact code-path mapping.

Required code-path map:

```text
raw physical KJJSSJ
-> physical_cubic_conventions
-> normalized real cubic coefficient
-> normal-form K3/K2_comm/K1_comm
-> analytic local derivatives
-> physical velocity
-> density closure
```

---

## Phase 3: Product-Derivative Primitives

Update:

```text
src/nlo_current/analytic_lie_derivatives.py
```

Add explicit product primitives if they do not already exist:

```python
def left_derivative_adjoint_product(
    SA_left,
    left_site,
    SA_right,
    right_site,
    derivative_site,
    derivative_color,
):
    """
    Return L_C [S_A(left_site) S_A(right_site)]
    with exact index ordering preserved.
    """

def second_left_derivative_adjoint_product(
    SA_left,
    left_site,
    SA_right,
    right_site,
    first_site,
    first_color,
    second_site,
    second_color,
):
    """
    Return ordered L_B L_C [S_A(left_site) S_A(right_site)].
    """
```

Requirements:

- explicit Kronecker locality;
- preserve matrix/index orientation;
- handle \(z=z'\) without special-casing away noncommutativity;
- use existing validated single-adjoint primitives;
- return complex-capable dtype;
- no finite differences;
- include per-term diagnostics:
  ```text
  LL_left
  cross_Cleft_Bright
  cross_Bleft_Cright
  LL_right
  ```

Add tests comparing these primitives directly against FD on:

1. \(z\neq z'\);
2. \(z=z'\);
3. same color;
4. different colors;
5. reversed derivative order.

---

## Phase 4: Analytic First Derivatives for \(K_{JJSSJ}\)

Update:

```text
src/nlo_current/analytic_cubic_derivatives.py
```

Implement:

```python
def analytic_LC_K3_KJJSSJ(
    U,
    physical_terms,
    *,
    sector_data=None,
    return_by_block=False,
):
    """
    Return (LC_K3)^{AB}=L_C K3^{ABC}
    for K_JJSSJ only.
    """

def analytic_LB_K3_KJJSSJ(
    U,
    physical_terms,
    *,
    sector_data=None,
    return_by_block=False,
):
    """
    Return (LB_K3)^{AC}=L_B K3^{ABC}
    for K_JJSSJ only.
    """
```

Requirements:

- differentiate the exact assembled \(K_3\);
- include all LLR/LRR/virtual blocks;
- include all pure-eight-kernel contributions;
- include any \(\widetilde K\)-routed pieces present in the normal form;
- preserve coordinate quadrature weights;
- preserve combined-index ordering;
- no FD;
- no silent fallback;
- cubic normalization applied exactly once upstream.

Return optional diagnostics:

```python
{
    "LLR": ...,
    "LRR": ...,
    "virtual_LLL": ...,
    "virtual_RRR": ...,
    "pure_eight": ...,
    "tilde_related": ...,
}
```

Use explicit expected-zero metadata for structurally zero blocks.

---

## Phase 5: Analytic Ordered Second Derivative

Update:

```text
src/nlo_current/analytic_cubic_derivatives.py
```

Implement:

```python
def analytic_d2K3_KJJSSJ(
    U,
    physical_terms,
    *,
    sector_data=None,
    return_by_block=False,
):
    """
    Return d2K3^A=L_B L_C K3^{ABC}
    for K_JJSSJ only.
    """
```

Requirements:

1. Use explicit ordered second-derivative product rules.
2. Include:
   - both derivatives on \(S_z\);
   - both derivatives on \(S_{z'}\);
   - \(C\) on \(S_z\), \(B\) on \(S_{z'}\);
   - \(B\) on \(S_z\), \(C\) on \(S_{z'}\).
3. Preserve \(z=z'\) ordering.
4. Preserve virtual and pure-eight contributions.
5. Do not finite-difference the analytic first derivative.
6. Do not symmetrize the cross terms.
7. Do not collapse ordered combined indices.
8. Return per-block and per-product-rule diagnostics.

A local-FD helper may exist only as:

```text
hybrid_local_fd_d2K3_KJJSSJ
```

and must never be called by `backend="analytic"`.

---

## Phase 6: Commutator-Correction Derivatives

Update:

```text
src/nlo_current/analytic_cubic_commutator_derivatives.py
```

Implement:

```python
def analytic_dK2_comm_KJJSSJ(...)
def analytic_dK1_comm_KJJSSJ(...)
```

Requirements:

- derive from the actual canonicalized commutator tensors;
- retain the virtual \(1/3\) structure;
- retain pure-eight-kernel contributions;
- classify each correction as:
  ```text
  structurally_zero
  diagnostically_zero
  nonzero
  ```
- compare nonzero pieces independently against FD;
- do not infer zero from one random configuration;
- preserve structure-constant signs and derivative order.

Create block-level outputs for:

```text
LL coincidence
RR coincidence
mixed coincidence
virtual correction
pure-eight correction
```

where applicable.

---

## Phase 7: Unified Backend Completion

Update:

```text
src/nlo_current/physical_coefficient_derivatives.py
```

After all \(K_{JJSSJ}\) checks pass, extend:

```python
backend="analytic"
```

to support all physical sectors:

```text
KJSJ
KJSSJ
Kqbarq
KJJSJ
KJJSSJ
```

Required behavior:

```python
compute_physical_coefficient_derivatives(
    ...,
    backend="analytic",
    sector_mask={"KJJSSJ"},
)
```

must return:

- \(LC_K3\);
- \(LB_K3\);
- \(d2K3\);
- commutator-correction derivative data;
- per-block metadata;
- zero global FD calls.

Once complete, the metadata should state:

```text
analytic_backend_complete_for_unbarred_physical_sectors: true
```

Do not set that flag before all acceptance criteria pass.

Keep:

```text
backend="finite_difference"
backend="hybrid_local_fd"
```

for reference and diagnostics.

No silent fallback is permitted.

---

## Phase 8: Independent Oracle Validation

Create:

```text
scripts/nlo_current/check_kjjssj_analytic_cubic_derivatives.py
```

Use:

- deterministic SU(3) configurations;
- smallest valid physical coordinate set;
- at least two random seeds;
- multiple FD steps;
- both distinct-site and coincident-site-sensitive setups.

Recommended FD scan:

```text
2e-3
1e-3
5e-4
2.5e-4
```

or the established stable window.

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

Also compare:

- \(dK2_{\rm comm}\);
- \(dK1_{\rm comm}\);
- each block;
- each product-rule contribution.

Report:

- max absolute residual;
- relative residual;
- norm of analytic result;
- norm of FD result;
- real residual;
- imaginary residual;
- FD convergence trend;
- runtime;
- coefficient evaluation count;
- expected-zero status.

Do not judge structurally zero tensors using relative error alone.

---

## Phase 9: Required Tests

Create:

```text
tests/nlo_current/test_kjjssj_analytic_cubic_derivatives.py
```

Required tests:

### Product primitives

1. First derivative of \(S_zS_{z'}\), \(z\neq z'\), agrees with FD.
2. First derivative of \(S_zS_z\) agrees with FD.
3. Ordered second derivative for \(z\neq z'\) agrees with FD.
4. Ordered second derivative for \(z=z'\) agrees with FD.
5. Reversing derivative order changes same-site results consistently with the Lie algebra.
6. Distinct-site derivative order commutes.

### Cubic blocks

7. LLR \(LC_K3\) agrees with FD.
8. LRR \(LC_K3\) agrees with FD.
9. Virtual \(LC_K3\) agrees with FD or is correctly tagged zero.
10. Pure-eight \(LC_K3\) agrees with FD.
11. Full \(LC_K3\) agrees with FD.
12. LLR \(LB_K3\) agrees with FD.
13. LRR \(LB_K3\) agrees with FD.
14. Virtual \(LB_K3\) status is correct.
15. Pure-eight \(LB_K3\) agrees with FD.
16. Full \(LB_K3\) agrees with FD.
17. LLR \(d2K3\) agrees with FD.
18. LRR \(d2K3\) agrees with FD.
19. Virtual \(d2K3\) status is correct.
20. Pure-eight \(d2K3\) agrees with FD.
21. Full \(d2K3\) agrees with FD.

### Commutator corrections

22. \(dK2_{\rm comm}\) agrees with FD.
23. \(dK1_{\rm comm}\) status is correctly classified.
24. Omitting commutator corrections produces a measurable mismatch in a configuration that exercises them.
25. Virtual \(1/3\) factor is tested.
26. Pure-eight correction is tested.

### Backend integrity

27. `backend="analytic"` performs no global FD calls.
28. No silent fallback occurs.
29. Cubic normalization is applied exactly once.
30. No `ComplexWarning`.
31. No imaginary information is silently lost.
32. Sector mask returns zero for unrelated sectors.
33. Full analytic backend supports all five physical sectors.

Use monkeypatching or explicit call counters to prove no FD path is invoked.

---

## Phase 10: Velocity Validation

Compare:

\[
v_{\rm analytic}^{K_{JJSSJ}}
\quad\text{vs}\quad
v_{\rm FD}^{K_{JJSSJ}}.
\]

Use test densities with:

- nonzero score;
- nonzero ordered Hessian-score;
- zero score;
- zero Hessian-score;
- constant density.

Report separately:

- \(d2K3\) drift contribution;
- \((L_CK_3)s_B\);
- \((L_BK_3)s_C\);
- \(K_3H_{BC}\);
- \(K_3s_Bs_C\);
- commutator-correction contributions.

Required toggles:

1. omit Hessian-score;
2. omit coefficient derivatives;
3. omit commutator corrections;
4. remove cubic normalization;
5. remove pure-eight contribution;
6. alter virtual \(1/3\) factor.

Each toggle must either worsen the chosen diagnostic or be explicitly documented as not exercised by that configuration. Do not claim a toggle test passed if the selected configuration gives an identically zero difference.

---

## Phase 11: Density-Side Closure

Use the existing independent density-side closure stack.

Required checks:

1. \(K_{JJSSJ}\)-only analytic closure;
2. sparse coincident-site-sensitive closure;
3. positive density with nonzero ordered Hessian-score;
4. constant-density limit;
5. analytic velocity versus FD velocity;
6. omit-Hessian toggle;
7. omit-coefficient-derivative toggle;
8. omit-commutator-correction toggle;
9. remove cubic-normalization toggle;
10. pure-eight omission toggle;
11. virtual \(1/3\) omission toggle.

Report:

\[
R_{\rm abs},
\qquad
R_{\rm rel}.
\]

After the sector-only checks pass, run one all-sector analytic closure:

\[
K_{JSJ}+K_{JSSJ}+K_{q\bar q}+K_{JJSJ}+K_{JJSSJ}.
\]

The direct density operator must remain independent of the analytic coefficient-derivative implementation. Do not reuse analytic derivative tensors on both sides of the closure identity.

---

## Phase 12: Test Runtime Control

The current suite is already slow. Do not place full FD scans inside default pytest.

Required structure:

### Default pytest

Use:

- one deterministic seed;
- one stable FD reference step;
- smallest valid lattice;
- block-level analytic tests;
- no benchmark loops.

### Validation script

Place:

- multiple seeds;
- multiple FD steps;
- full convergence scans;
- all-sector expensive closure.

### Optional slow marker

Add:

```python
@pytest.mark.slow
```

for tests that exceed a reasonable unit-test runtime.

Do not silently skip important tests. Document the command:

```bash
python3 -m pytest tests/nlo_current -q
```

and, if introduced:

```bash
python3 -m pytest tests/nlo_current -q -m slow
```

Update pytest configuration if needed so unknown-marker warnings do not occur.

Target: keep the default full suite materially below the projected runtime of embedding every KJJSSJ FD scan in pytest.

---

## Phase 13: Benchmark

Update:

```text
scripts/nlo_current/benchmark_coefficient_derivative_backends.py
```

Add:

- \(K_{JJSSJ}\)-only benchmark;
- full all-sector benchmark.

Compare:

```text
analytic
finite_difference
hybrid_local_fd
```

Record:

- wall time;
- coefficient evaluation count;
- observed speedup;
- peak memory if readily available;
- output residual against FD.

Create:

```text
reports/nlo_current/kjjssj_analytic_cubic_benchmark.md
```

Do not claim production asymptotic scaling from tiny-lattice timings.

---

## Phase 14: Documentation and Reports

Create:

```text
reports/nlo_current/kjjssj_analytic_cubic_validation_report.md
reports/nlo_current/kjjssj_analytic_cubic_failure_modes.md
```

The validation report must include:

- exact formulas;
- block inventory;
- product-rule decomposition;
- coincident-site tests;
- per-block residuals;
- full residuals;
- FD scan;
- velocity comparison;
- density closure;
- benchmark;
- backend routing;
- final analytic-completeness status.

The failure-mode report must include:

- wrong adjoint-index orientation;
- LLR/LRR transposition;
- omitted cross derivative;
- duplicate cross derivative;
- reversed derivative order;
- accidental symmetrization;
- wrong same-site commutator;
- missing virtual \(1/3\);
- missing pure-eight term;
- missing \(\widetilde K\)-routed term;
- duplicate or missing \((-i)\) normalization;
- hidden FD fallback;
- finite-grid policy mismatch;
- zero-tensor relative-error instability;
- complex-to-real cast.

Update:

```text
docs/nlo_current/analytic_physical_coefficient_derivative_derivation.md
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/physical_kernel_status.md

reports/nlo_current/analytic_coefficient_derivative_validation_report.md
reports/nlo_current/analytic_coefficient_derivative_benchmark.md
reports/nlo_current/physical_nlo_current_assembly_report.md
reports/nlo_current/physical_density_closure_report.md
reports/nlo_current/file_manifest.md
```

Only after all acceptance criteria pass may the documentation state:

```text
The analytic/local physical coefficient-derivative backend is complete
for all five implemented unbarred physical sectors on controlled
tiny-lattice diagnostics.
```

Do not call it production-ready.

---

## Phase 15: Acceptance Criteria

This workflow is complete only if:

1. Product derivatives for \(S_zS_{z'}\) pass FD checks.
2. Coincident-site ordered second derivatives pass.
3. \(K_{JJSSJ}\) analytic \(LC_K3\) passes FD checks.
4. \(K_{JJSSJ}\) analytic \(LB_K3\) passes FD checks.
5. \(K_{JJSSJ}\) analytic \(d2K3\) passes FD checks.
6. LLR and LRR blocks pass independently.
7. Virtual blocks are validated.
8. Virtual \(1/3\) is validated.
9. Pure-eight contributions are validated.
10. \(\widetilde K\)-routed contributions are validated if present.
11. Commutator-correction derivatives are validated.
12. No global FD call occurs under `backend="analytic"`.
13. No silent fallback occurs.
14. Cubic normalization is applied exactly once.
15. No complex warning or imaginary loss occurs.
16. Analytic and FD velocities agree.
17. \(K_{JJSSJ}\)-only density closure passes.
18. All-sector analytic density closure passes.
19. Toggle diagnostics are meaningfully exercised.
20. Runtime is benchmarked.
21. The full default suite passes:
    ```bash
    python3 -m pytest tests/nlo_current -q
    ```
22. Slow validation passes if a slow marker is introduced.
23. No production evolution is added.
24. No score/Hessian-score training is added.
25. No positivity or regulator-independence claim is added.

If any requirement fails, leave \(K_{JJSSJ}\) pending and report the exact failing block.

Do not mark the full analytic backend complete by partial cancellation in the velocity or closure.

---

## Final Codex Response Required

At the end, report:

1. Files created and modified.
2. Test commands and exact results.
3. Exact product-derivative primitives implemented.
4. LLR residuals.
5. LRR residuals.
6. Virtual residuals/status.
7. Pure-eight residuals/status.
8. \(\widetilde K\)-related residuals/status.
9. Full \(LC_K3\) residual.
10. Full \(LB_K3\) residual.
11. Full \(d2K3\) residual.
12. Commutator-correction derivative results.
13. Same-site ordered derivative validation.
14. Velocity analytic-vs-FD agreement.
15. \(K_{JJSSJ}\)-only closure.
16. All-sector analytic closure.
17. Toggle results.
18. Runtime and observed speedup.
19. Whether any global FD fallback remains.
20. Whether any complex warning or imaginary loss remains.
21. Whether the analytic backend is now complete for all five physical sectors.
22. Remaining blockers.
23. Confirmation that no production, training, positivity, or regulator claim was added.

Recommended next stage only after this workflow passes:

```text
matrix-free contracted Hessian-score evaluation
```

specifically:

\[
K_3^{ABC}H_{BC}
\]

without constructing the full dense ordered Hessian.
