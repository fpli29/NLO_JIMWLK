# Codex Workflow: Analytic/Local Physical Coefficient Derivatives for NLO JIMWLK

## Purpose

Continue from the completed physical density-side closure validation.

Current validated state:

1. All five unbarred physical NLO kernels are implemented:
   \[
   K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q},\quad K_{JJSJ},\quad K_{JJSSJ}.
   \]

2. All five physical dipole-sector checks pass.

3. The physical generalized-current assembly is implemented:
   \[
   \text{physical kernels}
   \rightarrow
   K_1,K_2,K_3
   \rightarrow
   \text{coefficient derivatives}
   \rightarrow
   v_{\rm NLO}.
   \]

4. The density-side closure has been validated on controlled tiny lattices:
   \[
   {\cal G}_{\rm direct}[W]
   =
   -L_A(v^A W).
   \]

5. The current coefficient-derivative backend is finite-difference and non-production.

The goal of this workflow is to replace the expensive global finite-difference coefficient derivatives with a source-grounded analytic/local backend while preserving the existing finite-difference implementation as the reference oracle.

The required contractions are:

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

The physical coordinate kernels depend only on transverse coordinates. The Lie derivatives act only on the Wilson-line/color coefficient structures and on commutator-correction tensors.

---

## Hard Constraints

Do **not**:

- modify production evolution code;
- train score or Hessian-score models;
- optimize large lattices before correctness is established;
- remove or weaken the finite-difference reference backend;
- invent generator-action signs or adjoint-index orientations;
- silently symmetrize ordered Lie derivatives;
- silently discard complex components;
- change physical kernel formulas;
- claim regulator independence;
- claim physical positivity;
- mark an analytic sector complete unless it agrees with the finite-difference oracle and preserves density closure.

If a convention is ambiguous, derive it from the existing generator-action tests and finite-difference checks. Stop that sector rather than inserting a guessed sign.

---

## Execution Mode

If `.git` metadata is unavailable, continue in no-git mode.

If Git is available:

```bash
git checkout -b nlo-current-analytic-coefficient-derivatives
```

If the worktree is dirty, report the dirty files before editing.

Maintain:

```text
reports/nlo_current/file_manifest.md
```

with every created or modified file.

---

## Phase 0: Baseline Audit

Inspect:

```text
src/nlo_current/physical_nlo_current.py
src/nlo_current/coefficient_derivatives.py
src/nlo_current/nlo_current_skeleton.py
src/nlo_current/nlo_velocity_evaluator.py
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py
src/nlo_current/cubic_commutator_terms.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/physical_kernel_adapter.py
src/nlo_current/physical_cubic_conventions.py
src/nlo_current/su3_adjoint.py
src/nlo_current/dipole_observable.py

tests/nlo_current/test_coefficient_derivatives.py
tests/nlo_current/test_physical_nlo_current.py
tests/nlo_current/test_physical_density_closure.py
tests/nlo_current/test_dipole_generator_actions.py

docs/nlo_current/coefficient_derivative_strategy.md
docs/nlo_current/coefficient_derivative_backend_limitations.md
docs/nlo_current/physical_nlo_current_assembly_plan.md
docs/nlo_current/physical_density_closure_derivation.md

reports/nlo_current/coefficient_derivative_backend_report.md
reports/nlo_current/physical_nlo_current_assembly_report.md
reports/nlo_current/physical_density_closure_report.md
reports/nlo_current/file_manifest.md
```

Create:

```text
reports/nlo_current/analytic_coefficient_derivative_start_status.md
```

Record:

- full test count;
- current finite-difference runtime;
- current closure residuals;
- current dtype conventions;
- current physical finite-grid policy;
- which derivative contractions are already available;
- that this workflow is non-production.

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

Expected baseline is approximately:

```text
129 passed
```

If all tests pass but the count differs, record the actual result and continue.

---

## Phase 1: Derivative Convention Note

Create:

```text
docs/nlo_current/analytic_lie_derivative_conventions.md
```

Document the exact left-generator convention used by the code.

Start from the already-tested definition:

\[
L_x^aF(U)
=
\left.
\frac{d}{d\epsilon}
F(e^{i\epsilon t^a}U_x)
\right|_{\epsilon=0}.
\]

Derive and verify, in the code convention:

\[
L_x^aU_y
=
i\,\delta_{xy}\,t^aU_y,
\]

\[
L_x^aU_y^\dagger
=
-i\,\delta_{xy}\,U_y^\dagger t^a.
\]

For the adjoint Wilson line, derive the exact index/sign rule from:

```text
src/nlo_current/su3_adjoint.py
tests/nlo_current/test_dipole_generator_actions.py
```

Do not assume the adjoint derivative sign from memory.

Document:

- the adjoint generator normalization;
- whether the derivative acts on the first or second adjoint index;
- left/right generator conversion;
- same-site commutator:
  \[
  [L_x^a,L_x^b]=f^{abc}L_x^c;
  \]
- distinct-site commutativity;
- ordered second derivatives;
- raw complex cubic kernel versus KLM-normalized real coefficients.

Add finite-difference convention-calibration tests before using the analytic rules elsewhere.

---

## Phase 2: Local Derivative Primitive Library

Create:

```text
src/nlo_current/analytic_lie_derivatives.py
```

Implement small, explicit primitives.

Suggested interfaces:

```python
def left_derivative_fundamental(U, site, color, target_site):
    """
    Return L_site^color U_target_site.
    """

def left_derivative_fundamental_dagger(U, site, color, target_site):
    """
    Return L_site^color U_target_site^\dagger.
    """

def left_derivative_adjoint(SA, site, color, target_site):
    """
    Return L_site^color S_A(target_site), with exact tested index convention.
    """

def second_left_derivative_adjoint(
    SA,
    first_site,
    first_color,
    second_site,
    second_color,
    target_site,
):
    """
    Return ordered L_first L_second S_A(target_site).
    """

def left_derivative_trace_word(...):
    """
    Product-rule derivative of an explicit fundamental trace word.
    """

def second_left_derivative_trace_word(...):
    """
    Ordered product-rule second derivative of a trace word.
    """
```

Requirements:

- explicit Kronecker locality;
- no full-configuration finite differences;
- preserve ordered derivative indices;
- support complex intermediate values;
- expose dtype;
- no silent casting;
- include docstrings with formulas.

Tests must compare every primitive against the existing finite-difference generator backend.

---

## Phase 3: Analytic Derivatives for Two-Generator Sectors

Create:

```text
src/nlo_current/analytic_two_generator_derivatives.py
```

Implement analytic/local derivatives for the physical \(K_2\) coefficients.

### 3.1 \(K_{JSJ}\)

The physical coordinate kernel is independent of Wilson lines. Differentiate only:

\[
(S_x-S_z)(S_y-S_z)
\]

or the exact coefficient representation used in the current assembly.

Implement:

```python
def analytic_dK2_KJSJ(...):
    """
    Return dK2^A = L_B K2^{AB} for K_JSJ.
    """
```

Check whether the expected coefficient-divergence cancellation holds in the project convention. Do not assume it; verify against finite differences.

### 3.2 \(K_{JSSJ}\)

Differentiate the exact ordered coefficient:

\[
A_{JSSJ}^{ab}
=
\int_{z,z'}
K_{JSSJ}
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}).
\]

Implement:

```python
def analytic_dK2_KJSSJ(...):
    """
    Return the full contracted L_B K2^{AB}.
    """
```

Preserve:

- \(z,z'\) ordering;
- adjoint index ordering;
- current basis conversion;
- finite-grid weights.

### 3.3 \(K_{q\bar q}\)

Differentiate the exact WORKNLO generator trace-product coefficient and the subtraction piece separately.

Implement:

```python
def analytic_dK2_Kqbarq_trace(...)
def analytic_dK2_Kqbarq_subtraction(...)
def analytic_dK2_Kqbarq(...)
```

Do not revert to the previously incorrect compact reduced trace-current expression.

Required diagnostics:

- trace contribution agreement;
- subtraction contribution agreement;
- full sum agreement;
- diagonal \(z'=z\) identity;
- finite-grid policy consistency.

---

## Phase 4: Analytic First Derivatives for Cubic Sectors

Create:

```text
src/nlo_current/analytic_cubic_derivatives.py
```

Implement the first derivatives:

\[
(LC\_K3)^{AB}=L_CK_3^{ABC},
\]

\[
(LB\_K3)^{AC}=L_BK_3^{ABC}.
\]

Work sector by sector.

### 4.1 \(K_{JJSJ}\)

Differentiate the exact coefficient blocks:

\[
A_{LLR}^{dea}
=
\int_z
K_{JJSJ}
f^{bde}S_z^{ba},
\]

\[
B_{LRR}^{ade}
=
-\int_z
K_{JJSJ}
f^{bde}S_z^{ab},
\]

plus virtual and commutator-correction contributions.

### 4.2 \(K_{JJSSJ}\)

Differentiate:

\[
A_{LLR}^{dea}
=
\int_{z,z'}
K_{JJSSJ}
f^{acb}S_z^{dc}S_{z'}^{eb},
\]

\[
B_{LRR}^{ade}
=
-\int_{z,z'}
K_{JJSSJ}
f^{acb}S_z^{cd}S_{z'}^{be},
\]

plus virtual, \(\widetilde K\)-related assembly pieces where they enter the normal-form tensors, and commutator corrections.

Requirements:

- raw physical cubic kernels may be complex;
- KLM-normalized normal-form coefficients should follow the existing explicit convention;
- differentiate the normalized coefficient used by the skeleton;
- no duplicate \((-i)\) factor;
- preserve derivative ordering;
- preserve same-site commutators.

---

## Phase 5: Analytic Second Derivatives for Cubic Sectors

Extend:

```text
src/nlo_current/analytic_cubic_derivatives.py
```

to implement:

\[
d2K3^A=L_BL_CK_3^{ABC}.
\]

Suggested staged approach:

1. Differentiate the analytic first-derivative representation.
2. Apply ordered product rules locally.
3. Handle:
   - derivatives acting twice on one adjoint Wilson line;
   - derivatives acting on different adjoint Wilson lines;
   - same-site noncommutativity;
   - coincident-site commutator correction tensors;
   - virtual terms.

Implement:

```python
def analytic_d2K3_KJJSJ(...)
def analytic_d2K3_KJJSSJ(...)
```

Do not obtain the “analytic” second derivative by finite-differencing the first derivative. A temporary hybrid diagnostic is allowed, but it must be named:

```text
hybrid_local_fd
```

and not marked analytic.

---

## Phase 6: Unified Backend API

Create or update:

```text
src/nlo_current/physical_coefficient_derivatives.py
```

Implement:

```python
def compute_physical_coefficient_derivatives(
    U,
    physical_terms,
    *,
    backend="analytic",
    sector_mask=None,
    fd_eps=None,
    return_by_sector=False,
):
    """
    Supported backends:
      - analytic
      - finite_difference
      - hybrid_local_fd
      - diagnostic
    """
```

Return a structured result:

```python
@dataclass
class PhysicalCoefficientDerivatives:
    dK2: np.ndarray
    LC_K3: np.ndarray
    LB_K3: np.ndarray
    d2K3: np.ndarray
    by_sector: dict
    backend: str
    metadata: dict
```

Metadata must include:

- implemented analytic sectors;
- pending sectors;
- fallback use;
- finite-difference step if any;
- dtype;
- cubic normalization convention;
- coordinate policy;
- evaluation counts.

Do not silently fall back from `analytic` to finite difference. Raise or mark the sector pending unless the caller explicitly selects a mixed backend.

---

## Phase 7: Sector-by-Sector Oracle Comparison

Create:

```text
scripts/nlo_current/check_analytic_coefficient_derivatives.py
```

Use:

- deterministic SU(3) configurations;
- the smallest valid physical coordinate set;
- at least two random seeds;
- multiple finite-difference reference steps.

For each sector report:

\[
R_{dK2}
=
\frac{\|dK2_{\rm analytic}-dK2_{\rm FD}\|}
{\|dK2_{\rm FD}\|+\epsilon},
\]

\[
R_{LC}
=
\frac{\|LC_{\rm analytic}-LC_{\rm FD}\|}
{\|LC_{\rm FD}\|+\epsilon},
\]

\[
R_{LB}
=
\frac{\|LB_{\rm analytic}-LB_{\rm FD}\|}
{\|LB_{\rm FD}\|+\epsilon},
\]

\[
R_{d2}
=
\frac{\|d2_{\rm analytic}-d2_{\rm FD}\|}
{\|d2_{\rm FD}\|+\epsilon}.
\]

Also report:

- maximum absolute residual;
- relative residual;
- imaginary residual;
- runtime;
- derivative evaluation count;
- memory use if readily available.

Do not compare a structurally zero tensor only by relative residual. Report absolute residual and expected-zero status.

---

## Phase 8: Required Tests

Create:

```text
tests/nlo_current/test_analytic_lie_derivatives.py
tests/nlo_current/test_analytic_coefficient_derivatives.py
```

Required tests:

### Primitive tests

1. \(L U\) agrees with finite differences.
2. \(L U^\dagger\) agrees with finite differences.
3. \(L S_A\) agrees with finite differences.
4. Ordered \(LL S_A\) agrees with finite differences.
5. Trace-word first derivatives agree.
6. Trace-word ordered second derivatives agree.
7. Same-site commutator identity passes.
8. Distinct-site derivatives commute.
9. No silent complex cast.

### Two-generator tests

10. \(K_{JSJ}\) analytic \(dK2\) agrees with FD.
11. \(K_{JSSJ}\) analytic \(dK2\) agrees with FD.
12. \(K_{q\bar q}\) trace derivative agrees with FD.
13. \(K_{q\bar q}\) subtraction derivative agrees with FD.
14. Full \(K_{q\bar q}\) derivative agrees with FD.

### Cubic tests

15. \(K_{JJSJ}\) analytic \(LC_K3\) agrees with FD.
16. \(K_{JJSJ}\) analytic \(LB_K3\) agrees with FD.
17. \(K_{JJSJ}\) analytic \(d2K3\) agrees with FD.
18. \(K_{JJSSJ}\) analytic \(LC_K3\) agrees with FD.
19. \(K_{JJSSJ}\) analytic \(LB_K3\) agrees with FD.
20. \(K_{JJSSJ}\) analytic \(d2K3\) agrees with FD.
21. Commutator-correction derivative contributions agree.
22. Raw versus KLM-normalized cubic convention is handled exactly once.
23. No `ComplexWarning`.
24. No imaginary information is silently lost.

### Backend tests

25. `backend="analytic"` uses no global finite-difference perturbations.
26. `backend="finite_difference"` preserves the reference behavior.
27. `backend="hybrid_local_fd"` is explicitly labeled.
28. No silent backend fallback.
29. Per-sector masks return correct zero/nonzero structures.

---

## Phase 9: Revalidate NLO Velocity and Density Closure

Update the physical NLO velocity path to accept:

```python
derivative_backend="analytic"
```

Run the existing density closure using analytic derivatives.

Required comparisons:

1. Analytic-derivative velocity versus FD-reference velocity.
2. Per-sector velocity differences.
3. Full projected density closure.
4. One full all-sector closure if runtime permits.
5. Constant-density limit.
6. \(K_3=0\) limit.
7. Cubic sparse closure with nonzero Hessian-score.
8. Toggle checks remain meaningful.

Target:

- preserve the established closure accuracy window;
- do not set a universal tolerance before inspecting scale;
- document both absolute and relative residuals.

If analytic derivatives disagree but density closure appears to pass because of cancellation, treat the backend as failed.

---

## Phase 10: Performance Benchmark

Create:

```text
scripts/nlo_current/benchmark_coefficient_derivative_backends.py
```

Benchmark:

- finite difference;
- analytic;
- hybrid local FD, if present.

Use at least:

- smallest valid setup;
- one slightly larger diagnostic setup if runtime permits.

Record:

- wall time;
- number of coefficient evaluations;
- speedup;
- peak memory if readily available;
- output residual against FD reference.

Create:

```text
reports/nlo_current/analytic_coefficient_derivative_benchmark.md
```

Do not claim asymptotic production scaling from two tiny cases. Report observed diagnostic speedup only.

---

## Phase 11: Documentation and Reports

Create:

```text
docs/nlo_current/analytic_physical_coefficient_derivative_derivation.md
reports/nlo_current/analytic_coefficient_derivative_validation_report.md
reports/nlo_current/analytic_coefficient_derivative_failure_modes.md
```

The derivation document must include:

- local generator rules;
- sector-by-sector derivative formulas;
- ordered second derivatives;
- commutator handling;
- cubic normalization handling;
- source-code mapping.

The validation report must include:

- implemented sectors;
- pending sectors;
- FD reference steps;
- residual tables;
- closure results;
- runtime comparisons;
- dtype checks;
- caveats.

The failure-mode report must include:

- wrong adjoint-index orientation;
- wrong generator sign;
- derivative-order reversal;
- accidental Hessian symmetrization;
- duplicate or missing \((-i)\) normalization;
- omitted commutator correction;
- finite-grid policy mismatch;
- silent backend fallback;
- zero-tensor relative-error instability;
- complex-to-real cast.

Update:

```text
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/physical_kernel_status.md
reports/nlo_current/physical_nlo_current_assembly_report.md
reports/nlo_current/physical_density_closure_report.md
reports/nlo_current/file_manifest.md
```

---

## Phase 12: Acceptance Criteria

This workflow is complete only if:

1. Analytic local generator primitives pass FD checks.
2. All implemented two-generator \(dK2\) derivatives pass FD checks.
3. All implemented cubic first derivatives pass FD checks.
4. All implemented cubic second derivatives pass FD checks.
5. Same-site ordered derivative conventions are preserved.
6. Commutator-correction derivatives are included.
7. Cubic normalization is applied exactly once.
8. No silent complex cast occurs.
9. `backend="analytic"` performs no global finite-difference perturbations.
10. Analytic and FD velocities agree on controlled tests.
11. Density-side closure still passes with the analytic backend.
12. Runtime improvement is measured.
13. The finite-difference backend remains available as a reference.
14. The full suite passes:

```bash
python3 -m pytest tests/nlo_current -q
```

15. No production evolution is implemented.
16. No score/Hessian-score training is implemented.
17. No physical positivity or regulator-independence claim is made.

If some sector cannot be completed analytically, leave it explicitly pending and provide a mixed-backend report. Do not mark the full analytic backend complete.

---

## Final Codex Response Required

At the end, report:

1. Files created and modified.
2. Test commands and results.
3. Exact analytic derivative primitives implemented.
4. Two-generator sectors completed.
5. Cubic sectors completed.
6. Any pending analytic sector.
7. FD reference step range.
8. Residuals for:
   - \(dK2\);
   - \(LC_K3\);
   - \(LB_K3\);
   - \(d2K3\).
9. Velocity agreement.
10. Density-closure agreement.
11. Runtime and observed speedup.
12. Whether any global FD fallback remains.
13. Whether any complex warning or imaginary loss remains.
14. Whether the analytic backend is fully justified or only partially complete.
15. Remaining blockers.

Recommended next stage only after this workflow passes:

```text
matrix-free contracted Hessian-score evaluation
```

followed by:

```text
physical-kernel Pawula/positivity diagnostics
```

and then:

```text
score/Hessian-score estimator design
```
