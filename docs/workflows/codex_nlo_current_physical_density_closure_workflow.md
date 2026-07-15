# Codex Workflow: Physical NLO JIMWLK Density-Side Closure Validation

## Purpose

Continue from the completed non-production physical NLO generalized-current assembly.

Current validated state:

1. All five Appendix A dipole sectors pass:
   \[
   K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q},\quad K_{JJSJ},\quad K_{JJSSJ}.
   \]

2. All five unbarred physical kernel interfaces are implemented and physically normalized.

3. The cubic convention is explicit:
   \[
   K_{\rm raw}^{\rm cubic}
   \longrightarrow
   (-i)K_{\rm raw}^{\rm cubic}
   \longrightarrow
   K_{\rm KLM\text{-}normalized}^{\rm real}.
   \]

4. The physical generalized-current assembly exists:
   - `assemble_physical_K1(...)`
   - `assemble_physical_K2(...)`
   - `assemble_physical_K3(...)`
   - `assemble_physical_terms(...)`
   - `compute_physical_coefficient_derivatives(...)`
   - `evaluate_physical_nlo_velocity(...)`

5. A diagnostic finite-difference coefficient-derivative backend exists.

The goal is to verify the full density-side closure:

\[
\boxed{
-L_A\!\left(v^A[W]\,W\right)
=
-L_A(K_1^A W)
+\frac12L_AL_B(K_2^{AB}W)
-\frac16L_AL_BL_C(K_3^{ABC}W)
}
\]

for controlled positive test densities \(W\) on the smallest valid physical-kernel lattice.

This is a dense, non-production, end-to-end closure check. It should prove numerically that the assembled physical generalized current reproduces the physical NLO density operator.

---

## Hard Constraints

Do **not**:

- implement production evolution;
- train score or Hessian-score models;
- optimize for large lattices;
- replace the finite-difference backend with an unvalidated analytic backend;
- claim physical positivity or regulator independence;
- loosen tolerances merely to force agreement;
- silently symmetrize the ordered Hessian-score;
- silently discard complex components;
- change already-passed Appendix A or physical-kernel formulas unless a test-proven bug is found.

This workflow is a tiny-lattice closure validation only.

---

## Execution Mode

If `.git` metadata is unavailable, continue in no-git mode.

If Git is available:

```bash
git checkout -b nlo-current-density-closure-validation
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
src/nlo_current/physical_nlo_current.py
src/nlo_current/nlo_current_skeleton.py
src/nlo_current/nlo_velocity_evaluator.py
src/nlo_current/coefficient_derivatives.py
src/nlo_current/physical_kernel_adapter.py
src/nlo_current/physical_kernels.py
src/nlo_current/physical_cubic_conventions.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/cubic_commutator_terms.py
src/nlo_current/dipole_observable.py

tests/nlo_current/test_physical_nlo_current.py
tests/nlo_current/test_coefficient_derivatives.py
tests/nlo_current/test_dipole_generator_actions.py

docs/nlo_current/physical_nlo_current_assembly_plan.md
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md

reports/nlo_current/physical_nlo_current_assembly_report.md
reports/nlo_current/physical_kernel_integration_report.md
reports/nlo_current/physical_cubic_dtype_audit.md
reports/nlo_current/file_manifest.md
```

Create:

```text
reports/nlo_current/density_closure_start_status.md
```

Record:

- full test count;
- physical assembly status;
- derivative backend status;
- cubic dtype convention;
- finite-grid policy;
- that this workflow is non-production.

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

Expected baseline is approximately:

```text
117 passed
```

If all tests pass but the count differs, continue and record the actual result.

---

## Phase 1: Closure Derivation Note

Create:

```text
docs/nlo_current/physical_density_closure_derivation.md
```

Derive the equality between the direct normal-form density operator and the generalized current.

Start from:

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W).
\]

Define:

\[
J^A
=
K_1^AW
-\frac12L_B(K_2^{AB}W)
+\frac16L_BL_C(K_3^{ABC}W).
\]

Then:

\[
\partial_YW=-L_AJ^A.
\]

Using:

\[
s_A=L_A\log W,
\qquad
H_{BC}=L_Bs_C,
\]

derive:

\[
v^A=\frac{J^A}{W},
\]

and explicitly reproduce:

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
\right].
\]

State these caveats:

1. \(H_{BC}=L_Bs_C\) is ordered and must not be assumed symmetric.
2. Coincident-site commutator corrections must already be included in \(K_1,K_2,K_3\).
3. The identity is local where \(W>0\).
4. Algebraic closure does not imply positivity preservation.
5. Physical finite-grid prescriptions remain diagnostic.

---

## Phase 2: Positive Test Density Module

Create:

```text
src/nlo_current/test_densities.py
```

Implement positive, exactly differentiable test densities.

Required dataclass:

```python
@dataclass
class TestDensityResult:
    log_weight: float
    weight: float
    score: np.ndarray
    hessian_score: np.ndarray
    metadata: dict
```

Implement at least three densities.

### Density A: Single-link trace density

\[
\Phi_A(U)=\lambda_1\,\mathrm{Re\,tr}(U_x),
\qquad
W_A(U)=e^{\Phi_A(U)}.
\]

### Density B: Dipole trace density

\[
\Phi_B(U)=
\lambda_1\,\mathrm{Re\,tr}(U_x)
+
\lambda_2\,\mathrm{Re\,tr}(U_x^\dagger U_y),
\]

\[
W_B(U)=e^{\Phi_B(U)}.
\]

### Density C: Multi-link nonlinear density

\[
\Phi_C(U)=
\lambda_1\,\mathrm{Re\,tr}(U_x)
+
\lambda_2\,\mathrm{Re\,tr}(U_x^\dagger U_y)
+
\lambda_3\left[\mathrm{Re\,tr}(U_y^\dagger U_z)\right]^2.
\]

\[
W_C(U)=e^{\Phi_C(U)}.
\]

Requirements:

- return unnormalized \(W\); normalization cancels in score and closure checks;
- \(W>0\) by construction;
- score and ordered Hessian-score must come from the same \(\log W\);
- use analytic derivatives where straightforward;
- otherwise use a high-accuracy finite-difference reference;
- expose the derivative backend in metadata;
- do not use a learned model.

Implement helpers:

```python
evaluate_test_density(...)
compute_test_density_score(...)
compute_test_density_hessian_score(...)
```

---

## Phase 3: Direct Density Operator

Create:

```text
src/nlo_current/physical_density_operator.py
```

Implement:

```python
def evaluate_direct_density_operator(
    U,
    W_fn,
    physical_terms,
    *,
    derivative_backend="finite_difference",
    fd_eps=...,
    sector_mask=None,
):
    """
    Evaluate:
      -L_A(K1^A W)
      + 1/2 L_A L_B(K2^{AB} W)
      - 1/6 L_A L_B L_C(K3^{ABC} W)
    at U.
    """
```

Also implement:

```python
def evaluate_direct_density_operator_by_sector(...):
    """
    Return contributions from:
      KJSJ
      KJSSJ
      Kqbarq
      KJJSJ
      KJJSSJ
      commutator corrections
    """
```

Requirements:

- preserve written derivative ordering;
- use the same Lie derivative convention as existing generator tests;
- preserve physical finite-grid and cubic normalization metadata;
- do not silently cast complex arrays to real;
- return complex diagnostics if needed, together with an expected-real residual;
- support masks for:
  - \(K_1\) only;
  - \(K_2\) only;
  - \(K_3\) only;
  - individual physical sectors;
  - commutator corrections on/off.

Correctness first; do not optimize.

---

## Phase 4: Current-Divergence Operator

Create:

```text
src/nlo_current/physical_current_divergence.py
```

Implement:

```python
def evaluate_current_divergence(
    U,
    W_fn,
    score,
    hessian_score,
    physical_terms,
    coefficient_derivatives,
    *,
    derivative_backend="finite_difference",
    fd_eps=...,
    sector_mask=None,
):
    """
    1. Evaluate physical NLO velocity v^A.
    2. Evaluate -L_A(v^A W).
    """
```

Also implement:

```python
def evaluate_current_divergence_by_sector(...):
    """
    Return per-sector divergence contributions.
    """
```

Important requirements:

- \(v^A\) depends on \(U\); the outer derivative must act on both \(v^A(U)\) and \(W(U)\);
- do not freeze the velocity during the divergence calculation;
- if the outer finite difference perturbs \(U\), recompute:
  - score;
  - Hessian-score;
  - physical terms;
  - coefficient derivatives;
  consistently at the perturbed point;
- use the same backend and Lie-derivative convention on both sides;
- record all step sizes in metadata.

---

## Phase 5: Closure Comparator

Create:

```text
src/nlo_current/physical_density_closure.py
```

Implement:

```python
@dataclass
class ClosureResult:
    direct_value: complex
    current_value: complex
    abs_residual: float
    rel_residual: float
    direct_by_sector: dict
    current_by_sector: dict
    sector_residuals: dict
    metadata: dict
```

Main API:

```python
def compare_physical_density_closure(
    U,
    density,
    *,
    physical_policy,
    derivative_backend="finite_difference",
    fd_eps=...,
    sector_mask=None,
):
    """
    Compare the direct density operator and generalized-current divergence.
    """
```

Compute:

\[
R_{\rm abs}
=
\left|
{\cal G}_{\rm direct}[W]
-
{\cal G}_{\rm current}[W]
\right|,
\]

\[
R_{\rm rel}
=
\frac{
\left|
{\cal G}_{\rm direct}[W]
-
{\cal G}_{\rm current}[W]
\right|
}{
|{\cal G}_{\rm direct}[W]|
+
|{\cal G}_{\rm current}[W]|
+
\epsilon_{\rm floor}
}.
\]

Report real and imaginary residuals separately.

---

## Phase 6: Finite-Difference Convergence Study

Create:

```text
scripts/nlo_current/check_physical_density_closure.py
```

Use:

- the smallest valid non-degenerate physical coordinate set;
- a deterministic SU(3) Wilson-line configuration;
- Density A, Density B, and Density C.

Use at least four finite-difference steps, for example:

```text
1e-3
5e-4
2.5e-4
1.25e-4
```

or a stable range established by existing derivative tests.

For each density and step, record:

- direct operator;
- current divergence;
- absolute residual;
- relative residual;
- real residual;
- imaginary residual;
- per-sector residuals;
- wall time;
- derivative evaluation count.

Look for a stable convergence window. Do not insist on monotonic convergence down to roundoff.

---

## Phase 7: Required Toggle Diagnostics

The script and tests must include controlled failures.

### Toggle 1: Omit coefficient derivatives

Compare:

```text
with coefficient derivatives
without coefficient derivatives
```

The omitted version should differ generically.

### Toggle 2: Omit Hessian-score

Set:

\[
H_{BC}=0.
\]

The cubic closure should fail generically.

### Toggle 3: Omit commutator corrections

Disable:

\[
K_{2,\rm comm},
\qquad
K_{1,\rm comm}.
\]

Use a configuration where coincident-site words are exercised. Closure should worsen.

### Toggle 4: Remove cubic normalization

Use raw cubic coefficients without the established \((-i)\) normalization. Closure should fail and/or have the wrong complex character.

### Toggle 5: Set \(K_3=0\)

Verify reduction to the second-order formula.

### Toggle 6: Set \(K_2=K_3=0\)

Verify:

\[
-L_A(K_1^AW)
=
-L_A(v^AW),
\qquad
v^A=K_1^A.
\]

### Toggle 7: Constant density

For:

\[
W=1,
\qquad
s=0,
\qquad
H=0,
\]

verify that only the expected coefficient-derivative and \(K_1\) pieces remain.

---

## Phase 8: Tests

Create:

```text
tests/nlo_current/test_physical_density_closure.py
```

Required tests:

1. Positive test densities return finite positive weights.
2. Score finite-difference consistency.
3. Ordered Hessian-score finite-difference consistency.
4. Hessian is not silently symmetrized.
5. First-order-only closure.
6. Second-order closure with \(K_3=0\).
7. Cubic synthetic closure with nonzero Hessian-score.
8. Full physical closure on the smallest valid setup.
9. Sector-by-sector closure.
10. Omitting coefficient derivatives fails generically.
11. Omitting Hessian-score fails for cubic sectors.
12. Omitting cubic normalization fails.
13. Omitting commutator corrections worsens closure when exercised.
14. No `ComplexWarning`.
15. No silent imaginary-part loss.
16. Constant-density limit.
17. Finite-difference convergence has at least one stable accuracy window.

Keep runtime controlled. Use the smallest valid grid and deterministic seeds.

---

## Phase 9: Reports

Create:

```text
reports/nlo_current/physical_density_closure_report.md
```

Include:

1. Coordinate set.
2. Wilson-line configuration.
3. Physical finite-grid policy.
4. Cubic normalization convention.
5. Test densities and parameters.
6. Score/Hessian backend.
7. Coefficient-derivative backend.
8. Finite-difference step table.
9. Direct/current values.
10. Absolute and relative residuals.
11. Per-sector residuals.
12. Toggle-failure results.
13. Runtime.
14. Warnings and caveats.

Add this boxed conclusion only if supported:

\[
\boxed{
\text{The non-production physical NLO density operator is numerically reproduced by the generalized score/Hessian-score current on the tested tiny lattice.}
}
\]

Immediately follow with:

```text
This is a controlled dense diagnostic closure test. It is not a production
evolution result, a regulator-independence proof, or a positivity proof.
```

Also create:

```text
reports/nlo_current/physical_density_closure_failure_modes.md
```

Document:

- finite-difference cancellation;
- inconsistent score/Hessian backends;
- frozen-velocity outer derivative errors;
- missing coefficient derivatives;
- missing commutator corrections;
- raw-versus-normalized cubic convention errors;
- singular coordinate entries;
- dtype/casting errors.

---

## Phase 10: Documentation Updates

Update:

```text
docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md
```

Add:

```text
Physical density-side closure validation
```

Include:

\[
{\cal G}_{\rm direct}[W]
=
{\cal G}_{\rm current}[W].
\]

State the tested densities, backend, residuals, and caveats.

Also update:

```text
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/physical_kernel_status.md
reports/nlo_current/physical_nlo_current_assembly_report.md
reports/nlo_current/file_manifest.md
```

Do not claim production readiness.

---

## Phase 11: Acceptance Criteria

The workflow is complete only if:

1. The closure comparator exists.
2. At least one positive nontrivial density passes full physical closure.
3. At least one cubic sector passes with nonzero Hessian-score.
4. The \(K_3=0\) second-order limit passes.
5. The first-order-only limit passes.
6. Per-sector residuals are reported.
7. Omitting coefficient derivatives causes a meaningful failure.
8. Omitting Hessian-score causes a meaningful cubic failure.
9. Omitting cubic normalization causes a meaningful failure.
10. Omitting commutator corrections worsens closure when relevant.
11. No silent complex cast remains.
12. A finite-difference convergence window is documented.
13. The full suite passes:

```bash
python3 -m pytest tests/nlo_current -q
```

14. No production evolution was added.
15. No score/Hessian-score training was added.
16. No physical positivity claim was made.

If full closure does not pass, do not force it. Classify the mismatch by:

- sector;
- derivative order;
- coefficient derivative;
- commutator correction;
- outer divergence;
- dtype;
- finite-difference step;
- finite-grid policy.

---

## Final Codex Response Required

At the end, report:

1. Files created and modified.
2. Test commands and results.
3. Positive test densities used.
4. Coordinate and Wilson-line setup.
5. Finite-difference backend and step range.
6. Best full closure residual.
7. Per-sector residuals.
8. Whether coefficient-derivative omission fails.
9. Whether Hessian omission fails.
10. Whether commutator omission fails.
11. Whether raw cubic normalization fails.
12. Whether any complex warning or imaginary loss remains.
13. Runtime.
14. Whether the boxed closure conclusion is justified.
15. Remaining blockers.

Recommended next stage only after closure passes:

```text
analytic/local physical coefficient derivatives
```

followed by:

```text
score/Hessian-score estimator design
```
