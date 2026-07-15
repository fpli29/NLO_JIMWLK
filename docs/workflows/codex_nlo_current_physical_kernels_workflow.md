# Codex Workflow: Non-Production Physical KLM Kernel Integration for NLO JIMWLK Current

## Purpose

Continue from the completed full Appendix A dipole validation workflow.

Current validated state:

1. All five NLO JIMWLK Hamiltonian sectors have passed Appendix A dipole validation:

   \[
   K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q},\quad K_{JJSJ},\quad K_{JJSSJ}.
   \]

2. The most delicate structures are validated:
   - \(K_{JJSSJ}\) \(\widetilde K\) contribution;
   - pure eight-kernel contribution;
   - cubic \(1/3\) virtual factors;
   - cubic \((-i)\) convention under the current Hermitian-generator direct-action code.

3. The non-production NLO current skeleton exists:

   \[
   \partial_YW
   =
   -L_A(K_1^AW)
   +
   \frac12L_AL_B(K_2^{AB}W)
   -
   \frac16L_AL_BL_C(K_3^{ABC}W).
   \]

4. The velocity evaluator supports score and Hessian-score:

   \[
   s_A=L_A\log W,
   \qquad
   H_{AB}=L_As_B.
   \]

The current goal is to add **non-production physical KLM kernel formulas** for small dense coordinate diagnostics.

This workflow should implement coordinate-space kernel functions and sanity tests only.

Do **not** implement production evolution.  
Do **not** optimize for large lattices.  
Do **not** train score/Hessian-score models.  
Do **not** replace the validated synthetic-kernel tests.  
Do **not** silently choose a regulator scheme as production-ready.

The output should be a small diagnostic physical-kernel module that can later feed:

```python
assemble_nlo_current_terms(...)
```

and the existing dipole/direct-action validation scripts.

---

## Execution Mode

This workspace may not be a Git repository. If `.git` metadata is unavailable, continue in no-git mode.

If no-git mode is used:

1. Do not modify existing production files unless absolutely necessary.
2. Prefer adding isolated files under:
   - `docs/nlo_current/`
   - `src/nlo_current/`
   - `tests/nlo_current/`
   - `scripts/nlo_current/`
   - `reports/nlo_current/`
3. Before editing an existing file, report the file and why.
4. Maintain/update:

```text
reports/nlo_current/file_manifest.md
```

listing every created or modified file.

If this is a Git repository, create a branch:

```bash
git checkout -b nlo-current-physical-kernels-diagnostic
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect:

```text
references/WORKNLO.tex
references/1405.0418v2.pdf

docs/nlo_current/dipole_validation_status.md
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/KLM_appendix_A_dipole_targets_notes.md
docs/nlo_current/two_generator_sector_summary.md
docs/nlo_current/three_generator_sector_summary.md

src/nlo_current/dipole_appendix_targets.py
src/nlo_current/dipole_hamiltonian_action.py
src/nlo_current/nlo_current_skeleton.py
src/nlo_current/synthetic_kernels.py
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py

reports/nlo_current/full_dipole_validation_report.md
reports/nlo_current/kjjssj_appendix_target_validation_report.md
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/physical_kernel_start_status.md
```

It should state:

- whether all five Appendix A dipole targets are marked passed;
- whether previous tests still pass;
- whether no-git mode is active;
- whether `WORKNLO.tex` is available;
- that this workflow is non-production kernel diagnostics only.

Run before changes:

```bash
python3 -m pytest tests/nlo_current -q
```

Expected current baseline is around:

```text
84 passed
```

If the number differs, report it, but continue if all tests pass.

---

## Phase 1: Documentation — Physical Kernel Integration Plan

Create:

```text
docs/nlo_current/physical_kernel_integration_plan.md
```

This document must include:

### 1. Goal

Implement non-production coordinate-space KLM kernel functions for diagnostic small lattices:

\[
K_{JSJ},
\quad
K_{JSSJ},
\quad
K_{q\bar q},
\quad
K_{JJSJ},
\quad
K_{JJSSJ}.
\]

### 2. Source of truth

Use `references/WORKNLO.tex` as primary source.

Use the PDF only as a cross-check.

Do not infer missing formulas. If a formula is not confidently found, implement a stub that raises `NotImplementedError` and document it.

### 3. Kernel types

Distinguish clearly:

- unbarred singlet kernels for dipole validation;
- barred/generalized kernels for nonsinglet/full configuration-level Hamiltonian.

This workflow should implement unbarred kernels first. If barred kernels are found and simple to add, add them only with explicit naming:

```python
KJSJ_unbarred(...)
KJSJ_barred(...)
```

Do not silently use barred kernels where unbarred kernels are expected.

### 4. Coordinate conventions

Use two-dimensional transverse coordinates:

```python
coords: np.ndarray, shape (Nsite, 2)
```

Define:

\[
X = x-z,
\qquad
Y = y-z,
\qquad
X' = x-z',
\qquad
Y' = y-z',
\]

or whatever notation KLM uses. Document the mapping exactly.

### 5. Singularity handling

KLM coordinate kernels have singular denominators when points coincide.

For this diagnostic workflow, implement a clear policy:

```python
singularity_policy = "raise" | "nan" | "eps"
```

Default should be:

```python
"raise"
```

for exact diagnostics.

If `"eps"` is used, require an explicit `eps` argument and record it in metadata.

Do not present any regulator as physically final.

### 6. Expected symmetries

Document and test known symmetries from KLM:

- \(K_{JSSJ}\) and \(K_{q\bar q}\) symmetric under appropriate \(z\leftrightarrow z'\) and/or \(x\leftrightarrow y\) exchanges.
- \(K_{JJSSJ}\) antisymmetric under simultaneous:
  \[
  x\leftrightarrow y,
  \qquad
  z\leftrightarrow z'.
  \]
- \(K_{JJSJ}\) has the antisymmetry convention already tested in synthetic diagnostics, if confirmed from source.

If a symmetry is uncertain, mark it uncertain and write a test as `xfail` or report-only, not as a hard pass.

---

## Phase 2: Extract Kernel Formula Notes

Create:

```text
docs/nlo_current/KLM_physical_kernel_formula_notes.md
```

This file should transcribe exact kernel formulas from `WORKNLO.tex`.

For each kernel, include:

```text
Kernel:
Formula:
WORKNLO.tex line range:
Arguments:
Symmetries:
Singular denominators:
Barred/unbarred status:
Implementation status:
Open questions:
```

Required kernels:

1. `K_JSJ`
2. `K_JSSJ`
3. `K_qbarq`
4. `K_JJSJ`
5. `K_JJSSJ`
6. `tilde_K` combination if used for diagnostics

Use `grep`/search to find definitions. Suggested search terms:

```bash
grep -n "K_{JSJ}\|K_{JSSJ}\|K_{q\\\\bar{q}}\|K_{JJSJ}\|K_{JJSSJ}\|tilde K\|\\bar K" references/WORKNLO.tex
```

Do not implement formulas before writing this notes file.

If a kernel formula is not found or appears in multiple conventions, stop that kernel implementation and mark it pending.

---

## Phase 3: Implement Coordinate Utilities

Create:

```text
src/nlo_current/coordinate_kernels.py
```

Start with coordinate utilities:

```python
import numpy as np

def vec(coords, i, j):
    """Return coords[i] - coords[j]."""

def norm2(v):
    """Return squared Euclidean norm."""

def dot(a, b):
    """Return Euclidean dot product."""

def cross2(a, b):
    """
    Return 2D scalar cross product:
        a_x b_y - a_y b_x
    if needed by kernel formulas.
    """

def safe_inv(x, singularity_policy="raise", eps=None, name="denominator"):
    """
    Return 1/x with explicit singularity policy.
    """

def validate_coords(coords):
    """
    Validate shape (Nsite,2), finite values.
    """

def pairwise_dist2(coords):
    """
    Return matrix r_ij^2.
    """
```

Tests should cover:

- shape validation;
- zero-distance detection;
- safe inverse behavior for `"raise"`, `"nan"`, and `"eps"` policies.

---

## Phase 4: Implement Physical Kernel Module

Create:

```text
src/nlo_current/physical_kernels.py
```

Implement dense diagnostic kernel builders.

Required interface:

```python
def KJSJ_unbarred_value(coords, x, y, z, *, Nc=3, nf=0, alpha_s=1.0,
                        singularity_policy="raise", eps=None):
    """
    Return K_JSJ(x,y;z) from the KLM formula.
    """

def KJSSJ_unbarred_value(coords, x, y, z, zp, *, Nc=3, nf=0, alpha_s=1.0,
                         singularity_policy="raise", eps=None):
    """
    Return K_JSSJ(x,y;z,z').
    """

def Kqbarq_unbarred_value(coords, x, y, z, zp, *, Nc=3, nf=0, alpha_s=1.0,
                          singularity_policy="raise", eps=None):
    """
    Return K_qbarq(x,y;z,z').
    """

def KJJSJ_unbarred_value(coords, w, x, y, z, *, Nc=3, nf=0, alpha_s=1.0,
                         singularity_policy="raise", eps=None):
    """
    Return K_JJSJ(w;x,y;z).
    """

def KJJSSJ_unbarred_value(coords, w, x, y, z, zp, *, Nc=3, nf=0, alpha_s=1.0,
                          singularity_policy="raise", eps=None):
    """
    Return K_JJSSJ(w;x,y;z,z').
    """
```

Also provide dense array builders:

```python
def build_KJSJ_unbarred(coords, **params):
    """Return array shape (N,N,N): K[x,y,z]."""

def build_KJSSJ_unbarred(coords, **params):
    """Return array shape (N,N,N,N): K[x,y,z,zp]."""

def build_Kqbarq_unbarred(coords, **params):
    """Return array shape (N,N,N,N)."""

def build_KJJSJ_unbarred(coords, **params):
    """Return array shape (N,N,N,N): K[w,x,y,z]."""

def build_KJJSSJ_unbarred(coords, **params):
    """Return array shape (N,N,N,N,N): K[w,x,y,z,zp]."""

def build_all_unbarred_physical_kernels(coords, **params):
    """
    Return dict:
      {
        "KJSJ": ...,
        "KJSSJ": ...,
        "Kqbarq": ...,
        "KJJSJ": ...,
        "KJJSSJ": ...,
      }
    plus metadata if useful.
    """
```

Important:

- If a formula is pending, the corresponding function should raise `NotImplementedError`.
- Do not silently return zeros for unimplemented kernels.
- Keep kernels real or complex exactly as formula requires.
- Preserve \(i\) factors if present.
- Expose `alpha_s`, `Nc`, `nf`, and any scheme constants only if present in source formulas.
- Record if formulas are LO-substitution-compatible, e.g. \(K_{JSJ}\to -M/2\), in docs not in code unless explicitly needed.

---

## Phase 5: Kernel Formula Tests

Create:

```text
tests/nlo_current/test_physical_kernels.py
```

Tests:

### Test 1: coordinate utility behavior

Verify coordinate utilities and singularity policies.

### Test 2: dense builder shapes

For a non-degenerate coordinate set with \(N=4\), verify shapes:

```text
KJSJ:   (N,N,N)
KJSSJ:  (N,N,N,N)
Kqbarq: (N,N,N,N)
KJJSJ:  (N,N,N,N)
KJJSSJ: (N,N,N,N,N)
```

Only run this test for kernels that are implemented. Pending kernels should be explicitly skipped with a message.

### Test 3: finite values away from singularities

For a coordinate set with all relevant points distinct, implemented kernels should be finite or explicitly complex-finite.

### Test 4: singularity policy

For coincident coordinate input, default `"raise"` should raise a clear error. `"nan"` should produce `np.nan`. `"eps"` should produce finite values with metadata or clear docstring.

### Test 5: symmetry checks

For implemented kernels with known symmetries:

- \(K_{JSSJ}\) symmetry as documented.
- \(K_{q\bar q}\) symmetry as documented.
- \(K_{JJSSJ}\) KLM simultaneous antisymmetry:
  \[
  K(w;x,y;z,z')=-K(w;y,x;z',z).
  \]
- \(K_{JJSJ}\) antisymmetry if confirmed.

Do not hard-code uncertain symmetries as tests. If uncertain, write report-only diagnostics.

### Test 6: tilde-K consistency if implemented

If `tilde_K` helper exists, compare direct helper output to the definition from notes.

---

## Phase 6: Integration with Existing Skeleton

Create:

```text
src/nlo_current/physical_kernel_adapter.py
```

Implement:

```python
def physical_kernels_for_skeleton(coords, *, Nc=3, nf=0, alpha_s=1.0,
                                  singularity_policy="raise", eps=None):
    """
    Build physical unbarred kernels and return a dict compatible with:
        assemble_nlo_current_terms(...)
    """

def physical_kernel_metadata(coords, params):
    """
    Return metadata:
      - coordinate count
      - singularity policy
      - implemented kernels
      - pending kernels
      - parameter values
    """
```

Create test:

```text
tests/nlo_current/test_physical_kernel_adapter.py
```

Tests:

1. Adapter returns expected keys for implemented kernels.
2. Metadata records pending kernels if any.
3. Adapter output can be passed to `assemble_nlo_current_terms(...)` in metadata-only mode.
4. If all kernels are implemented, adapter output can be passed to dense assembly for \(N=2\) or \(3\) using random Wilson lines.

Do not require full dense NLO current assembly if kernel singularities make the tiny grid invalid. Use a non-degenerate coordinate set and/or metadata-only mode.

---

## Phase 7: Diagnostic Script

Create:

```text
scripts/nlo_current/check_physical_kernel_integration.py
```

It should:

1. Build a small non-degenerate coordinate set, e.g. \(N=4\):

```python
coords = np.array([
    [0.0, 0.0],
    [1.0, 0.2],
    [0.3, 1.1],
    [1.4, 1.3],
])
```

2. Build all implemented physical kernels.
3. Compute and report:
   - shapes;
   - norms;
   - finite/nonfinite counts;
   - known symmetry residuals;
   - singularity policy used;
   - implemented vs pending kernels.

4. If possible, pass kernels to existing non-production skeleton adapter or `assemble_nlo_current_terms(...)` in metadata-only mode.

5. Save:

```text
reports/nlo_current/physical_kernel_integration_report.md
```

The report must include:

- source formula notes path;
- random seed if random data used;
- coordinates used;
- parameter values;
- per-kernel implementation status;
- per-kernel norm;
- per-kernel symmetry residuals;
- warnings about non-production status and singularity handling.

---

## Phase 8: Update Status Documentation

Create:

```text
docs/nlo_current/physical_kernel_status.md
```

Include:

1. Which physical kernels are implemented.
2. Which formulas were sourced from `WORKNLO.tex` line ranges.
3. Which kernels remain pending, if any.
4. Symmetry checks passed.
5. Singularity policy.
6. How to call the diagnostic builder.
7. What remains before production use:
   - proper UV/IR regularization or subtraction strategy;
   - barred/nonsinglet kernels if needed;
   - analytic coefficient derivatives for physical kernels;
   - performance/sparse locality;
   - physical validation beyond synthetic coordinate tests.

Also update:

```text
docs/nlo_current/nlo_current_map_status.md
```

Add a section:

```text
Physical kernel integration status
```

---

## Phase 9: Optional Dipole Recheck with Physical Kernels

If all five kernels are implemented and a non-singular coordinate set can be chosen, add a script section to:

```text
scripts/nlo_current/full_dipole_validation.py
```

or create:

```text
scripts/nlo_current/full_dipole_validation_physical_kernels.py
```

This should:

1. Use physical kernels instead of synthetic kernels.
2. Run the already-passed Appendix A direct-vs-target comparisons.
3. Save:

```text
reports/nlo_current/full_dipole_validation_physical_kernel_report.md
```

If singularities or pending kernels prevent this, skip and document why.

Do not force this phase if formulas are pending.

---

## Phase 10: Acceptance Criteria

Stop when all are true:

1. Formula notes exist:
   ```text
   docs/nlo_current/KLM_physical_kernel_formula_notes.md
   ```

2. Coordinate utilities exist:
   ```text
   src/nlo_current/coordinate_kernels.py
   ```

3. Physical kernel module exists:
   ```text
   src/nlo_current/physical_kernels.py
   ```

4. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

5. Diagnostic report exists:
   ```text
   reports/nlo_current/physical_kernel_integration_report.md
   ```

6. Status doc exists:
   ```text
   docs/nlo_current/physical_kernel_status.md
   ```

7. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

8. No production evolution code modified.
9. No score/Hessian-score model training implemented.
10. Any unimplemented formula is explicitly marked pending and not silently replaced by zeros.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Which physical kernels were implemented:
   - \(K_{JSJ}\)
   - \(K_{JSSJ}\)
   - \(K_{q\bar q}\)
   - \(K_{JJSJ}\)
   - \(K_{JJSSJ}\)
4. Which formulas remain pending, if any.
5. WORKNLO.tex line ranges used for each implemented formula.
6. Symmetry residuals for implemented kernels.
7. Singularity policy behavior.
8. Whether physical kernels can be passed into the non-production skeleton.
9. Whether optional physical-kernel dipole recheck ran.
10. Remaining blockers before production:
    - regulator/subtraction strategy;
    - barred/nonsinglet kernels;
    - analytic coefficient derivatives;
    - score/Hessian-score estimator design;
    - performance and sparse/local implementation.

Do not claim production readiness. This workflow integrates physical kernels for small dense diagnostics only.
