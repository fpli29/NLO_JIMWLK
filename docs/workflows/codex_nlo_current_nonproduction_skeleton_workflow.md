# Codex Workflow: Assemble a Non-Production NLO JIMWLK Current Skeleton

## Purpose

Continue from the completed workflows:

1. Two-generator sector:
   \[
   K_{JSJ}: \text{LO-like score current},
   \]
   \[
   K_{JSSJ},\ K_{q\bar q}: \text{generic ordered }J_LAJ_R\text{ currents with coefficient drift}.
   \]

2. Cubic sector, distinct-site:
   \[
   K_{JJSJ},\ K_{JJSSJ}: \text{score + Hessian-score currents}.
   \]

3. Cubic coincident-site commutators:
   \[
   [L_x^a,L_x^b]=f^{abc}L_x^c
   \]
   were canonicalized and shown to generate lower-order corrections:
   \[
   K_{2,\rm comm},
   \qquad
   K_{1,\rm comm}.
   \]

The current goal is to assemble a **non-production skeleton** for the NLO generalized probability current.

This skeleton should organize the already validated pieces into a unified diagnostic interface:

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W),
\]

with current

\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W),
\]

and velocity

\[
v^A
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
K_3^{ABC}(H_{BC}+s_Bs_C)
\right],
\]

where

\[
s_A=L_A\log W,
\qquad
H_{BC}=L_Bs_C.
\]

Do **not** implement production evolution.  
Do **not** train score or Hessian-score models.  
Do **not** optimize for large lattices.  
Do **not** claim the NLO flow is complete for production use.  
This is a dense small-lattice diagnostic skeleton only.

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
git checkout -b nlo-current-nonproduction-skeleton
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect the existing artifacts from previous workflows:

```text
docs/nlo_current/two_generator_sector_summary.md
docs/nlo_current/three_generator_sector_summary.md
docs/nlo_current/cubic_current_with_commutator_corrections.md
docs/nlo_current/cubic_coincident_site_commutators.md

src/nlo_current/su3_adjoint.py
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py
src/nlo_current/finite_difference_scores.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/cubic_commutator_terms.py

reports/nlo_current/kjssj_symmetry_report.md
reports/nlo_current/kqbarq_symmetry_report.md
reports/nlo_current/kjjsj_cubic_requirements_report.md
reports/nlo_current/kjjssj_cubic_requirements_report.md
reports/nlo_current/cubic_commutator_corrections_report.md
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/nlo_current_skeleton_start_status.md
```

It should state:

- whether previous workflow files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- whether the commutator workflow completed;
- that this workflow is non-production only.

Run before changes:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 1: Documentation — NLO Current Skeleton Design

Create:

```text
docs/nlo_current/nlo_current_skeleton_design.md
```

This document must include the following.

---

### 1. Normal form

Use the density normal form:

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W).
\]

The corresponding current is:

\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

The corresponding velocity is:

\[
v^A
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
K_3^{ABC}(H_{BC}+s_Bs_C)
\right].
\]

where

\[
s_A=L_A\log W,
\qquad
H_{BC}=L_Bs_C.
\]

---

### 2. Sector assignment

State the current sector map:

#### \(K_{JSJ}\)

Contributes to \(K_2\) only.

\[
K_{JSJ}: K_2.
\]

LO-like score current.

#### \(K_{JSSJ}\)

Contributes to generic second-order \(K_2\)-type current.

\[
K_{JSSJ}: K_2.
\]

Requires score and coefficient derivatives, but not Hessian-score.

#### \(K_{q\bar q}\)

Contributes to generic second-order \(K_2\)-type current.

\[
K_{q\bar q}: K_2.
\]

Requires score and coefficient derivatives, but not Hessian-score.

#### \(K_{JJSJ}\)

Contributes to:

\[
K_3^{JJSJ}
+
K_{2,\rm comm}^{JJSJ}
+
K_{1,\rm comm}^{JJSJ}.
\]

The previous diagnostic found:

\[
K_{2,\rm comm}^{JJSJ}\neq0,
\]

and

\[
K_{1,\rm comm}^{JJSJ}=0
\]

in that synthetic diagnostic, but the interface must allow a nonzero \(K_1\).

#### \(K_{JJSSJ}\)

Contributes to:

\[
K_3^{JJSSJ}
+
K_{2,\rm comm}^{JJSSJ}
+
K_{1,\rm comm}^{JJSSJ}.
\]

The previous diagnostic found both nonzero:

\[
K_{2,\rm comm}^{JJSSJ}\neq0,
\qquad
K_{1,\rm comm}^{JJSSJ}\neq0.
\]

---

### 3. Data model

Define a small dense data structure for non-production diagnostics:

```python
@dataclass
class NLOCurrentTerms:
    K1: np.ndarray
    K2: np.ndarray
    K3: np.ndarray
    metadata: dict
```

Suggested shapes for \(N_{\rm site}=N\), \(N_c^2-1=8\), and combined dimension:

\[
D=8N.
\]

```text
K1: shape (D,)
K2: shape (D,D)
K3: shape (D,D,D)
```

The metadata should record sector contributions, norms, symmetry diagnostics, and whether commutator corrections were included.

---

### 4. Velocity evaluator scope

The skeleton may evaluate:

\[
v^A
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
K_3^{ABC}(H_{BC}+s_Bs_C)
\right].
\]

For this workflow, coefficient derivatives may be handled in one of two ways:

1. finite-difference diagnostic callbacks, or
2. explicitly passed dense arrays:
   \[
   L_BK_2^{AB},
   \quad
   L_BL_CK_3^{ABC},
   \quad
   L_CK_3^{ABC},
   \quad
   L_BK_3^{ABC}.
   \]

Do not attempt production automatic differentiation.

---

### 5. Non-production warning

State clearly:

\[
\boxed{
\text{This skeleton is not a production NLO evolution implementation.}
}
\]

It is only for validating interfaces, signs, term assembly, shapes, and diagnostic velocity evaluation.

---

## Phase 2: Implement Dense Current-Term Data Structures

Create:

```text
src/nlo_current/nlo_current_skeleton.py
```

Implement:

```python
from dataclasses import dataclass, field
import numpy as np

@dataclass
class NLOCurrentTerms:
    K1: np.ndarray
    K2: np.ndarray
    K3: np.ndarray
    metadata: dict = field(default_factory=dict)

    @property
    def dim(self) -> int:
        ...

    def validate_shapes(self) -> None:
        ...

    def norms(self) -> dict:
        ...

def combined_dim(nsite: int, n_color: int = 8) -> int:
    ...

def flatten_index(site: int, color: int, n_color: int = 8) -> int:
    ...

def unflatten_index(index: int, n_color: int = 8) -> tuple[int, int]:
    ...

def empty_terms(nsite: int, n_color: int = 8) -> NLOCurrentTerms:
    ...

def add_terms(lhs: NLOCurrentTerms, rhs: NLOCurrentTerms, label: str | None = None) -> NLOCurrentTerms:
    """
    Return a new NLOCurrentTerms with K1/K2/K3 added.
    Merge metadata safely.
    """

def symmetrize_K2(K2: np.ndarray) -> np.ndarray:
    """
    Optional helper for diagnostics only.
    Do not force generic two-generator terms to be symmetric.
    """

def symmetrize_K3(K3: np.ndarray) -> np.ndarray:
    """
    Optional helper for canonical normal-form diagnostics only.
    Do not silently symmetrize ordered terms before commutator canonicalization.
    """
```

Important:

- `validate_shapes` should catch mismatched dimensions.
- Do not symmetrize by default.
- Metadata should preserve sector labels and warnings.

---

## Phase 3: Sector Assembly Functions

In:

```text
src/nlo_current/nlo_current_skeleton.py
```

add non-production sector assembly functions.

These may use existing utilities from:

```text
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py
src/nlo_current/cubic_commutator_terms.py
```

Required functions:

```python
def assemble_kjsj_terms(U_fund, S_adj, KJSJ, metadata_only=False) -> NLOCurrentTerms:
    """
    Build the K_JSJ contribution to K2.
    If full dense assembly is not available, create a clear TODO path and metadata entry.
    """

def assemble_kjssj_terms(U_fund, S_adj, KJSSJ, f, metadata_only=False) -> NLOCurrentTerms:
    """
    Build the K_JSSJ contribution to K2 using the existing A/C-left utilities.
    """

def assemble_kqbarq_terms(U_fund, S_adj, Kqbarq, gens, metadata_only=False) -> NLOCurrentTerms:
    """
    Build the K_qbarq contribution to K2 using existing utilities.
    """

def assemble_kjjsj_terms(U_fund, S_adj, KJJSJ, f, include_commutators=True, metadata_only=False) -> NLOCurrentTerms:
    """
    Build K_JJSJ contribution:
        K3 canonical part
        K2_comm
        K1_comm
    using existing cubic and commutator utilities.
    """

def assemble_kjjssj_terms(U_fund, S_adj, KJJSSJ, f, include_commutators=True, metadata_only=False) -> NLOCurrentTerms:
    """
    Build K_JJSSJ contribution:
        K3 canonical part
        K2_comm
        K1_comm
    using existing cubic and commutator utilities.
    """

def assemble_nlo_current_terms(
    U_fund,
    S_adj,
    kernels: dict,
    gens,
    f,
    include_commutators: bool = True,
    metadata_only: bool = False,
) -> NLOCurrentTerms:
    """
    Assemble all currently validated NLO current pieces into K1, K2, K3.

    Expected kernel keys:
        "KJSJ"
        "KJSSJ"
        "Kqbarq"
        "KJJSJ"
        "KJJSSJ"

    Missing kernels should be allowed and recorded in metadata.
    """
```

If a complete dense assembly of a sector is too risky, implement `metadata_only=True` cleanly and document what remains. But at minimum, the skeleton must be able to assemble synthetic dense small-lattice contributions for all five sectors.

---

## Phase 4: Velocity Evaluator

Create:

```text
src/nlo_current/nlo_velocity_evaluator.py
```

Implement:

```python
def evaluate_velocity_from_terms(
    terms,
    score,
    hessian_score,
    dK2=None,
    dK3_first=None,
    d2K3=None,
):
    """
    Evaluate the diagnostic velocity:

    v^A =
        K1^A
        - 1/2 [dK2^A + K2^{AB} s_B]
        + 1/6 [d2K3^A
               + dK3_first_from_C^AB? * s_B
               + dK3_first_from_B^AC? * s_C
               + K3^{ABC}(H_BC + s_B s_C)]

    For this diagnostic skeleton, use explicit derivative arrays.

    Suggested shapes:
        score: (D,)
        hessian_score: (D,D)
        dK2: (D,), representing L_B K2^{AB}
        dK3_first: tuple or dict containing two arrays:
            "LC_K3_ABC": shape (D,D), representing (L_C K3^{ABC})
            "LB_K3_ABC": shape (D,D), representing (L_B K3^{ABC})
        d2K3: (D,), representing L_B L_C K3^{ABC}

    If derivative arrays are None, treat them as zero but add a warning to metadata or return diagnostics.
    """

def evaluate_velocity_score_only(terms, score):
    """
    Convenience diagnostic:
        v = K1 - 1/2 K2 @ score + 1/6 K3 contracted with score*score
    ignoring coefficient derivatives and Hessian-score.

    This is not the full NLO velocity.
    It is only for shape/smoke tests.
    """

def cubic_density_contraction(K3, score, hessian_score):
    """
    Return vector:
        K3^{ABC}(H_BC + s_B s_C).
    """
```

Use `numpy.einsum` carefully and add shape validation.

Important:

- If coefficient derivative arrays are omitted, issue a diagnostic warning in the returned metadata or via a companion function.
- Do not silently present derivative-free velocity as complete.

---

## Phase 5: Synthetic Kernel Factory

Create:

```text
src/nlo_current/synthetic_kernels.py
```

Implement small dense synthetic kernels for tests:

```python
def synthetic_kernels_all(nsite, rng):
    """
    Return a dict with synthetic kernels:
        KJSJ: shape (nsite,nsite,nsite)
        KJSSJ: shape (nsite,nsite,nsite,nsite)
        Kqbarq: shape (nsite,nsite,nsite,nsite)
        KJJSJ: shape (nsite,nsite,nsite,nsite)
        KJJSSJ: shape (nsite,nsite,nsite,nsite,nsite)

    Impose the same diagnostic symmetries used in previous workflows where appropriate:
        KJSSJ: x/y and optional z/zp symmetry
        Kqbarq: x/y and optional z/zp symmetry
        KJJSJ: tested antisymmetric convention in x/y
        KJJSSJ: K(w;x,y;z,z') = -K(w;y,x;z',z)
    """
```

Also include individual helpers if needed.

Do not attempt physical kernels.

---

## Phase 6: Tests for Data Structures and Velocity Evaluator

Create:

```text
tests/nlo_current/test_nlo_current_skeleton.py
```

Tests:

### Test 1: empty terms shape

For \(N_{\rm site}=3\), verify:

\[
D=24
\]

and shapes:

```text
K1: (24,)
K2: (24,24)
K3: (24,24,24)
```

### Test 2: add_terms

Create two random `NLOCurrentTerms` objects and verify addition and metadata merge.

### Test 3: velocity evaluator shape

Use random small `K1,K2,K3`, score, and Hessian-score. Verify output shape:

```text
(D,)
```

### Test 4: cubic contraction identity

Compare `cubic_density_contraction` against explicit loops for small \(D=3\).

### Test 5: missing derivative warning

Call `evaluate_velocity_from_terms` without coefficient derivative arrays. Verify that the result is produced and that diagnostics/warnings indicate derivative terms were treated as zero.

---

## Phase 7: Tests for Sector Assembly

Create:

```text
tests/nlo_current/test_nlo_sector_assembly.py
```

Use:

```text
src/nlo_current/synthetic_kernels.py
```

and random SU(3) Wilson lines from existing utilities.

Tests:

### Test 1: assemble all sectors metadata

Call:

```python
assemble_nlo_current_terms(...)
```

with all five synthetic kernels.

Verify metadata contains entries for:

```text
KJSJ
KJSSJ
Kqbarq
KJJSJ
KJJSSJ
commutators
```

### Test 2: output shapes

Verify assembled `K1,K2,K3` shapes are correct.

### Test 3: commutator toggle

Call with:

```python
include_commutators=True
```

and:

```python
include_commutators=False
```

Verify that the metadata records the toggle.

If possible, verify that `K2` and/or `K1` differ when commutators are included.

### Test 4: missing kernel handling

Pass only a subset of kernels. Verify no crash, zero terms for missing sectors, and metadata notes missing kernels.

### Test 5: no silent symmetrization

For generic ordered terms, verify the assembled `K2` is not forcibly symmetrized unless explicitly requested.

---

## Phase 8: Diagnostic Script

Create:

```text
scripts/nlo_current/build_nlo_current_skeleton_demo.py
```

It should:

1. Generate random SU(3) Wilson lines for \(N_{\rm site}=2\) or \(3\).
2. Generate all synthetic kernels.
3. Assemble:

```python
terms = assemble_nlo_current_terms(...)
```

with commutators enabled.
4. Print/write:
   - norms of \(K_1,K_2,K_3\);
   - sector contribution norms;
   - whether commutators were included;
   - whether coefficient derivative arrays were omitted in velocity evaluation.
5. Create a random score and Hessian-score.
6. Evaluate diagnostic velocity using:
   - full interface with derivative arrays omitted but warning recorded;
   - score-only convenience evaluator.
7. Save:

```text
reports/nlo_current/nlo_current_skeleton_demo_report.md
```

The report must include:

- random seed;
- \(N_{\rm site}\);
- combined dimension \(D\);
- sector norm table;
- total \(K_1,K_2,K_3\) norms;
- velocity norm;
- explicit warnings:
  - physical kernels not used;
  - coefficient derivatives omitted or synthetic;
  - non-production only.

---

## Phase 9: Documentation Update

Create/update:

```text
docs/nlo_current/nlo_current_map_status.md
```

This should be a high-level status document suitable for later paper notes.

It must include:

### Current map

\[
\partial_YW=-L_AJ^A_{\rm NLO},
\]

\[
J^A_{\rm NLO}
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

### Dependency table

| sector | normal-form contribution | density derivatives needed |
|---|---|---|
| \(K_{JSJ}\) | \(K_2\) | score |
| \(K_{JSSJ}\) | \(K_2\) | score + coefficient drift |
| \(K_{q\bar q}\) | \(K_2\) | score + coefficient drift |
| \(K_{JJSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |
| \(K_{JJSSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |

### Important caveats

- This is not a production implementation.
- Physical kernels and barred/unbarred choices still need final integration.
- Coefficient derivatives need a robust implementation strategy.
- Score/Hessian-score estimation is not implemented.
- Dipole validation remains skeleton-only for several terms.

---

## Phase 10: Acceptance Criteria

Stop when all are true:

1. Design document exists:
   ```text
   docs/nlo_current/nlo_current_skeleton_design.md
   ```

2. Skeleton module exists:
   ```text
   src/nlo_current/nlo_current_skeleton.py
   ```

3. Velocity evaluator exists:
   ```text
   src/nlo_current/nlo_velocity_evaluator.py
   ```

4. Synthetic kernels module exists:
   ```text
   src/nlo_current/synthetic_kernels.py
   ```

5. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

6. Demo report exists:
   ```text
   reports/nlo_current/nlo_current_skeleton_demo_report.md
   ```

7. Status doc exists:
   ```text
   docs/nlo_current/nlo_current_map_status.md
   ```

8. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

9. No production evolution code modified.
10. No score/Hessian-score model training implemented.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether `NLOCurrentTerms` shape validation passed.
4. Whether all five sectors can be represented in the skeleton metadata:
   \[
   K_{JSJ},K_{JSSJ},K_{q\bar q},K_{JJSJ},K_{JJSSJ}.
   \]
5. Whether commutator corrections are included in the skeleton interface.
6. Whether diagnostic velocity evaluation works with supplied score and Hessian-score.
7. Which coefficient derivative terms are currently explicit inputs or omitted.
8. Any remaining issue before production design.
9. Recommended next step:
   - coefficient-derivative implementation strategy;
   - physical-kernel integration;
   - dipole validation;
   - score/Hessian-score estimator design.

Do not claim the NLO current implementation is production-ready. This workflow assembles a non-production diagnostic skeleton only.
