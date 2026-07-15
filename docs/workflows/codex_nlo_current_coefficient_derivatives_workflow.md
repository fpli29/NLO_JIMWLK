# Codex Workflow: Coefficient-Derivative Diagnostic Backend for the NLO JIMWLK Current Skeleton

## Purpose

Continue from the completed non-production NLO current skeleton workflow.

Current status:

1. The NLO current skeleton exists and represents the density-side normal form

   \[
   \partial_YW
   =
   -L_A(K_1^AW)
   +
   \frac12L_AL_B(K_2^{AB}W)
   -
   \frac16L_AL_BL_C(K_3^{ABC}W).
   \]

2. The current and velocity are represented as

   \[
   J^A
   =
   K_1^AW
   -
   \frac12L_B(K_2^{AB}W)
   +
   \frac16L_BL_C(K_3^{ABC}W),
   \]

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

3. The current skeleton accepts coefficient-derivative inputs:

   ```text
   dK2
   dK3_first["LC_K3_ABC"]
   dK3_first["LB_K3_ABC"]
   d2K3
   ```

   but currently treats them as optional external arrays and uses zero with warnings when omitted.

The goal of this workflow is to add a **dense small-lattice diagnostic backend** for coefficient derivatives. This backend should compute the derivative arrays by finite differences from coefficient callbacks:

\[
K_2(U),\qquad K_3(U).
\]

Do **not** implement production automatic differentiation.  
Do **not** implement production NLO evolution.  
Do **not** train score/Hessian-score models.  
Do **not** optimize for large lattices.  
This is diagnostic finite-difference infrastructure only.

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
git checkout -b nlo-current-coefficient-derivatives
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect:

```text
docs/nlo_current/nlo_current_skeleton_design.md
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/cubic_current_with_commutator_corrections.md

src/nlo_current/nlo_current_skeleton.py
src/nlo_current/nlo_velocity_evaluator.py
src/nlo_current/synthetic_kernels.py
src/nlo_current/finite_difference_scores.py
src/nlo_current/su3_adjoint.py

tests/nlo_current/test_nlo_current_skeleton.py
tests/nlo_current/test_nlo_sector_assembly.py

scripts/nlo_current/build_nlo_current_skeleton_demo.py
reports/nlo_current/nlo_current_skeleton_demo_report.md
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/coefficient_derivatives_start_status.md
```

It should state:

- whether previous skeleton files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- how `evaluate_velocity_from_terms(...)` currently expects coefficient-derivative arrays;
- that this workflow is diagnostic finite-difference only.

Run before changes:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 1: Documentation — Coefficient Derivative Strategy

Create:

```text
docs/nlo_current/coefficient_derivative_strategy.md
```

This document must include the following.

### 1. Why coefficient derivatives are needed

From the NLO velocity formula,

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

the coefficient-derivative arrays are:

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

These are derivatives of known Wilson-line coefficient functions, not learned density derivatives.

### 2. Diagnostic finite-difference backend

For small dense tests, allow coefficient callbacks:

```python
def K2_callback(U_fund, S_adj) -> np.ndarray:
    # returns K2 with shape (D,D)

def K3_callback(U_fund, S_adj) -> np.ndarray:
    # returns K3 with shape (D,D,D)
```

Then compute finite-difference derivatives using left perturbations.

For combined index

\[
B=(x,b),
\]

use central differences:

\[
L_BF(U)
\approx
\frac{F(e^{+\epsilon t^b}U_x)-F(e^{-\epsilon t^b}U_x)}{2\epsilon}.
\]

The backend should produce:

```text
dK2: shape (D,)
LC_K3_ABC: shape (D,D)
LB_K3_ABC: shape (D,D)
d2K3: shape (D,)
```

### 3. Definitions and contractions

For

\[
K_2^{AB},
\]

\[
dK2^A=\sum_B L_BK_2^{AB}.
\]

For

\[
K_3^{ABC},
\]

\[
(LC\_K3)^{AB}=\sum_C L_CK_3^{ABC},
\]

\[
(LB\_K3)^{AC}=\sum_B L_BK_3^{ABC},
\]

\[
d2K3^A=\sum_{B,C}L_BL_CK_3^{ABC}.
\]

Note that \(L_BL_C\) is ordered. For same-site \(B,C\), commutator issues belong to the canonicalized \(K_3,K_2,K_1\) representation, but the finite-difference backend should apply derivatives in the order specified by the velocity formula.

### 4. Non-production warning

State clearly:

\[
\boxed{
\text{This finite-difference backend is not a production coefficient-derivative implementation.}
}
\]

It is only for dense diagnostics, product-rule tests, and velocity sensitivity studies.

---

## Phase 2: Implement Coefficient Derivative Backend

Create:

```text
src/nlo_current/coefficient_derivatives.py
```

Implement:

```python
import numpy as np

def validate_coefficient_shapes(K2, K3):
    """
    Validate K2 shape (D,D) and K3 shape (D,D,D), with matching D.
    Return D.
    """

def left_perturbed_copy(U_fund, site, color, gens, eps):
    """
    Return a copy of U_fund with U_site left-perturbed by exp(eps * t^color)
    using the same convention as existing finite-difference utilities.

    Reuse existing left_perturb/random SU(3) utilities if available.
    """

def fd_left_derivative_array(callback, U_fund, S_adj_builder, site, color, gens, eps):
    """
    Central finite difference derivative of an array-valued callback.

    callback(U_fund, S_adj) -> np.ndarray
    S_adj_builder(U_fund) -> S_adj

    Return:
        L_{site,color} callback(U)
    """

def compute_dK2_fd(K2_callback, U_fund, S_adj_builder, gens, eps=1e-5):
    """
    Compute dK2^A = sum_B L_B K2^{A B}.

    Return:
        dK2 shape (D,)
    """

def compute_dK3_first_fd(K3_callback, U_fund, S_adj_builder, gens, eps=1e-5):
    """
    Compute the two first-derivative contractions:

        LC_K3_ABC[A,B] = sum_C L_C K3[A,B,C]
        LB_K3_ABC[A,C] = sum_B L_B K3[A,B,C]

    Return:
        {
          "LC_K3_ABC": array shape (D,D),
          "LB_K3_ABC": array shape (D,D),
        }
    """

def compute_d2K3_fd(K3_callback, U_fund, S_adj_builder, gens, eps=1e-4):
    """
    Compute d2K3^A = sum_{B,C} L_B L_C K3^{A B C}.

    Use nested central finite differences.

    This is expensive and only intended for D <= O(24).
    Return:
        d2K3 shape (D,)
    """

def compute_all_coefficient_derivatives_fd(
    K2_callback,
    K3_callback,
    U_fund,
    S_adj_builder,
    gens,
    eps_first=1e-5,
    eps_second=1e-4,
):
    """
    Convenience wrapper returning:
        dK2, dK3_first, d2K3
    """
```

Important implementation notes:

- Do not assume \(N_{\rm site}=3\); infer from `U_fund`.
- Use shape validation.
- For `S_adj_builder`, reuse existing adjoint conversion utilities.
- Keep explicit warnings/comments that nested finite difference is expensive.
- Keep code clear rather than optimized.
- If existing finite-difference convention uses anti-Hermitian generators or `expm`, follow that convention exactly.

---

## Phase 3: Add Product-Rule Utilities

Create or update:

```text
src/nlo_current/coefficient_derivatives.py
```

Add diagnostic helper functions:

```python
def product_rule_K2_rhs(K2, dK2, score):
    """
    Return vector:
        dK2^A + K2^{AB}s_B
    """

def product_rule_K3_rhs(K3, dK3_first, d2K3, score, hessian_score):
    """
    Return vector:
        d2K3^A
        + (LC_K3_ABC)^{AB}s_B
        + (LB_K3_ABC)^{AC}s_C
        + K3^{ABC}(H_BC+s_Bs_C)
    """

def velocity_from_coeff_derivative_backend(terms, score, hessian_score, derivatives):
    """
    Thin wrapper around evaluate_velocity_from_terms using the computed derivative arrays.
    """
```

These helpers should match the contractions used in:

```text
src/nlo_current/nlo_velocity_evaluator.py
```

---

## Phase 4: Tests for Coefficient Derivative Shapes and Contractions

Create:

```text
tests/nlo_current/test_coefficient_derivatives.py
```

Tests:

### Test 1: zero coefficient derivatives for constant callbacks

Define constant \(K_2,K_3\) callbacks independent of \(U\). Verify:

\[
dK2=0,
\qquad
LC\_K3=0,
\qquad
LB\_K3=0,
\qquad
d2K3=0.
\]

### Test 2: shape validation

For \(N_{\rm site}=2\), \(D=16\), verify:

```text
dK2: (D,)
LC_K3_ABC: (D,D)
LB_K3_ABC: (D,D)
d2K3: (D,)
```

### Test 3: linear callback gives nonzero first derivative

Define a simple coefficient callback depending on one Wilson line, e.g.

\[
K_2^{AB}(U)=M^{AB}\,\mathrm{ReTr}(Q U_0).
\]

Verify `compute_dK2_fd` is nonzero and stable over eps.

### Test 4: product-rule K2 contraction

For synthetic \(K_2(U)\), toy log density \(\log W\), and score \(s_B\), verify finite-difference approximation of

\[
\frac{1}{W}L_B(K_2^{AB}W)
\]

matches

\[
L_BK_2^{AB}+K_2^{AB}s_B.
\]

This can be tested for a fixed \(A\) or all \(A\) on a small \(D\).

### Test 5: product-rule K3 contraction

For synthetic \(K_3(U)\), toy log density, score, and Hessian-score, verify finite-difference approximation of

\[
\frac{1}{W}L_BL_C(K_3^{ABC}W)
\]

matches

\[
L_BL_CK_3^{ABC}
+
(L_CK_3^{ABC})s_B
+
(L_BK_3^{ABC})s_C
+
K_3^{ABC}(H_{BC}+s_Bs_C).
\]

Use small \(N_{\rm site}=1\) or \(2\) and loose finite-difference tolerances.  
If full all-index contraction is too expensive, test selected indices and document it.

Run:

```bash
python3 -m pytest tests/nlo_current/test_coefficient_derivatives.py -q
```

---

## Phase 5: Integrate with Velocity Evaluator Tests

Update:

```text
tests/nlo_current/test_nlo_current_skeleton.py
```

or create:

```text
tests/nlo_current/test_nlo_velocity_with_coefficient_derivatives.py
```

Tests:

### Test 1: derivative arrays remove omission warnings

Call `evaluate_velocity_from_terms(...)` with computed derivative arrays and verify that the diagnostics no longer say coefficient derivatives were treated as zero.

### Test 2: omitting derivatives changes velocity

For a nonconstant synthetic coefficient callback, compare:

```python
v_without_derivatives
```

against:

```python
v_with_derivatives
```

Verify:

\[
\|v_{\rm with}-v_{\rm without}\|>0
\]

above a small threshold.

### Test 3: constant coefficients agree with derivative-free path

For constant \(K_2,K_3\), computed derivative arrays are zero, so velocity with derivatives should match velocity with omitted derivative arrays, up to warnings.

---

## Phase 6: Diagnostic Script

Create:

```text
scripts/nlo_current/check_coefficient_derivative_backend.py
```

It should:

1. Generate random SU(3) Wilson lines for \(N_{\rm site}=1\) or \(2\).
2. Build synthetic nonconstant \(K_2(U),K_3(U)\) callbacks.
3. Compute:

   ```text
   dK2
   LC_K3_ABC
   LB_K3_ABC
   d2K3
   ```

4. Assemble a small `NLOCurrentTerms`.
5. Generate toy score and Hessian-score.
6. Evaluate velocity with and without coefficient derivatives.
7. Save:

```text
reports/nlo_current/coefficient_derivative_backend_report.md
```

The report must include:

- random seed;
- \(N_{\rm site}\);
- combined dimension \(D\);
- eps values;
- norms:
  - \(\|dK2\|\)
  - \(\|LC\_K3\|\)
  - \(\|LB\_K3\|\)
  - \(\|d2K3\|\)
  - \(\|v_{\rm with}\|\)
  - \(\|v_{\rm without}\|\)
  - \(\|v_{\rm with}-v_{\rm without}\|\)
- product-rule residuals if available;
- warning that this backend is finite-difference diagnostic only.

---

## Phase 7: Update Skeleton Demo

Update:

```text
scripts/nlo_current/build_nlo_current_skeleton_demo.py
```

Add an optional path:

```python
--with-coefficient-derivatives
```

or just add a second demo section that computes finite-difference coefficient derivatives for very small \(N_{\rm site}=1\) or \(2\).

Update:

```text
reports/nlo_current/nlo_current_skeleton_demo_report.md
```

to include:

- derivative-enabled velocity norm;
- derivative-omitted velocity norm;
- their difference;
- note that finite-difference derivatives are only feasible for tiny diagnostics.

Do not make the default demo too slow.

---

## Phase 8: Documentation Update

Create/update:

```text
docs/nlo_current/nlo_current_map_status.md
```

Add a section:

```text
Coefficient derivative status
```

Include:

- coefficient derivative backend exists;
- finite-difference diagnostic only;
- velocity evaluator can now be run with explicit derivative arrays;
- omitted derivative path remains available only for smoke tests and warns;
- production strategy still unresolved.

Also create:

```text
docs/nlo_current/coefficient_derivative_backend_limitations.md
```

Include limitations:

- \(O(D)\) first derivatives of full arrays;
- \(O(D^2)\) nested second derivatives for \(K_3\);
- finite-difference roundoff at small eps;
- not suitable for realistic lattice sizes;
- future production options:
  - analytic coefficient derivatives;
  - automatic differentiation;
  - sparse/local kernel structure;
  - stochastic trace/derivative estimators.

---

## Phase 9: Acceptance Criteria

Stop when all are true:

1. Strategy doc exists:
   ```text
   docs/nlo_current/coefficient_derivative_strategy.md
   ```

2. Backend exists:
   ```text
   src/nlo_current/coefficient_derivatives.py
   ```

3. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

4. Diagnostic report exists:
   ```text
   reports/nlo_current/coefficient_derivative_backend_report.md
   ```

5. Skeleton demo report updated:
   ```text
   reports/nlo_current/nlo_current_skeleton_demo_report.md
   ```

6. Limitations doc exists:
   ```text
   docs/nlo_current/coefficient_derivative_backend_limitations.md
   ```

7. Map status doc updated:
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
3. Whether constant coefficient callbacks give zero derivatives.
4. Whether nonconstant coefficient callbacks give nonzero derivatives.
5. Whether \(K_2\) product-rule tests passed.
6. Whether \(K_3\) product-rule tests passed or were limited to selected indices.
7. Norms of:
   \[
   dK2,\quad LC\_K3,\quad LB\_K3,\quad d2K3.
   \]
8. Whether velocity with coefficient derivatives differs from velocity with omitted derivatives.
9. Any finite-difference stability or eps sensitivity issues.
10. Recommended next step:
    - physical-kernel integration,
    - full dipole validation,
    - analytic/AD coefficient derivative design,
    - score/Hessian-score estimator design.

Do not claim production readiness. This workflow only adds a finite-difference diagnostic coefficient-derivative backend.
