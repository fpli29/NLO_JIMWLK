# Codex Workflow: Derive and Test the \(K_{JJSJ}\) Cubic Ordered Current

## Purpose

Continue from the completed two-generator NLO JIMWLK current workflows.

The two-generator sector is now stable:

\[
K_{JSJ}: \text{LO-like symmetric score current},
\]

\[
K_{JSSJ},\ K_{q\bar q}: \text{generic ordered }J_LAJ_R\text{ currents with coefficient drift}.
\]

The current goal is to start the three-generator sector with the simpler cubic term:

\[
K_{JJSJ}.
\]

Do **not** implement the full NLO flow.  
Do **not** start \(K_{JJSSJ}\) yet.  
Only complete the \(K_{JJSJ}\)-specific derivation, coefficient builders, finite-difference checks, and reports.

The key scientific objective is to confirm that \(K_{JJSJ}\) requires both:

\[
s_A=L_A\log W
\]

and

\[
H_{AB}=L_As_B=L_AL_B\log W,
\]

i.e. score and Hessian-score.

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
3. Before editing an existing file, report the file and the reason.
4. Maintain/update:
   ```text
   reports/nlo_current/file_manifest.md
   ```
   listing every created or modified file.

If this is a Git repository, create a branch:

```bash
git checkout -b nlo-current-kjjsj-cubic
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect the artifacts from the two-generator workflows:

```text
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
docs/nlo_current/Kqbarq_ordered_current.md
docs/nlo_current/two_generator_sector_summary.md
src/nlo_current/su3_adjoint.py
src/nlo_current/two_generator_terms.py
tests/nlo_current/test_su3_adjoint_conventions.py
tests/nlo_current/test_kqbarq_coefficient.py
scripts/nlo_current/check_kjssj_symmetry.py
scripts/nlo_current/check_kqbarq_symmetry.py
reports/nlo_current/kjssj_symmetry_report.md
reports/nlo_current/kqbarq_symmetry_report.md
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/kjjsj_start_status.md
```

It should state:

- whether previous workflow files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- the currently implemented conventions:
  \[
  S_A^{ab}=2\mathrm{ReTr}(t^a U t^b U^\dagger),
  \]
  \[
  J_R^a=S_A^{ba}J_L^b,
  \]
  and left perturbation convention.

Run before starting changes:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 1: Add the \(K_{JJSJ}\) Derivation Document

Create:

```text
docs/nlo_current/KJJSJ_cubic_ordered_current.md
```

This document must contain the following sections.

---

### 1. Hamiltonian structure

Use the five-kernel KLM form. The \(K_{JJSJ}\) term is:

\[
H_{JJSJ}
=
\int_{w,x,y,z}
K_{JJSJ}(w;x,y;z) f^{bde}
\left[
J_L^d(x)J_L^e(y)S_z^{ba}J_R^a(w)
-
J_L^a(w)S_z^{ab}J_R^d(x)J_R^e(y)
+
\frac13
\left(
J_L^d(x)J_L^e(y)J_L^b(w)
-
J_R^d(x)J_R^e(y)J_R^b(w)
\right)
\right].
\]

This has three pieces:

1. LLR real-like piece:
   \[
   J_LJ_LSJ_R.
   \]
2. LRR real-like piece:
   \[
   -J_LSJ_RJ_R.
   \]
3. \(1/3\) virtual piece:
   \[
   \frac13(J_LJ_LJ_L-J_RJ_RJ_R).
   \]

This term is cubic in charge generators, so it should generate score and Hessian-score contributions.

---

### 2. Ordered cubic current lemma: LLR block

Define a generic ordered LLR block:

\[
H_{LLR}[A]
=
\int A^{dea}(x,y,w;U)
J_L^d(x)J_L^e(y)J_R^a(w).
\]

With observable evolution

\[
\partial_Y{\cal O}=-H_{LLR}{\cal O},
\]

the density-side operator is

\[
(\partial_YW)_{LLR}
=
+J_R^a(w)J_L^e(y)J_L^d(x)
\left[
A^{dea}W
\right],
\]

because three integrations by parts give a minus sign from the adjoint, and the outer KLM minus gives another sign. Check the sign numerically.

Using

\[
J_R^a(w)=S_w^{ha}L_w^h,
\]

and assuming the verified divergence-free relation

\[
L_w^hS_w^{ha}=0,
\]

a natural current component is

\[
\boxed{
J_{LLR}^{(w,h)}
=
- S_w^{ha}
L_y^eL_x^d[A^{dea}W].
}
\]

This is because

\[
(\partial_YW)_{LLR}
=
-L_w^h J_{LLR}^{(w,h)}.
\]

Then the velocity is

\[
\boxed{
v_{LLR}^{(w,h)}
=
- S_w^{ha}
\frac{1}{W}L_y^eL_x^d[A^{dea}W].
}
\]

Expand:

\[
\frac{1}{W}L_y^eL_x^d[A^{dea}W]
=
L_y^eL_x^dA^{dea}
+
(L_x^dA^{dea})s_y^e
+
(L_y^eA^{dea})s_x^d
+
A^{dea}\left(L_y^es_x^d+s_y^es_x^d\right).
\]

Thus

\[
\boxed{
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
A^{dea}\left(H_{yx}^{ed}+s_y^es_x^d\right)
\right],
}
\]

where

\[
H_{yx}^{ed}=L_y^es_x^d.
\]

This is the first place Hessian-score enters.

**Important:** Codex must check this sign and derivative ordering with finite differences. If a sign differs, update this document and the report.

---

### 3. Ordered cubic current lemma: LRR block

Define

\[
H_{LRR}[B]
=
\int B^{ade}(w,x,y;U)
J_L^a(w)J_R^d(x)J_R^e(y).
\]

Convert right generators:

\[
J_R^d(x)=S_x^{pd}L_x^p,
\qquad
J_R^e(y)=S_y^{qe}L_y^q.
\]

A safe approach is to express the density-side current directly from the adjoint normal ordering and verify by finite differences.

The expected structure is a current with component at the left-generator coordinate \(w\):

\[
J_{LRR}^{(w,a)}
\sim
+\ \frac{1}{W}
\text{two derivatives acting on }
\left[
B^{ade}S_x^{pd}S_y^{qe}W
\right].
\]

Do not assume the final sign. Derive and test.

The final expanded form should include:

- coefficient second derivatives;
- coefficient first derivative times score;
- Hessian-score:
  \[
  L_x^p s_y^q
  \quad\text{or}\quad
  L_y^q s_x^p;
  \]
- score product:
  \[
  s_x^p s_y^q.
  \]

Document the final tested formula.

---

### 4. Ordered cubic current lemma: LLL and RRR virtual blocks

The virtual block is:

\[
H_{\rm virt}
=
\frac13
\int V^{deb}(x,y,w;U)
\left[
J_L^d(x)J_L^e(y)J_L^b(w)
-
J_R^d(x)J_R^e(y)J_R^b(w)
\right].
\]

For \(K_{JJSJ}\),

\[
V^{deb}(x,y,w;U)
=
K_{JJSJ}(w;x,y;z)f^{bde}
\]

integrated over \(z\).

Build current lemmas for:

\[
H_{LLL}[V]=V^{ABC}L_AL_BL_C,
\]

and for the RRR term after converting all right generators to left generators.

The generic cubic density normal form is:

\[
\partial_YW\supset -\frac16L_AL_BL_C(K_3^{ABC}W),
\]

with current

\[
J^A_{\rm cubic}
=
\frac16L_BL_C(K_3^{ABC}W).
\]

Then

\[
v^A_{\rm cubic}
=
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

Use this formula only after mapping the ordered KLM coefficient and sign to a canonical \(K_3^{ABC}\). Do not guess this mapping; verify with finite differences.

---

### 5. Apply to \(K_{JJSJ}\)

Define the three observable-side coefficient blocks:

#### LLR block

\[
A_{LLR}^{dea}(x,y,w;U)
=
\int_z
K_{JJSJ}(w;x,y;z)
f^{bde}S_z^{ba}.
\]

#### LRR block

The Hamiltonian contains

\[
- K_{JJSJ}(w;x,y;z)f^{bde}
J_L^a(w)S_z^{ab}J_R^d(x)J_R^e(y).
\]

Define

\[
B_{LRR}^{ade}(w,x,y;U)
=
-\int_z
K_{JJSJ}(w;x,y;z)
f^{bde}S_z^{ab}.
\]

#### Virtual LLL/RRR block

\[
V^{deb}(x,y,w;U)
=
\frac13\int_z
K_{JJSJ}(w;x,y;z)f^{bde}.
\]

The virtual contribution is

\[
V^{deb}
\left[
J_L^d(x)J_L^e(y)J_L^b(w)
-
J_R^d(x)J_R^e(y)J_R^b(w)
\right].
\]

Codex must keep the \(1/3\) factor.

---

### 6. State expected ML object requirements

The \(K_{JJSJ}\) term requires:

\[
s_A=L_A\log W,
\]

and

\[
H_{AB}=L_As_B=L_AL_B\log W.
\]

It may also require coefficient derivatives:

\[
L_AA,\qquad L_AL_BA.
\]

These are derivatives of known Wilson-line coefficient functions, not learned density derivatives.

State clearly:

\[
\boxed{
K_{JJSJ}\text{ is the first term that requires Hessian-score.}
}
\]

---

## Phase 2: Add Cubic Utility Module

Create:

```text
src/nlo_current/three_generator_terms.py
```

Add small dense validation utilities only. Do not implement production-volume kernels.

Required functions:

```python
def kjjsj_A_LLR_from_kernel(S_adj, KJJSJ, f):
    """
    Build A_LLR[x,y,w,d,e,a] =
        sum_z KJJSJ[w,x,y,z] f[b,d,e] S_adj[z,b,a].
    """

def kjjsj_B_LRR_from_kernel(S_adj, KJJSJ, f):
    """
    Build B_LRR[w,x,y,a,d,e] =
        - sum_z KJJSJ[w,x,y,z] f[b,d,e] S_adj[z,a,b].
    """

def kjjsj_V_virtual_from_kernel(KJJSJ, f):
    """
    Build V[x,y,w,d,e,b] =
        (1/3) sum_z KJJSJ[w,x,y,z] f[b,d,e].
    """

def flatten_cubic_index(site, color, n_color=8):
    """Return combined index A=(site,color)."""

def dense_cubic_tensor_from_blocks(...):
    """
    Optional helper: convert LLR/LRR/LLL/RRR blocks to dense combined-index tensors
    for small-lattice diagnostics only.
    """
```

If useful, add:

```python
def synthetic_kjjsj_kernel(nsite, rng, antisym_xy=True):
    """
    Generate a synthetic KJJSJ[w,x,y,z] with the expected antisymmetry in x,y:
        K(w;x,y;z) = -K(w;y,x;z),
    because f^{bde} is antisymmetric in d,e and the Hamiltonian structure should
    preserve ordering/antisymmetry.
    """
```

If the exact symmetry of \(K_{JJSJ}\) in \(x,y\) is already known from the formula, implement it. If not, produce both unconstrained and antisymmetrized synthetic kernels and report both.

---

## Phase 3: Add Finite-Difference Derivative Tools for Score and Hessian

Update or create:

```text
src/nlo_current/finite_difference_scores.py
```

Required functions:

```python
def toy_log_density(U_list, params):
    """
    Return log W for a toy density.
    Example:
        log W = lambda * sum_i ReTr(Q_i U_i)
              + eta * sum_{i<j} ReTr(U_i U_j^\dagger)
    Use a coupled term so Hessian-score is nonzero.
    """

def fd_score(logW_func, U_list, site, color, gens, eps):
    """Finite-difference L_A log W."""

def fd_hessian_score(logW_func, U_list, site_a, color_a, site_b, color_b, gens, eps):
    """Finite-difference L_A s_B = L_A L_B log W."""

def fd_left_derivative_scalar(func, U_list, site, color, gens, eps):
    """General scalar left derivative."""

def fd_left_second_derivative_scalar(func, U_list, site_a, color_a, site_b, color_b, gens, eps):
    """General scalar second left derivative."""
```

Use central differences where possible.

The toy log density must include a coupled term so that at least some off-diagonal Hessian-score components are nonzero.

---

## Phase 4: Tests for Cubic Lemmas

Create:

```text
tests/nlo_current/test_cubic_ordered_current.py
```

Include the following tests.

### Test 1: Hessian-score is nonzero for toy density

Verify that the toy log density produces nonzero Hessian-score for at least one pair of distinct combined indices.

This prevents false confidence from testing a linear density.

### Test 2: LLR cubic expansion identity

For random SU(3) Wilson lines and a small random coefficient \(A^{dea}(U)\), verify numerically:

\[
\frac1W L_y^eL_x^d[A^{dea}W]
\]

matches the expanded expression:

\[
L_y^eL_x^dA^{dea}
+
(L_x^dA^{dea})s_y^e
+
(L_y^eA^{dea})s_x^d
+
A^{dea}(H_{yx}^{ed}+s_y^es_x^d).
\]

Use finite differences.

Acceptance:

- residual decreases as \(\epsilon\) decreases over a moderate range;
- use loose tolerances initially, e.g. \(10^{-4}\) to \(10^{-5}\);
- report residuals for \(\epsilon=10^{-3},10^{-4},10^{-5}\) if stable.

### Test 3: LLR current sign check

Check the sign of the proposed LLR current by comparing:

\[
(\partial_YW)_{LLR}
=
+J_RJ_LJ_L(AW)
\]

against

\[
-L_wJ^{(w)}_{LLR}.
\]

If the sign in the document is wrong, update the document and implementation. Do not force the test.

### Test 4: LRR current lemma smoke test

Derive and test the LRR current sign/ordering using a simple random \(B^{ade}(U)\).

This may be implemented as a finite-difference comparison between:

\[
(\partial_YW)_{LRR}
=
\text{direct adjoint density-side operator}
\]

and the proposed left-divergence current form.

If a clean current form is not obtained in this workflow, write a report and stop before claiming completion.

### Test 5: virtual LLL/RRR mapping smoke test

Build small random \(V\) and verify the canonical cubic normal-form expansion for LLL.

For RRR, convert right generators to left generators and check at least shape/sign consistency. If full finite-difference validation is too much, document the remaining TODO explicitly.

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 5: \(K_{JJSJ}\) Coefficient Builder Tests

Create:

```text
tests/nlo_current/test_kjjsj_coefficients.py
```

Tests:

### Test 1: shapes

For \(N_{\rm site}=3\), verify:

\[
A_{LLR}[x,y,w,d,e,a]
\]

has shape:

```text
(nsite,nsite,nsite,8,8,8)
```

and similarly for \(B_{LRR}\), \(V\).

### Test 2: realness

For real synthetic \(K_{JJSJ}\), real adjoint Wilson lines, and real \(f^{abc}\), verify the coefficient blocks are real.

### Test 3: antisymmetry in \(d,e\)

Because \(f^{bde}\) is antisymmetric,

\[
A_{LLR}^{dea}=-A_{LLR}^{eda}
\]

if the kernel is symmetric under \(x\leftrightarrow y\), or check the correct combined antisymmetry if the kernel has an explicit \(x,y\) antisymmetry.

Do not hard-code an assumption. Test and report both:

- synthetic symmetric \(K(w;x,y;z)=K(w;y,x;z)\);
- synthetic antisymmetric \(K(w;x,y;z)=-K(w;y,x;z)\).

Document which convention makes the full ordered coefficient consistent.

### Test 4: virtual \(1/3\) factor

Verify numerically that removing the \(1/3\) factor changes the virtual coefficient norm by exactly factor 3.

This prevents losing the virtual factor.

---

## Phase 6: Diagnostic Script and Report

Create:

```text
scripts/nlo_current/check_kjjsj_cubic_requirements.py
```

It should:

1. Generate random SU(3) Wilson lines on \(N_{\rm site}=3\).
2. Generate synthetic \(K_{JJSJ}\) kernels under both candidate \(x,y\) symmetry conventions:
   - symmetric;
   - antisymmetric.
3. Build LLR, LRR, and virtual coefficient blocks.
4. Evaluate the norms of:
   - coefficient second-derivative term;
   - score-linear terms;
   - Hessian-score term;
   - score-product term;
   on a toy log density.
5. Save:

```text
reports/nlo_current/kjjsj_cubic_requirements_report.md
```

The report must include:

- random seed;
- kernel symmetry convention tested;
- evidence that Hessian-score terms are nonzero;
- which coefficient blocks were built;
- whether LLR/LRR/virtual current signs passed tests;
- remaining uncertainties.

The report should conclude either:

\[
\text{K}_{JJSJ}\text{ current can be represented as score + Hessian-score with tested signs}
\]

or

\[
\text{K}_{JJSJ}\text{ current still has unresolved ordering/sign issue}
\]

Do not overclaim.

---

## Phase 7: Update Sector Summary

Create or update:

```text
docs/nlo_current/three_generator_sector_start.md
```

It should state:

1. Two-generator sector status:
   - \(K_{JSJ}\): LO-like score current.
   - \(K_{JSSJ}\), \(K_{q\bar q}\): generic ordered-current with score and coefficient drift.
2. First three-generator term:
   - \(K_{JJSJ}\) introduces Hessian-score.
3. Working formula for LLR block after tests.
4. Working formula for LRR block after tests.
5. Treatment of \(1/3\) virtual term.
6. What remains for \(K_{JJSSJ}\).

If signs remain unresolved, make this explicit in the document.

---

## Phase 8: Dipole Validation Skeleton Update

Update:

```text
scripts/nlo_current/validate_dipole_two_generator_terms.py
```

or create a new file:

```text
scripts/nlo_current/validate_dipole_kjjsj_skeleton.py
```

Add the KLM Appendix A target for \(K_{JJSJ}\):

The dipole action should involve combinations like

\[
iN_c^2
\left[
K_{JJSJ}(v;u,v;z)
+
K_{JJSJ}(u;v,u;z)
\right]
\left[
s(u,z)s(z,v)-s(u,v)
\right],
\]

and the more detailed expression from the appendix should be listed as a TODO target.

This is a skeleton only. Do not block workflow completion on full dipole validation.

---

## Phase 9: Acceptance Criteria

Stop when all are true:

1. Derivation doc exists:
   ```text
   docs/nlo_current/KJJSJ_cubic_ordered_current.md
   ```

2. Cubic utility module exists:
   ```text
   src/nlo_current/three_generator_terms.py
   ```

3. Finite-difference score/Hessian utility exists:
   ```text
   src/nlo_current/finite_difference_scores.py
   ```

4. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

5. Report exists:
   ```text
   reports/nlo_current/kjjsj_cubic_requirements_report.md
   ```

6. Summary exists:
   ```text
   docs/nlo_current/three_generator_sector_start.md
   ```

7. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

8. No full NLO flow implemented.
9. No \(K_{JJSSJ}\) implementation started.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether the LLR cubic current sign passed.
4. Whether the LRR cubic current sign passed.
5. Whether the \(1/3\) virtual term was included and tested.
6. Whether Hessian-score terms were nonzero in the toy-density diagnostics.
7. Which \(K_{JJSJ}\) kernel symmetry convention was used/tested.
8. Whether any sign/order convention issue remains.
9. Recommended next step:
   - if \(K_{JJSJ}\) is stable, move to \(K_{JJSSJ}\);
   - otherwise, resolve the remaining cubic ordering/sign issue first.

Do not claim that the NLO current is complete. This workflow only starts and validates the first three-generator term \(K_{JJSJ}\).
