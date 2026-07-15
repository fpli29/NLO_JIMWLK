# Codex Workflow: Derive and Test the \(K_{JJSSJ}\) Cubic Ordered Current

## Purpose

Continue from the completed \(K_{JJSJ}\) cubic ordered-current workflow.

Current status:

1. The two-generator sector is stable:
   \[
   K_{JSJ}: \text{LO-like symmetric score current},
   \]
   \[
   K_{JSSJ},\ K_{q\bar q}: \text{generic ordered }J_LAJ_R\text{ currents with coefficient drift}.
   \]

2. The first cubic term \(K_{JJSJ}\) has been validated in the distinct-site small-lattice setting:
   - LLR cubic current sign passed.
   - LRR cubic current sign passed.
   - \(1/3\) virtual factor passed.
   - Virtual LLL/RRR sign smoke tests passed.
   - Hessian-score contributions were nonzero.
   - The antisymmetric synthetic kernel convention gave the expected combined \((x,d)\leftrightarrow(y,e)\) behavior.

3. Coincident-site commutator handling is still outside scope and must **not** be silently ignored in production.

The current goal is to derive and validate the final cubic NLO term:

\[
K_{JJSSJ}.
\]

Do **not** implement the full NLO flow.  
Do **not** handle coincident-site commutators in this workflow.  
Do **not** modify production evolution code.  
Only complete the distinct-site \(K_{JJSSJ}\) cubic ordered-current derivation, coefficient builders, finite-difference checks, diagnostics, and reports.

The key scientific objective is to confirm that \(K_{JJSSJ}\), like \(K_{JJSJ}\), requires:

\[
s_A=L_A\log W
\]

and

\[
H_{AB}=L_As_B=L_AL_B\log W.
\]

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
git checkout -b nlo-current-kjjssj-cubic
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect the artifacts from previous workflows:

```text
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
docs/nlo_current/Kqbarq_ordered_current.md
docs/nlo_current/two_generator_sector_summary.md
docs/nlo_current/KJJSJ_cubic_ordered_current.md
docs/nlo_current/three_generator_sector_start.md
src/nlo_current/su3_adjoint.py
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py
src/nlo_current/finite_difference_scores.py
tests/nlo_current/test_su3_adjoint_conventions.py
tests/nlo_current/test_kqbarq_coefficient.py
tests/nlo_current/test_cubic_ordered_current.py
tests/nlo_current/test_kjjsj_coefficients.py
scripts/nlo_current/check_kjssj_symmetry.py
scripts/nlo_current/check_kqbarq_symmetry.py
scripts/nlo_current/check_kjjsj_cubic_requirements.py
reports/nlo_current/kjssj_symmetry_report.md
reports/nlo_current/kqbarq_symmetry_report.md
reports/nlo_current/kjjsj_cubic_requirements_report.md
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/kjjssj_start_status.md
```

It should state:

- whether previous workflow files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- currently implemented conventions:
  \[
  S_A^{ab}=2\mathrm{ReTr}(t^a U t^b U^\dagger),
  \]
  \[
  J_R^a=S_A^{ba}J_L^b,
  \]
  and the left perturbation convention;
- the \(K_{JJSJ}\) distinct-site caveat:
  coincident-site commutators remain unresolved.

Run before starting changes:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 1: Add the \(K_{JJSSJ}\) Derivation Document

Create:

```text
docs/nlo_current/KJJSSJ_cubic_ordered_current.md
```

This document must contain the following sections.

---

### 1. Hamiltonian structure

Use the five-kernel KLM form. The \(K_{JJSSJ}\) term is:

\[
H_{JJSSJ}
=
\int_{w,x,y,z,z'}
K_{JJSSJ}(w;x,y;z,z')f^{acb}
\left[
J_L^d(x)J_L^e(y)S_z^{dc}S_{z'}^{eb}J_R^a(w)
-
J_L^a(w)S_z^{cd}S_{z'}^{be}J_R^d(x)J_R^e(y)
+
\frac13
\left(
J_L^c(x)J_L^b(y)J_L^a(w)
-
J_R^c(x)J_R^b(y)J_R^a(w)
\right)
\right].
\]

This term has three pieces:

1. LLR real-like piece:
   \[
   J_LJ_LS(z)S(z')J_R.
   \]
2. LRR real-like piece:
   \[
   -J_LS(z)S(z')J_RJ_R.
   \]
3. \(1/3\) virtual piece:
   \[
   \frac13(J_LJ_LJ_L-J_RJ_RJ_R).
   \]

This term is cubic in charge generators, so it should require score and Hessian-score.

---

### 2. Reuse the tested cubic current lemmas

State that this workflow reuses the distinct-site cubic lemmas validated for \(K_{JJSJ}\):

#### LLR block

For

\[
H_{LLR}[A]
=
A^{DEA}J_L^D J_L^E J_R^A,
\]

with KLM observable convention

\[
\partial_Y{\cal O}=-H{\cal O},
\]

the tested distinct-site current structure is:

\[
v_{LLR}^{(w,h)}
=
-S_w^{ha}
\frac{1}{W}L_y^eL_x^d[A^{dea}W].
\]

Expanded:

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

where

\[
H_{yx}^{ed}=L_y^es_x^d.
\]

#### LRR block

Reuse the LRR current lemma and sign that passed in the \(K_{JJSJ}\) workflow. Document the explicit formula used in code.

#### Virtual LLL/RRR block

Reuse the canonical cubic normal-form expansion and the tested \(1/3\) virtual factor handling from the \(K_{JJSJ}\) workflow.

Do not invent new signs; test all signs again for the \(K_{JJSSJ}\) coefficient blocks.

---

### 3. Define \(K_{JJSSJ}\) coefficient blocks

#### LLR block

The LLR part is

\[
K_{JJSSJ}(w;x,y;z,z')f^{acb}
J_L^d(x)J_L^e(y)S_z^{dc}S_{z'}^{eb}J_R^a(w).
\]

Define

\[
\boxed{
A_{LLR}^{dea}(x,y,w;U)
=
\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}
S_z^{dc}S_{z'}^{eb}.
}
\]

#### LRR block

The LRR part is

\[
-
K_{JJSSJ}(w;x,y;z,z')f^{acb}
J_L^a(w)S_z^{cd}S_{z'}^{be}J_R^d(x)J_R^e(y).
\]

Define

\[
\boxed{
B_{LRR}^{ade}(w,x,y;U)
=
-\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}
S_z^{cd}S_{z'}^{be}.
}
\]

#### Virtual LLL/RRR block

The virtual piece is

\[
\frac13K_{JJSSJ}(w;x,y;z,z')f^{acb}
\left[
J_L^c(x)J_L^b(y)J_L^a(w)
-
J_R^c(x)J_R^b(y)J_R^a(w)
\right].
\]

Define

\[
\boxed{
V^{cba}(x,y,w)
=
\frac13\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')f^{acb}.
}
\]

The \(1/3\) factor must be included and tested.

---

### 4. Kernel symmetry and antisymmetry

KLM states that \(K_{JJSSJ}\) is antisymmetric under simultaneous interchange

\[
x\leftrightarrow y,
\qquad
z\leftrightarrow z'.
\]

Therefore the synthetic diagnostic kernels should include at least:

\[
K(w;x,y;z,z')=-K(w;y,x;z',z).
\]

Also test unconstrained and alternative symmetries if useful, but the KLM-like simultaneous antisymmetry must be the primary diagnostic.

The report must state which symmetry convention was used.

---

### 5. Expected density-derivative requirements

Like \(K_{JJSJ}\), the \(K_{JJSSJ}\) term is cubic in charge generators. Thus it requires:

\[
s_A=L_A\log W
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
K_{JJSSJ}\text{ requires Hessian-score.}
}
\]

---

### 6. Scope limitation

This workflow validates distinct-site ordered-current structure only.

Do not claim coincident-site ordering is resolved.  
Explicitly state that coincident-site commutators must be handled in a separate workflow:

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

---

## Phase 2: Extend Cubic Utility Module

Update:

```text
src/nlo_current/three_generator_terms.py
```

Add small dense validation utilities only. Do not implement production-volume kernels.

Required functions:

```python
def kjjssj_A_LLR_from_kernel(S_adj, KJJSSJ, f):
    """
    Build A_LLR[x,y,w,d,e,a] =
        sum_{z,z'} KJJSSJ[w,x,y,z,zp] f[a,c,b]
            S_adj[z,d,c] S_adj[zp,e,b].
    """

def kjjssj_B_LRR_from_kernel(S_adj, KJJSSJ, f):
    """
    Build B_LRR[w,x,y,a,d,e] =
        - sum_{z,z'} KJJSSJ[w,x,y,z,zp] f[a,c,b]
            S_adj[z,c,d] S_adj[zp,b,e].
    """

def kjjssj_V_virtual_from_kernel(KJJSSJ, f):
    """
    Build V[x,y,w,c,b,a] =
        (1/3) sum_{z,z'} KJJSSJ[w,x,y,z,zp] f[a,c,b].
    """

def synthetic_kjjssj_kernel(nsite, rng, klm_antisym=True):
    """
    Generate synthetic KJJSSJ[w,x,y,z,zp].

    If klm_antisym=True, impose:
        K(w;x,y;z,z') = -K(w;y,x;z',z).
    """
```

If `three_generator_terms.py` already has generic helpers from \(K_{JJSJ}\), reuse them.

Do not alter the \(K_{JJSJ}\) functions except for harmless shared refactors. If refactoring is needed, keep all existing \(K_{JJSJ}\) tests passing.

---

## Phase 3: Add \(K_{JJSSJ}\) Coefficient Tests

Create:

```text
tests/nlo_current/test_kjjssj_coefficients.py
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

and similarly for:

- \(B_{LRR}[w,x,y,a,d,e]\);
- \(V[x,y,w,c,b,a]\).

### Test 2: realness

For real synthetic \(K_{JJSSJ}\), real adjoint Wilson lines, and real \(f^{abc}\), verify all coefficient blocks are real.

### Test 3: KLM-like simultaneous antisymmetry

For synthetic KLM-like kernel:

\[
K(w;x,y;z,z')=-K(w;y,x;z',z),
\]

verify the generated kernel satisfies this relation numerically.

Then measure and report the induced combined symmetry/antisymmetry of \(A_{LLR}\) under:

\[
(x,d)\leftrightarrow(y,e).
\]

Do not force a guessed sign; report the measured behavior.

### Test 4: virtual \(1/3\) factor

Verify numerically that removing the \(1/3\) factor changes the virtual coefficient norm by exactly factor 3.

### Test 5: zero behavior for degenerate synthetic kernel

If \(K_{JJSSJ}=0\), all blocks should be zero.

---

## Phase 4: Add Cubic Current Tests for \(K_{JJSSJ}\)

Create:

```text
tests/nlo_current/test_kjjssj_cubic_current.py
```

Reuse finite-difference score/Hessian utilities from:

```text
src/nlo_current/finite_difference_scores.py
```

Tests:

### Test 1: LLR sign check with \(K_{JJSSJ}\)-like coefficient

Use random SU(3) Wilson lines, toy log density with nonzero Hessian-score, and the \(K_{JJSSJ}\) LLR coefficient.

Compare direct density-side operator:

\[
(\partial_YW)_{LLR}
=
+J_RJ_LJ_L(AW)
\]

against the proposed current divergence:

\[
-L_w J^{(w)}_{LLR}.
\]

Acceptance:

- finite-difference residual decreases with \(\epsilon\);
- use moderate eps values, e.g. \(10^{-3},10^{-4},10^{-5}\);
- tolerate \(10^{-4}\) to \(10^{-5}\) initially;
- report max residual.

### Test 2: LRR sign check with \(K_{JJSSJ}\)-like coefficient

Use the LRR coefficient and reuse the tested \(K_{JJSJ}\) LRR current lemma.

Compare direct adjoint density-side operator with proposed current divergence.

Report residual.

### Test 3: virtual LLL/RRR sign smoke test

Use \(V^{cba}\), include the \(1/3\) factor, and verify LLL/RRR signs using the existing cubic virtual test pattern from \(K_{JJSJ}\).

Report residuals.

### Test 4: Hessian-score contribution is nonzero

Using the toy coupled log density, evaluate the norm of the Hessian-score contribution for \(K_{JJSSJ}\) LLR and/or LRR blocks.

Assert it is nonzero above a small threshold, e.g.

```text
> 1e-10
```

This confirms \(K_{JJSSJ}\) genuinely requires Hessian-score.

---

## Phase 5: Diagnostic Script and Report

Create:

```text
scripts/nlo_current/check_kjjssj_cubic_requirements.py
```

It should:

1. Generate random SU(3) Wilson lines on \(N_{\rm site}=3\).
2. Generate synthetic \(K_{JJSSJ}\) kernels:
   - KLM-like simultaneous antisymmetry:
     \[
     K(w;x,y;z,z')=-K(w;y,x;z',z).
     \]
   - optionally unconstrained random kernel for stress comparison.
3. Build LLR, LRR, and virtual coefficient blocks.
4. Evaluate diagnostic norms of:
   - coefficient second-derivative term;
   - score-linear terms;
   - Hessian-score term;
   - score-product term;
   using a toy log density with nonzero Hessian-score.
5. Run or summarize sign-check residuals.
6. Save:

```text
reports/nlo_current/kjjssj_cubic_requirements_report.md
```

The report must include:

- random seed;
- \(N_{\rm site}\);
- kernel symmetry convention tested;
- coefficient block norms;
- evidence Hessian-score terms are nonzero;
- whether LLR/LRR/virtual signs passed tests;
- explicit statement that coincident-site commutators are not resolved;
- remaining uncertainties.

Conclusion should be one of:

\[
K_{JJSSJ}\text{ distinct-site current can be represented as score + Hessian-score with tested signs}.
\]

or

\[
K_{JJSSJ}\text{ still has unresolved distinct-site ordering/sign issues}.
\]

Do not overclaim.

---

## Phase 6: Update Three-Generator Summary

Update:

```text
docs/nlo_current/three_generator_sector_start.md
```

or create:

```text
docs/nlo_current/three_generator_sector_summary.md
```

It should contain:

1. \(K_{JJSJ}\) status:
   - LLR/LRR/virtual signs passed.
   - Hessian-score nonzero.
   - coincident-site commutators unresolved.
2. \(K_{JJSSJ}\) status:
   - coefficient blocks;
   - sign-test results;
   - Hessian-score result;
   - KLM-like kernel antisymmetry tested;
   - coincident-site commutators unresolved.
3. Shared conclusion:
   \[
   \text{NLO three-generator sector requires score + Hessian-score}.
   \]
4. What remains:
   - coincident-site commutator workflow;
   - dipole validation beyond skeleton;
   - eventual production-flow design.

---

## Phase 7: Dipole Validation Skeleton Update

Create:

```text
scripts/nlo_current/validate_dipole_kjjssj_skeleton.py
```

Add KLM Appendix A target notes for \(K_{JJSSJ}\), including the appearance of the \(\tilde K\) combination and the operator structures:

\[
N_c^3s(z,v)s(z',z)s(u,z')
\]

and traces like

\[
\mathrm{tr}\left[
S(v)S^\dagger(z)S(z')S^\dagger(u)S(z)S^\dagger(z')
\right].
\]

This is a skeleton only. Do not block workflow completion on full dipole validation.

---

## Phase 8: Acceptance Criteria

Stop when all are true:

1. Derivation doc exists:
   ```text
   docs/nlo_current/KJJSSJ_cubic_ordered_current.md
   ```

2. Cubic utility functions exist in:
   ```text
   src/nlo_current/three_generator_terms.py
   ```

3. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

4. Report exists:
   ```text
   reports/nlo_current/kjjssj_cubic_requirements_report.md
   ```

5. Three-generator summary exists or is updated:
   ```text
   docs/nlo_current/three_generator_sector_summary.md
   ```
   or
   ```text
   docs/nlo_current/three_generator_sector_start.md
   ```

6. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

7. No full NLO flow implemented.
8. No coincident-site commutator solution attempted beyond clearly documenting it as unresolved.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether the \(K_{JJSSJ}\) LLR sign passed.
4. Whether the \(K_{JJSSJ}\) LRR sign passed.
5. Whether the \(1/3\) virtual term was included and tested.
6. Whether Hessian-score terms were nonzero.
7. Which \(K_{JJSSJ}\) kernel symmetry convention was tested.
8. Whether any distinct-site sign/order issue remains.
9. Explicit caveat on coincident-site commutators.
10. Recommended next step:
    - if \(K_{JJSSJ}\) is stable, move to a dedicated coincident-site commutator workflow for \(K_{JJSJ}+K_{JJSSJ}\);
    - otherwise, resolve the remaining \(K_{JJSSJ}\) ordering/sign issue first.

Do not claim the NLO current is complete. This workflow only validates the distinct-site \(K_{JJSSJ}\) cubic ordered current.
