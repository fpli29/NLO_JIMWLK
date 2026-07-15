# Codex Workflow: Add \(K_{q\bar q}\) Two-Generator Current to NLO JIMWLK Generalized Probability Current

## Purpose

Continue from the completed ordered \(J_L A J_R\) workflow.

The previous workflow established:

1. The KLM convention is
   \[
   \frac{d}{dY}{\cal O}=-H{\cal O},
   \qquad
   \partial_YW=-H^\dagger W.
   \]
2. The corrected \(K_{JSJ}\) contribution is LO-like:
   \[
   v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
   \]
3. The ordered \(J_L A J_R\) current lemma passed finite-difference checks.
4. Numerically,
   \[
   L_y^hS_y^{hb}\approx0
   \]
   under the implemented convention.
5. \(K_{JSSJ}\) has an order-one antisymmetric component:
   \[
   \frac{\|C-C^T\|}{\|C\|}\sim1.3,
   \]
   so it must **not** be simplified to a pure symmetric score current. It should be kept as a generic ordered-current / second-order current with coefficient drift and possible commutator drift.

The current goal is to add the next two-generator NLO term:

\[
K_{q\bar q}.
\]

Do **not** implement the full NLO flow yet.  
Do **not** start the three-generator terms yet.  
Only complete the \(K_{q\bar q}\) derivation, utilities, tests, and symmetry report.

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
git checkout -b nlo-current-kqbarq-two-generator
```

If the worktree is dirty, stop and report the dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect these existing files from the previous workflow:

```text
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
src/nlo_current/su3_adjoint.py
src/nlo_current/two_generator_terms.py
tests/nlo_current/test_su3_adjoint_conventions.py
scripts/nlo_current/check_kjssj_symmetry.py
reports/nlo_current/kjssj_symmetry_report.md
reports/nlo_current/file_manifest.md
```

Then create/update a short status note:

```text
reports/nlo_current/kqbarq_start_status.md
```

It should state:

- whether the previous files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- the convention currently used for:
  \[
  S_A^{ab}=2\mathrm{ReTr}(t^a U t^b U^\dagger),
  \]
  and
  \[
  J_R^a=S_A^{ba}J_L^b.
  \]

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

before starting new changes.

---

## Phase 1: Add Derivation Notes for \(K_{q\bar q}\)

Create or update:

```text
docs/nlo_current/Kqbarq_ordered_current.md
```

This document must include the following derivation.

### 1. Hamiltonian structure

The quark-pair contribution in the KLM Hamiltonian has the two-generator structure

\[
H_{q\bar q}
=
\int_{x,y,z,z'}
K_{q\bar q}(x,y;z,z')
\left[
2J_L^a(x)
\mathrm{tr}\!\left(S^\dagger(z)t^aS(z')t^b\right)
J_R^b(y)
-
J_L^a(x)S_A^{ab}(z)J_R^b(y)
\right].
\]

For configuration-level current use the generalized/nonsinglet barred kernel:

\[
K_{q\bar q}\rightarrow \bar K_{q\bar q}.
\]

For singlet dipole validation, use the unbarred kernel.

### 2. Define the ordered \(J_L A J_R\) coefficient

Define

\[
A_{q\bar q}^{ab}(x,y;U)
=
\int_{z,z'}
\bar K_{q\bar q}(x,y;z,z')
\left[
2\mathrm{tr}\!\left(S^\dagger(z)t^aS(z')t^b\right)
-
S_A^{ab}(z)
\right].
\]

Then

\[
H_{q\bar q}
=
\int_{x,y}
A_{q\bar q}^{ab}(x,y;U)
J_L^a(x)J_R^b(y).
\]

### 3. Apply the ordered \(J_L A J_R\) lemma

From the previous workflow, for

\[
H_{LR}[A]=A^{ab}J_L^aJ_R^b
\]

with

\[
\partial_Y{\cal O}=-H_{LR}{\cal O},
\]

and with

\[
J_R^b(y)=S_y^{hb}L_y^h,
\]

the density-side current is

\[
J_{LR}^{(y,h)}
=
S_y^{hb}
L_x^a[A^{ab}W],
\]

assuming

\[
L_y^hS_y^{hb}=0
\]

under the repo convention.

Therefore

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

So for the quark term,

\[
\boxed{
v_{q\bar q}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{q\bar q}^{ab}
+
A_{q\bar q}^{ab}s_x^a
\right].
}
\]

### 4. State what this term needs

Because \(K_{q\bar q}\) is second order in charge generators,

\[
\boxed{
K_{q\bar q}\text{ needs the score }s_A,
\text{ but not the Hessian-score }L_As_B.
}
\]

However, as with \(K_{JSSJ}\), do **not** assume the coefficient is symmetric.  
The symmetry/asymmetry of the left-basis coefficient must be measured.

### 5. Left-basis coefficient

Convert to left basis:

\[
C_{q\bar q}^{(x,a)(y,h)}
=
A_{q\bar q}^{ab}(x,y;U)S_y^{hb}.
\]

Then decompose:

\[
C^{AB}=C^{(AB)}+C^{[AB]}.
\]

If \(C^{[AB]}\neq0\), keep generic second-order current / commutator drift.

---

## Phase 2: Add \(K_{q\bar q}\) Utilities

Update:

```text
src/nlo_current/two_generator_terms.py
```

Add functions:

```python
def kqbarq_A_from_kernel(U_fund, S_adj, Kqbarq, gens):
    """
    Build A_qbarq^{ab}(x,y):
    A^{ab}(x,y) = sum_{z,z'} K_qbarq(x,y;z,z')
        [2 Tr(S_fund(z)^dagger t^a S_fund(z') t^b) - S_adj[z,a,b]].

    Inputs:
        U_fund: array/list of fundamental Wilson lines, shape (Nsite,3,3)
        S_adj: adjoint Wilson lines, shape (Nsite,8,8)
        Kqbarq: dense kernel, shape (Nsite,Nsite,Nsite,Nsite)
        gens: fundamental SU(3) generators t^a

    Output:
        A: array, shape (Nsite,Nsite,8,8), indexed A[x,y,a,b].
    """

def kqbarq_C_left_from_A(A, S_adj):
    """
    Convert A^{ab}(x,y) J_L^a(x) J_R^b(y)
    to left-basis coefficient:
        C^{(x,a)(y,h)} = A[x,y,a,b] S_adj[y,h,b].
    Return dense matrix/tensor in the same format used for K_JSSJ.
    """
```

If `kjssj_C_left_from_A(A, S_adj)` already does exactly this, reuse it and add a documented alias rather than duplicating logic.

Also add a helper if useful:

```python
def qbarq_trace_block(Uz, Uzp, gens):
    """
    Return B^{ab}(z,z') = 2 Tr(Uz^\dagger t^a Uzp t^b).
    """
```

Use complex trace carefully. The coefficient should be real up to numerical roundoff. Return `np.real_if_close`.

---

## Phase 3: Algebraic Tests for \(K_{q\bar q}\)

Create:

```text
tests/nlo_current/test_kqbarq_coefficient.py
```

### Test 1: \(z'=z\) subtraction

For random SU(3) \(U_z\), verify:

\[
2\mathrm{tr}\!\left(S^\dagger(z)t^aS(z)t^b\right)
-
S_A^{ab}(z)
=
0
\]

within double-precision tolerance.

This checks consistency between the fundamental trace block and adjoint Wilson line convention.

Acceptance:

```text
max_abs_error < 1e-10
```

or justified tolerance.

### Test 2: coefficient reality

For random Wilson lines and a real synthetic kernel, verify

\[
A_{q\bar q}^{ab}
\]

is real up to roundoff.

### Test 3: left-basis conversion shape

Build \(A\) on \(N_{\rm site}=3\) and verify that

\[
C^{(x,a)(y,h)}
=
A^{ab}(x,y)S_y^{hb}
\]

has the expected dense matrix shape or tensor shape consistent with `sym_asym_parts`.

### Test 4: reuse ordered LR lemma

If the previous ordered LR finite-difference test was generic enough, ensure it is still passing.  
If not, add a small smoke test using a \(q\bar q\)-like coefficient \(A(U)\).

Run:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 4: Symmetry Report for \(C_{q\bar q}\)

Create:

```text
scripts/nlo_current/check_kqbarq_symmetry.py
```

It should:

1. Generate random SU(3) Wilson lines for \(N_{\rm site}=3\) or \(4\).
2. Generate a synthetic real kernel satisfying:
   \[
   K(x,y;z,z')=K(y,x;z,z')
   \]
   and optionally also:
   \[
   K(x,y;z,z')=K(x,y;z',z).
   \]
3. Build
   \[
   A_{q\bar q}^{ab}(x,y).
   \]
4. Convert to left-basis
   \[
   C^{(x,a)(y,h)}.
   \]
5. Compute:
   \[
   r_{\rm asym}=\frac{\|C-C^T\|}{\|C\|}.
   \]
6. Save:

```text
reports/nlo_current/kqbarq_symmetry_report.md
```

The report must include:

- random seed;
- \(N_{\rm site}\);
- kernel symmetries imposed;
- asymmetry ratio with \(x\leftrightarrow y\) symmetry only;
- asymmetry ratio with both \(x\leftrightarrow y\) and \(z\leftrightarrow z'\) symmetry, if tested;
- conclusion:
  - if order-one asymmetry, keep generic current / commutator drift;
  - if near-zero asymmetry, state that it appears symmetric in the synthetic test but still requires analytic confirmation.

Do not present the numerical test as a proof.

---

## Phase 5: Update Two-Generator Sector Summary

Create or update:

```text
docs/nlo_current/two_generator_sector_summary.md
```

It should summarize the current status:

### \(K_{JSJ}\)

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
(S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}),
\]

\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
\]

Status:

\[
\text{LO-like symmetric score current.}
\]

### \(K_{JSSJ}\)

\[
A_{JSSJ}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{JSSJ}
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}),
\]

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

Status:

\[
\text{score + coefficient drift; order-one antisymmetric component measured.}
\]

### \(K_{q\bar q}\)

\[
A_{q\bar q}^{ab}(x,y)
=
\int_{z,z'}
\bar K_{q\bar q}
\left[
2\mathrm{tr}(S^\dagger(z)t^aS(z')t^b)-S_A^{ab}(z)
\right],
\]

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

Status:

- fill in after symmetry report.

### Shared conclusion

All two-generator NLO pieces require at most the score \(s_A\).  
None of them require Hessian-score.

Hessian-score first appears in the three-generator sector:

\[
K_{JJSJ},
\qquad
K_{JJSSJ}.
\]

---

## Phase 6: Dipole Validation Skeleton Update

Update or create:

```text
scripts/nlo_current/validate_dipole_two_generator_terms.py
```

Add the \(K_{q\bar q}\) target expression from KLM Appendix A:

\[
-H_{q\bar q}s(u,v)
=
-\int_{z,z'}K_{q\bar q}(u,v;z,z')
\left(
N_c s(u,z')s(z,v)
-
\frac{1}{N_c^2}
\mathrm{tr}[S^\dagger(u)S(v)S^\dagger(z)S(z')]
\right.
\]

\[
\left.
-
\frac{1}{N_c^2}
\mathrm{tr}[S^\dagger(u)S(v)S^\dagger(z')S(z)]
+
\frac{1}{N_c}
s(u,v)s(z,z')
\right)
\]

up to the precise normalization/sign from the KLM appendix.

If the existing dipole utilities are not ready, leave it as a clear TODO scaffold.  
Do not block this workflow on full dipole validation.

---

## Phase 7: Acceptance Criteria

Stop when all are true:

1. The \(K_{q\bar q}\) derivation document exists:
   ```text
   docs/nlo_current/Kqbarq_ordered_current.md
   ```

2. The utilities exist or documented aliases exist:
   - `kqbarq_A_from_kernel`
   - `kqbarq_C_left_from_A` or alias to the generic \(A\to C\) converter
   - optional `qbarq_trace_block`

3. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

4. A report exists:
   ```text
   reports/nlo_current/kqbarq_symmetry_report.md
   ```

5. The two-generator sector summary exists:
   ```text
   docs/nlo_current/two_generator_sector_summary.md
   ```

6. The manifest is updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

7. No full NLO flow is implemented.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether the \(z'=z\) subtraction identity passed:
   \[
   2\mathrm{tr}(S^\dagger t^a S t^b)-S_A^{ab}=0.
   \]
4. The measured \(K_{q\bar q}\) asymmetry ratio.
5. Whether \(K_{q\bar q}\) appears symmetric or must be kept as a generic ordered-current / commutator-drift term.
6. Whether any sign or left/right convention issue remains.
7. Recommended next step:
   - if two-generator sector is stable, move to \(K_{JJSJ}\);
   - otherwise, fix the two-generator current representation first.

Do not claim that the NLO current is complete. This workflow only completes the \(K_{q\bar q}\) two-generator term and the two-generator sector summary.
