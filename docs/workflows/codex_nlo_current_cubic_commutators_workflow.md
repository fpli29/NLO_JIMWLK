# Codex Workflow: Coincident-Site Commutator Corrections for \(K_{JJSJ}+K_{JJSSJ}\)

## Purpose

Continue from the completed distinct-site cubic workflows for:

\[
K_{JJSJ},
\qquad
K_{JJSSJ}.
\]

Previous workflows established, in the distinct-site small-lattice scope:

1. \(K_{JJSJ}\) LLR/LRR/virtual signs passed finite-difference checks.
2. \(K_{JJSSJ}\) LLR/LRR/virtual signs passed finite-difference checks.
3. Both terms have nonzero Hessian-score contributions.
4. No distinct-site sign/order issue remains.
5. Coincident-site commutators remain unresolved.

This workflow handles the missing algebraic piece:

\[
[L_x^a,L_y^b]=0 \quad (x\neq y),
\]

but

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

In full lattice sums, cubic terms include coincident coordinate sectors such as:

\[
x=y,\qquad x=w,\qquad y=w,\qquad x=y=w.
\]

These sectors can generate lower-order commutator corrections when ordered cubic derivatives are mapped into a canonical current form.

The goal is to derive, implement, and test a safe canonicalization layer for coincident-site ordered Lie derivatives in the cubic sector.

Do **not** implement the full NLO flow.  
Do **not** train or add score/Hessian-score models.  
Do **not** modify production evolution code.  
Only implement commutator algebra utilities, tests, diagnostics, reports, and documentation for \(K_{JJSJ}+K_{JJSSJ}\).

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
git checkout -b nlo-current-cubic-commutators
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect:

```text
docs/nlo_current/KJJSJ_cubic_ordered_current.md
docs/nlo_current/KJJSSJ_cubic_ordered_current.md
docs/nlo_current/three_generator_sector_summary.md
src/nlo_current/su3_adjoint.py
src/nlo_current/three_generator_terms.py
src/nlo_current/finite_difference_scores.py
tests/nlo_current/test_cubic_ordered_current.py
tests/nlo_current/test_kjjsj_coefficients.py
tests/nlo_current/test_kjjssj_coefficients.py
tests/nlo_current/test_kjjssj_cubic_current.py
scripts/nlo_current/check_kjjsj_cubic_requirements.py
scripts/nlo_current/check_kjjssj_cubic_requirements.py
reports/nlo_current/kjjsj_cubic_requirements_report.md
reports/nlo_current/kjjssj_cubic_requirements_report.md
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/cubic_commutator_start_status.md
```

It should state:

- whether previous cubic workflow files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- the implemented Lie derivative convention;
- the unresolved issue:
  \[
  [L_x^a,L_x^b]=f^{abc}L_x^c
  \]
  for coincident sites.

Run before changes:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 1: Documentation — Coincident-Site Problem Statement

Create:

```text
docs/nlo_current/cubic_coincident_site_commutators.md
```

This document must include the following.

### 1. Why distinct-site validation is insufficient

The previous workflows validated cubic ordered currents assuming distinct sites, so that

\[
[L_x^a,L_y^b]=0
\quad\text{for}\quad x\neq y.
\]

But full lattice sums include coincident coordinates. On the same site,

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

Therefore ordered products such as

\[
L_x^aL_x^bL_w^c,
\qquad
L_x^aL_w^cL_x^b,
\qquad
L_x^aL_x^bL_x^c
\]

cannot be freely symmetrized or reordered without generating lower-order terms.

### 2. Canonical ordering

Define a canonical ordering for combined indices:

\[
A=(x,a).
\]

A simple canonical order is lexicographic:

\[
(x,a)<(y,b)
\]

if either \(x<y\), or \(x=y\) and \(a<b\).

The canonicalization rule for neighboring same-site derivatives is:

\[
L_x^aL_x^b
=
L_x^bL_x^a
+
f^{abc}L_x^c
\]

if the pair is swapped from \((a,b)\) to \((b,a)\).

For different sites:

\[
L_x^aL_y^b=L_y^bL_x^a,
\qquad x\neq y.
\]

### 3. Operator-level normal ordering

For an ordered derivative word

\[
L_A L_B L_C
\]

acting on a coefficient-density product \(F[U]=K[U]W[U]\), canonicalization should return:

\[
L_A L_B L_C F
=
\sum_{\alpha} c_\alpha L_{\alpha_1}L_{\alpha_2}L_{\alpha_3}F
+
\sum_{\beta} d_\beta L_{\beta_1}L_{\beta_2}F
+
\sum_{\gamma} e_\gamma L_{\gamma_1}F.
\]

The cubic piece contributes to the Hessian-score current.  
The induced quadratic and linear pieces contribute to lower-order current corrections.

### 4. Classification of commutator-induced terms

The output of canonicalization must be classified into:

\[
K_3^{ABC}
\quad\text{cubic},
\]

\[
K_{2,\rm comm}^{AB}
\quad\text{quadratic commutator correction},
\]

\[
K_{1,\rm comm}^{A}
\quad\text{linear commutator correction}.
\]

Then the density current contains:

\[
J^A_{\rm cubic}
=
\frac16L_BL_C(K_3^{ABC}W),
\]

\[
J^A_{\rm quad,comm}
=
-\frac12L_B(K_{2,\rm comm}^{AB}W),
\]

\[
J^A_{\rm lin,comm}
=
K_{1,\rm comm}^{A}W,
\]

with signs depending on the canonical density-side normal form. The sign convention must be verified by finite differences.

### 5. Scope

This workflow should produce a canonicalization and testing layer. It does not need to fold these corrections into production flow.

Explicitly state:

\[
\boxed{
\text{This workflow validates algebraic commutator corrections only.}
}
\]

---

## Phase 2: Implement Lie Derivative Word Algebra

Create:

```text
src/nlo_current/lie_word_algebra.py
```

Use small symbolic objects, not heavy CAS.

Represent a Lie derivative label as:

```python
(site: int, color: int)
```

Represent a word as a tuple of labels:

```python
((site1, color1), (site2, color2), ...)
```

Implement:

```python
def is_canonical_word(word):
    """
    Return True if the derivative word is in canonical lexicographic order.
    """

def canonicalize_word(word, f, n_color=8):
    """
    Canonicalize a word of length 1, 2, or 3 using:
        L_x^a L_x^b = L_x^b L_x^a + f^{abc} L_x^c
    when swapping same-site adjacent derivatives.

    For different sites, swapping produces no commutator.

    Return a dictionary:
        {canonical_word: coefficient}

    Example:
        canonicalize_word(((0,2),(0,1)), f)
    should return:
        {
          ((0,1),(0,2)): 1,
          ((0,c),): f[2,1,c] summed over c
        }
    with the sign consistent with:
        L^a L^b = L^b L^a + [L^a,L^b]
                = L^b L^a + f^{abc}L^c.
    """

def canonicalize_terms(terms, f, n_color=8):
    """
    Input:
        terms: dict {word: coefficient}
    Output:
        canonicalized dict {word: coefficient}
    Combine like terms and drop tiny coefficients.
    """

def split_by_order(terms):
    """
    Return:
        cubic_terms, quadratic_terms, linear_terms, scalar_terms
    based on word length.
    """
```

Important:

- The algorithm must recursively canonicalize until all words are ordered.
- Avoid infinite loops.
- Drop coefficients with abs < 1e-14.
- Include docstrings explaining the sign convention.

---

## Phase 3: Unit Tests for Lie Word Algebra

Create:

```text
tests/nlo_current/test_lie_word_algebra.py
```

Tests:

### Test 1: different-site commute

For \(x\neq y\),

\[
L_x^aL_y^b
=
L_y^bL_x^a.
\]

Canonicalization should only reorder, with no lower-order term.

### Test 2: same-site two-word commutator

For \(a>b\), check:

\[
L_x^aL_x^b
=
L_x^bL_x^a
+
f^{abc}L_x^c.
\]

Compare symbolic output.

### Test 3: Jacobi consistency on three same-site words

Canonicalize the same three-letter word through different swap paths and verify the result is path-independent within numerical tolerance. This is a nontrivial consistency check of the SU(3) structure constants.

### Test 4: already canonical word unchanged

A canonical word returns itself with coefficient 1 and no lower-order correction.

### Test 5: split_by_order

Input a mixed dictionary and verify correct separation into cubic/quadratic/linear terms.

Run:

```bash
python3 -m pytest tests/nlo_current/test_lie_word_algebra.py -q
```

---

## Phase 4: Finite-Difference Validation of Canonicalization

Create:

```text
tests/nlo_current/test_lie_word_canonicalization_fd.py
```

Use existing finite-difference tools from:

```text
src/nlo_current/finite_difference_scores.py
src/nlo_current/su3_adjoint.py
```

Define a nontrivial scalar test function, for example:

\[
F(U_0,U_1)=
\mathrm{ReTr}(Q_0U_0)
+
0.3\,\mathrm{ReTr}(U_0U_1^\dagger)
+
0.2\,\mathrm{ReTr}(Q_1U_1).
\]

Tests:

### Test 1: same-site second derivative canonicalization

Numerically compare:

\[
L_x^aL_x^bF
\]

computed directly by nested finite differences against canonicalized expression:

\[
L_x^bL_x^aF+f^{abc}L_x^cF.
\]

### Test 2: same-site third derivative canonicalization

Numerically compare:

\[
L_x^aL_x^bL_x^cF
\]

against the canonicalized combination of third-, second-, and first-derivative terms.

### Test 3: mixed coincident pattern \(x=y\neq w\)

Compare:

\[
L_x^aL_x^bL_w^cF
\]

against canonicalized expression.

### Test 4: mixed coincident pattern \(x=w\neq y\)

Compare:

\[
L_x^aL_y^bL_x^cF
\]

against canonicalized expression.

Use loose finite-difference tolerances initially. Report residuals for:

```text
eps = 1e-3, 1e-4, 1e-5
```

Accept if residuals are small and/or decrease over a stable range. If finite differences are noisy at \(1e-5\), document it and use the stable range.

---

## Phase 5: Map Cubic Blocks to Canonical Terms

Create:

```text
src/nlo_current/cubic_commutator_terms.py
```

This module should convert cubic coefficient blocks from \(K_{JJSJ}\) and \(K_{JJSSJ}\) into canonical word terms.

Required functions:

```python
def cubic_block_terms_from_LLR(C_left):
    """
    Convert an already-left-basis LLR coefficient C[x,y,w,d,e,h] into:
        coefficient * word
    word = ((x,d), (y,e), (w,h))
    """

def cubic_block_terms_from_LRR(C_left):
    """
    Convert an already-left-basis LRR coefficient C[w,x,y,a,p,q] into:
        coefficient * word
    word = ((w,a), (x,p), (y,q))
    """

def cubic_block_terms_from_virtual_LLL(V):
    """
    Convert virtual LLL coefficient to word terms:
        word = ((x,c), (y,b), (w,a))
    """

def cubic_block_terms_from_virtual_RRR(V_left):
    """
    Convert virtual RRR after right-to-left conversion to word terms.
    """

def canonicalize_cubic_block_terms(terms, f):
    """
    Apply lie_word_algebra.canonicalize_terms, then split by order.
    Return:
        cubic_terms, quadratic_comm_terms, linear_comm_terms
    """
```

If previous utilities do not yet provide fully left-basis cubic coefficients, add small helper functions to multiply by adjoint matrices for right generators.

Do not implement production memory-heavy tensors. Use dense small-lattice diagnostics only.

---

## Phase 6: Coefficient-Level Diagnostics for \(K_{JJSJ}\) and \(K_{JJSSJ}\)

Create:

```text
scripts/nlo_current/check_cubic_commutator_corrections.py
```

It should:

1. Generate random SU(3) Wilson lines for \(N_{\rm site}=2\) or \(3\).
2. Generate synthetic kernels:
   - \(K_{JJSJ}\) with the tested antisymmetric convention;
   - \(K_{JJSSJ}\) with KLM-like simultaneous antisymmetry.
3. Build left-basis cubic terms using existing coefficient builders.
4. Canonicalize them with `lie_word_algebra`.
5. Split into:
   - cubic canonical terms;
   - quadratic commutator terms;
   - linear commutator terms.
6. Save:

```text
reports/nlo_current/cubic_commutator_corrections_report.md
```

The report must include:

- random seed;
- \(N_{\rm site}\);
- number of raw cubic terms;
- number of canonical cubic terms;
- norm of canonical cubic coefficients;
- norm of quadratic commutator coefficients;
- norm of linear commutator coefficients;
- separate sections for \(K_{JJSJ}\) and \(K_{JJSSJ}\);
- explicit statement whether commutator corrections are nonzero in coincident sectors.

Expected conclusion:

\[
\boxed{
\text{commutator-induced lower-order terms are present unless numerically zero.}
}
\]

If they are zero in synthetic diagnostics, do not assume proof; report it as diagnostic.

---

## Phase 7: Current Formula Update Document

Create:

```text
docs/nlo_current/cubic_current_with_commutator_corrections.md
```

This document should explain how the final cubic current should be represented after canonicalization.

Use schematic form:

\[
\partial_YW
=
-L_A(K_1^AW)
+
\frac12L_AL_B(K_2^{AB}W)
-
\frac16L_AL_BL_C(K_3^{ABC}W).
\]

After canonicalization of cubic Hamiltonian pieces:

\[
K_3\rightarrow K_{3,\rm canonical},
\]

and commutators induce:

\[
K_2\rightarrow K_2+K_{2,\rm comm},
\]

\[
K_1\rightarrow K_1+K_{1,\rm comm}.
\]

Then

\[
J^A
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

Therefore

\[
v^A
=
K_1^A
-
\frac12\left[
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

State clearly:

- cubic canonical terms require Hessian-score;
- quadratic commutator terms require score;
- linear commutator terms require no density derivatives beyond \(W\);
- coefficient derivatives are derivatives of known Wilson-line functions.

---

## Phase 8: Finite-Difference End-to-End Coincident Tests

Create:

```text
tests/nlo_current/test_cubic_commutator_end_to_end.py
```

This is the most important test.

Use a small manually constructed coefficient and scalar \(F=K[U]W[U]\). Compare:

1. Direct ordered cubic operator on \(F\), including coincident sites.
2. Canonicalized expression produced by the symbolic word algebra.

Test the four patterns:

### Pattern A: \(x=y\neq w\)

\[
L_x^aL_x^bL_w^cF.
\]

### Pattern B: \(x=w\neq y\)

\[
L_x^aL_y^bL_x^cF.
\]

### Pattern C: \(y=w\neq x\)

\[
L_x^aL_y^bL_y^cF.
\]

### Pattern D: \(x=y=w\)

\[
L_x^aL_x^bL_x^cF.
\]

Acceptance:

- residuals are small for a stable finite-difference range;
- report max residuals;
- if one pattern is numerically unstable, mark as warning and include details in the report.

---

## Phase 9: Update Three-Generator Summary

Update:

```text
docs/nlo_current/three_generator_sector_summary.md
```

Add a section:

```text
Coincident-site commutator status
```

Include:

- canonical Lie-word ordering rule;
- tests passed;
- whether commutator-induced \(K_2\) and \(K_1\) corrections were nonzero;
- remaining limitations;
- statement that production flow still requires assembling these corrections into the final NLO current implementation.

---

## Phase 10: Acceptance Criteria

Stop when all are true:

1. Documentation exists:
   ```text
   docs/nlo_current/cubic_coincident_site_commutators.md
   docs/nlo_current/cubic_current_with_commutator_corrections.md
   ```

2. Algebra module exists:
   ```text
   src/nlo_current/lie_word_algebra.py
   ```

3. Commutator mapping module exists:
   ```text
   src/nlo_current/cubic_commutator_terms.py
   ```

4. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```

5. Diagnostic report exists:
   ```text
   reports/nlo_current/cubic_commutator_corrections_report.md
   ```

6. End-to-end coincident tests exist and pass or explicitly report numerical instability:
   ```text
   tests/nlo_current/test_cubic_commutator_end_to_end.py
   ```

7. Three-generator summary updated:
   ```text
   docs/nlo_current/three_generator_sector_summary.md
   ```

8. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```

9. No full NLO flow implemented.
10. No score/Hessian model training implemented.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether same-site two-word commutator tests passed:
   \[
   L_x^aL_x^b=L_x^bL_x^a+f^{abc}L_x^c.
   \]
4. Whether same-site and mixed-site cubic canonicalization tests passed.
5. Whether \(K_{JJSJ}\) commutator-induced quadratic/linear corrections were nonzero.
6. Whether \(K_{JJSSJ}\) commutator-induced quadratic/linear corrections were nonzero.
7. End-to-end coincident finite-difference residuals.
8. Any remaining sign/order ambiguity.
9. Recommended next step:
   - if commutator algebra is stable, assemble a non-production NLO current skeleton;
   - otherwise, fix the failed commutator/canonicalization cases first.

Do not claim the NLO current implementation is complete. This workflow validates coincident-site commutator algebra and its lower-order corrections only.
