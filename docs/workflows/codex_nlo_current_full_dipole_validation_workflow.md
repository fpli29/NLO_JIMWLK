# Codex Workflow: Full Dipole Validation for All NLO JIMWLK Current Sectors

## Purpose

Continue from the completed NLO current skeleton and coefficient-derivative diagnostic workflows.

The goal is to validate the observable-side Hamiltonian action on the dipole against the KLM Appendix A results, sector by sector:

\[
K_{JSJ},\quad K_{JSSJ},\quad K_{q\bar q},\quad K_{JJSJ},\quad K_{JJSSJ}.
\]

The validation observable is

\[
s(u,v)=\frac{1}{N_c}\mathrm{tr}\left[U^\dagger(u)U(v)\right],\qquad N_c=3.
\]

The target is

\[
-H_{\rm sector}s(u,v).
\]

Use **unbarred KLM kernels** for singlet dipole validation. Barred/nonsinglet kernels are not the target of this workflow.

Do **not** implement production evolution. Do **not** train score/Hessian-score models. Do **not** optimize for large lattices. This is a non-production validation workflow only.

---

## Critical Rule

Do **not** invent Appendix A formulas.

If the exact KLM Appendix A expression for a sector is not available in the repo notes, existing skeletons, or local reference material, create a TODO target and mark that sector as:

```text
internal-consistency-only, Appendix A target missing
```

A sector may pass internal checks, but it must not be reported as “Appendix A passed” unless the exact Appendix A expression is implemented and compared.

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

If this is a Git repository, create a branch:

```bash
git checkout -b nlo-current-full-dipole-validation
```

If the worktree is dirty, stop and report dirty files before editing.

---

## Phase 0: Read Existing Artifacts

Inspect:

```text
docs/nlo_current/two_generator_sector_summary.md
docs/nlo_current/three_generator_sector_summary.md
docs/nlo_current/nlo_current_map_status.md
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
docs/nlo_current/Kqbarq_ordered_current.md
docs/nlo_current/KJJSJ_cubic_ordered_current.md
docs/nlo_current/KJJSSJ_cubic_ordered_current.md
docs/nlo_current/cubic_current_with_commutator_corrections.md
docs/nlo_current/coefficient_derivative_strategy.md

src/nlo_current/su3_adjoint.py
src/nlo_current/two_generator_terms.py
src/nlo_current/three_generator_terms.py
src/nlo_current/lie_word_algebra.py
src/nlo_current/cubic_commutator_terms.py
src/nlo_current/nlo_current_skeleton.py
src/nlo_current/nlo_velocity_evaluator.py
src/nlo_current/coefficient_derivatives.py
src/nlo_current/synthetic_kernels.py

scripts/nlo_current/validate_dipole_two_generator_terms.py
scripts/nlo_current/validate_dipole_kjjsj_skeleton.py
scripts/nlo_current/validate_dipole_kjjssj_skeleton.py
reports/nlo_current/file_manifest.md
```

Create/update:

```text
reports/nlo_current/dipole_validation_start_status.md
```

It should state:

- whether previous workflow files exist;
- whether previous tests still pass;
- whether no-git mode is active;
- where current dipole validation skeletons are;
- whether exact Appendix A target formulas are already available locally.

Run before changes:

```bash
python3 -m pytest tests/nlo_current -q
```

---

## Phase 1: Documentation — Dipole Validation Plan

Create:

```text
docs/nlo_current/full_dipole_validation_plan.md
```

Include:

### Observable and target

\[
s(u,v)=\frac{1}{N_c}\mathrm{tr}[U^\dagger(u)U(v)]
\]

and

\[
\frac{d}{dY}s(u,v)=-H_{\rm NLO}s(u,v).
\]

### Kernel choice

For singlet dipole validation, use unbarred kernels:

\[
K_{JSJ},K_{JSSJ},K_{q\bar q},K_{JJSJ},K_{JJSSJ}.
\]

### Two independent paths

Path A: direct observable-side generator action.

Path B: closed-form Appendix A expression.

A sector is fully validated only if Path A and Path B agree within numerical tolerance.

### Internal fallback

If exact Appendix A target is unavailable, implement internal checks only:

- direct analytic generator action versus finite-difference action;
- zero-kernel test;
- linearity in kernel;
- symmetry/antisymmetry stress tests;
- known subtraction identities such as the \(q\bar q\) \(z'=z\) vanishing identity.

Mark those sectors as pending exact Appendix A transcription.

---

## Phase 2: Implement Dipole Observable and Generator Actions

Create:

```text
src/nlo_current/dipole_observable.py
```

Implement:

```python
def dipole(U_fund, u, v):
    """Return s(u,v) = (1/Nc) Tr(U_u^dagger U_v)."""

def left_generator_action_on_dipole(U_fund, site, color, u, v, gens):
    """Analytic J_L action on s(u,v)."""

def right_generator_action_on_dipole(U_fund, site, color, u, v, gens):
    """Analytic J_R action on s(u,v)."""

def apply_generator_word_to_dipole(U_fund, word, u, v, gens, side_labels):
    """
    Apply an ordered word of generators to the dipole.
    word: tuple of (site, color)
    side_labels: tuple/list of "L" or "R", same length as word.
    Respect ordering exactly as provided.
    """

def fd_generator_action_on_dipole(U_fund, site, color, u, v, gens, side="L", eps=1e-6):
    """Finite-difference check for one generator action."""

def fd_generator_word_action_on_dipole(U_fund, word, u, v, gens, side_labels, eps=1e-5):
    """Nested finite-difference action for small validation only."""
```

Important:

- Follow existing left/right perturbation conventions in `su3_adjoint.py`.
- Be explicit about Hermitian versus anti-Hermitian generator convention.
- Analytic and finite-difference conventions must match existing tests.

---

## Phase 3: Tests for Dipole Generator Actions

Create:

```text
tests/nlo_current/test_dipole_generator_actions.py
```

Tests:

1. Dipole normalization: \(s(u,u)=1\).
2. Single \(J_L\) analytic action matches finite differences.
3. Single \(J_R\) analytic action matches finite differences.
4. Two-generator words match finite differences for `LL`, `LR`, `RL`, `RR`.
5. Three-generator words match finite differences for representative `LLR`, `LRR`, `LLL`, `RRR`.

Use loose but meaningful tolerances for nested finite differences and document eps sensitivity.

---

## Phase 4: Implement Direct Hamiltonian Action on the Dipole

Create:

```text
src/nlo_current/dipole_hamiltonian_action.py
```

Implement direct observable-side action functions returning \(-H_{\rm sector}s(u,v)\):

```python
def action_KJSJ_direct(U_fund, S_adj, KJSJ, u, v, gens):
    """
    Apply J_L J_L + J_R J_R - 2 J_L S J_R directly to the dipole.
    Return -H_KJSJ s(u,v).
    """

def action_KJSSJ_direct(U_fund, S_adj, KJSSJ, u, v, f, gens):
    """
    Apply the ordered KLM term:
        f f J_L S S J_R - Nc J_L S J_R.
    Return -H_KJSSJ s(u,v).
    """

def action_Kqbarq_direct(U_fund, S_adj, Kqbarq, u, v, gens):
    """Return -H_Kqbarq s(u,v) by direct generator action."""

def action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens):
    """
    Direct ordered action for LLR, LRR, and 1/3 virtual LLL-RRR pieces.
    Return -H_KJJSJ s(u,v).
    """

def action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens):
    """
    Direct ordered action for LLR, LRR, and 1/3 virtual LLL-RRR pieces.
    Return -H_KJJSSJ s(u,v).
    """

def action_all_sectors_direct(U_fund, S_adj, kernels, u, v, f, gens):
    """Return dict of sector direct actions and total."""
```

Important:

- Use explicit loops for clarity.
- Respect KLM observable-side ordering.
- Include the global sign \(dO/dY=-HO\).
- Use \(N_c=3\).
- Use unbarred kernels.

---

## Phase 5: Implement Appendix A Target Functions

Create:

```text
src/nlo_current/dipole_appendix_targets.py
```

Required interface:

```python
def target_KJSJ_appendix(U_fund, KJSJ, u, v):
    """Appendix A target for -H_KJSJ s(u,v)."""

def target_KJSSJ_appendix(U_fund, KJSSJ, u, v):
    """Raise NotImplementedError if exact formula is unavailable."""

def target_Kqbarq_appendix(U_fund, Kqbarq, u, v):
    """Raise NotImplementedError unless exact sign/normalization is confirmed."""

def target_KJJSJ_appendix(U_fund, KJJSJ, u, v):
    """Raise NotImplementedError unless exact formula is available."""

def target_KJJSSJ_appendix(U_fund, KJJSSJ, u, v):
    """Raise NotImplementedError unless exact formula is available."""

def appendix_target_available(sector_name):
    """Return bool."""
```

### Known KJSJ target

Implement the known target:

\[
-H_{JSJ}s(u,v)
=
2N_c\int_z
K_{JSJ}(u,v;z)
\left[
s(u,z)s(z,v)-s(u,v)
\right].
\]

Use the sign convention already established in previous docs/tests. If the direct action reveals an overall sign mismatch, stop and report.

### Other sectors

Only implement exact Appendix A formulas if they are available locally and confirmed. Otherwise raise `NotImplementedError` with a precise message.

Do not implement partial formulas as passing targets.

---

## Phase 6: Tests for Direct Action and Appendix Targets

Create:

```text
tests/nlo_current/test_dipole_hamiltonian_action.py
```

Tests:

### Test 1: zero kernels

For all sectors, zero kernels must give zero direct action.

### Test 2: linearity in kernel

For each sector:

\[
A[K_1+K_2]=A[K_1]+A[K_2].
\]

### Test 3: KJSJ Appendix A match

Compare:

```python
action_KJSJ_direct(...)
```

against:

```python
target_KJSJ_appendix(...)
```

for random SU(3), several dipoles, and synthetic kernels.

### Test 4: available Appendix targets

For each implemented target, compare direct action against target.

For unavailable targets, assert `NotImplementedError` is raised and ensure the sector is marked pending, not passed.

### Test 5: \(q\bar q\) subtraction consistency

Use the validated identity:

\[
2\mathrm{tr}(U^\dagger(z)t^aU(z)t^b)-S_A^{ab}(z)=0
\]

to build a synthetic case where the \(q\bar q\) coefficient vanishes and verify direct action vanishes.

### Test 6: virtual \(1/3\) sensitivity

For \(K_{JJSJ}\) and \(K_{JJSSJ}\), verify that removing the \(1/3\) virtual factor changes the direct action for generic kernels/Wilson lines.

---

## Phase 7: Full Dipole Validation Script

Create:

```text
scripts/nlo_current/full_dipole_validation.py
```

It should:

1. Generate random SU(3) Wilson lines for \(N_{site}=3\) or \(4\).
2. Generate synthetic unbarred kernels for all five sectors.
3. Choose several dipoles, for example `(0,1)`, `(1,2)`, `(0,2)`.
4. Compute direct actions for all sectors.
5. For each sector with Appendix A target available:
   - compute target;
   - compute absolute and relative residual.
6. For missing targets:
   - run internal consistency checks;
   - mark as `pending-target`.
7. Save:

```text
reports/nlo_current/full_dipole_validation_report.md
```

The report must include a table with:

```text
sector | direct action norm | appendix target available | residual | status
```

Statuses should be:

```text
passed | failed | pending-target
```

---

## Phase 8: Update Existing Dipole Skeletons

Update existing skeleton scripts if present:

```text
scripts/nlo_current/validate_dipole_two_generator_terms.py
scripts/nlo_current/validate_dipole_kjjsj_skeleton.py
scripts/nlo_current/validate_dipole_kjjssj_skeleton.py
```

They should either delegate to the new shared modules or clearly state they are superseded by:

```text
scripts/nlo_current/full_dipole_validation.py
```

Do not leave contradictory formulas across scripts.

---

## Phase 9: Documentation Update

Create/update:

```text
docs/nlo_current/dipole_validation_status.md
```

Summarize:

1. Direct generator-action validation status.
2. Appendix A target status per sector.
3. Which sectors fully passed Appendix A comparison.
4. Which sectors are pending exact target transcription.
5. Any sign/normalization issues found.
6. Whether the \(1/3\) virtual factor sensitivity test passed.
7. Remaining work before physical-kernel integration.

Also update:

```text
docs/nlo_current/nlo_current_map_status.md
```

Add a `Dipole validation status` section.

---

## Phase 10: Acceptance Criteria

Stop when all are true:

1. Dipole observable module exists:
   ```text
   src/nlo_current/dipole_observable.py
   ```
2. Direct Hamiltonian action module exists:
   ```text
   src/nlo_current/dipole_hamiltonian_action.py
   ```
3. Appendix target module exists:
   ```text
   src/nlo_current/dipole_appendix_targets.py
   ```
4. Tests pass:
   ```bash
   python3 -m pytest tests/nlo_current -q
   ```
5. Full validation report exists:
   ```text
   reports/nlo_current/full_dipole_validation_report.md
   ```
6. Dipole validation status doc exists:
   ```text
   docs/nlo_current/dipole_validation_status.md
   ```
7. Manifest updated:
   ```text
   reports/nlo_current/file_manifest.md
   ```
8. No production evolution code modified.
9. No score/Hessian-score model training implemented.
10. Any unavailable Appendix A target is clearly marked pending, not passed.

---

## Final Codex Response Required

At the end, summarize:

1. Files created/modified.
2. Tests run and results.
3. Whether analytic dipole generator actions match finite differences.
4. Whether direct Hamiltonian actions pass zero-kernel and linearity tests.
5. For each sector:
   - \(K_{JSJ}\)
   - \(K_{JSSJ}\)
   - \(K_{q\bar q}\)
   - \(K_{JJSJ}\)
   - \(K_{JJSSJ}\)

   report:
   - Appendix target available yes/no;
   - pass/fail/pending;
   - max residual if compared.

6. Whether the \(q\bar q\) subtraction consistency test passed.
7. Whether the \(1/3\) virtual sensitivity test passed for cubic sectors.
8. Any sign/normalization ambiguity.
9. Recommended next step:
   - exact transcription of missing Appendix A targets if any;
   - physical-kernel integration if all targets pass;
   - otherwise fix failed sector first.

Do not claim production readiness. This workflow validates dipole action only.
