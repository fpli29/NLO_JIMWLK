# K_JJSSJ Workflow Start Status

This workflow is running in local no-git mode. The workspace has no `.git`
metadata, so no branch or clean-worktree gate is available.

## Previous artifacts

- `docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md`: exists.
- `docs/nlo_current/Kqbarq_ordered_current.md`: exists.
- `docs/nlo_current/two_generator_sector_summary.md`: exists.
- `docs/nlo_current/KJJSJ_cubic_ordered_current.md`: exists.
- `docs/nlo_current/three_generator_sector_start.md`: exists.
- `src/nlo_current/su3_adjoint.py`: exists.
- `src/nlo_current/two_generator_terms.py`: exists.
- `src/nlo_current/three_generator_terms.py`: exists.
- `src/nlo_current/finite_difference_scores.py`: exists.
- `tests/nlo_current/test_su3_adjoint_conventions.py`: exists.
- `tests/nlo_current/test_kqbarq_coefficient.py`: exists.
- `tests/nlo_current/test_cubic_ordered_current.py`: exists.
- `tests/nlo_current/test_kjjsj_coefficients.py`: exists.
- `scripts/nlo_current/check_kjssj_symmetry.py`: exists.
- `scripts/nlo_current/check_kqbarq_symmetry.py`: exists.
- `scripts/nlo_current/check_kjjsj_cubic_requirements.py`: exists.
- `reports/nlo_current/kjssj_symmetry_report.md`: exists.
- `reports/nlo_current/kqbarq_symmetry_report.md`: exists.
- `reports/nlo_current/kjjsj_cubic_requirements_report.md`: exists.
- `reports/nlo_current/file_manifest.md`: exists.

## Pre-change tests

Command:

```bash
python3 -m pytest tests/nlo_current -q
```

Result:

```text
19 passed in 0.41s
```

## Active conventions

The small-lattice utilities use
\[
S_A^{ab}=2\,{\rm ReTr}(t^a U t^b U^\dagger),
\]
with \({\rm tr}(t^a t^b)=\delta^{ab}/2\).

The verified right-to-left conversion is
\[
J_R^a=S_A^{ba}J_L^b.
\]

The left perturbation convention is
\[
L^aF(U)=\frac{d}{d\epsilon}F(e^{i\epsilon t^a}U)\bigg|_{\epsilon=0}.
\]

## K_JJSJ caveat carried forward

The \(K_{JJSJ}\) workflow validated distinct-site ordered signs only.
Coincident-site commutators remain unresolved:
\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

This \(K_{JJSSJ}\) workflow keeps the same distinct-site scope.

