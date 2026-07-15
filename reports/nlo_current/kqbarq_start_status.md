# K_qbarq Workflow Start Status

This workflow is running in local no-git mode. The workspace has no `.git`
metadata, so no branch or clean-worktree gate is available.

## Previous artifacts

- `docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md`: exists.
- `src/nlo_current/su3_adjoint.py`: exists.
- `src/nlo_current/two_generator_terms.py`: exists.
- `tests/nlo_current/test_su3_adjoint_conventions.py`: exists.
- `scripts/nlo_current/check_kjssj_symmetry.py`: exists.
- `reports/nlo_current/kjssj_symmetry_report.md`: exists.
- `reports/nlo_current/file_manifest.md`: exists.

## Pre-change tests

Command:

```bash
python3 -m pytest tests/nlo_current -q
```

Result:

```text
7 passed in 0.38s
```

## Active conventions

The current small-lattice utilities use fundamental generators normalized by
\[
{\rm tr}(t^a t^b)=\delta^{ab}/2.
\]

The adjoint Wilson line convention is
\[
S_A^{ab}=2\,{\rm ReTr}(t^a U t^b U^\dagger).
\]

The left/right generator convention verified by finite differences is
\[
J_R^a=S_A^{ba}J_L^b.
\]

The previous ordered-current check also verified
\[
\sum_h L_y^h S_y^{hb}\approx 0
\]
under this convention.

