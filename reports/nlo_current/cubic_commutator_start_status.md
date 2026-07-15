# Cubic Commutator Workflow Start Status

This workflow is running in local no-git mode. The workspace has no `.git`
metadata, so no branch or clean-worktree gate is available.

## Previous cubic artifacts

- `docs/nlo_current/KJJSJ_cubic_ordered_current.md`: exists.
- `docs/nlo_current/KJJSSJ_cubic_ordered_current.md`: exists.
- `docs/nlo_current/three_generator_sector_summary.md`: exists.
- `src/nlo_current/su3_adjoint.py`: exists.
- `src/nlo_current/three_generator_terms.py`: exists.
- `src/nlo_current/finite_difference_scores.py`: exists.
- `tests/nlo_current/test_cubic_ordered_current.py`: exists.
- `tests/nlo_current/test_kjjsj_coefficients.py`: exists.
- `tests/nlo_current/test_kjjssj_coefficients.py`: exists.
- `tests/nlo_current/test_kjjssj_cubic_current.py`: exists.
- `scripts/nlo_current/check_kjjsj_cubic_requirements.py`: exists.
- `scripts/nlo_current/check_kjjssj_cubic_requirements.py`: exists.
- `reports/nlo_current/kjjsj_cubic_requirements_report.md`: exists.
- `reports/nlo_current/kjjssj_cubic_requirements_report.md`: exists.
- `reports/nlo_current/file_manifest.md`: exists.

## Pre-change tests

Command:

```bash
python3 -m pytest tests/nlo_current -q
```

Result:

```text
27 passed in 0.74s
```

## Implemented Lie derivative convention

The left perturbation convention is
\[
L^aF(U)=\frac{d}{d\epsilon}F(e^{i\epsilon t^a}U)\bigg|_{\epsilon=0}.
\]

With \([t^a,t^b]=if^{abc}t^c\), a direct finite-difference commutator check
confirmed
\[
[L^a,L^b]=f^{abc}L^c
\]
under this convention.

The adjoint convention remains
\[
S_A^{ab}=2\,{\rm ReTr}(t^a U t^b U^\dagger),
\]
and the verified right-to-left conversion is
\[
J_R^a=S_A^{ba}J_L^b.
\]

## Unresolved issue

Previous cubic workflows validated distinct-site signs only. Full lattice sums
include coincident sectors where
\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

This workflow validates the symbolic algebra and induced lower-order
commutator corrections. It does not implement a production NLO current.

