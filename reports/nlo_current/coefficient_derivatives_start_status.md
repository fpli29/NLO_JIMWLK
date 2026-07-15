# Coefficient-Derivative Backend Start Status

## Mode

- no_git_mode: true
- production_evolution_scope: not modified
- score_or_hessian_training_scope: not implemented
- backend_scope: dense small-lattice finite-difference diagnostics only

## Previous Skeleton Files

| file | exists |
|---|---:|
| `docs/nlo_current/nlo_current_skeleton_design.md` | yes |
| `docs/nlo_current/nlo_current_map_status.md` | yes |
| `docs/nlo_current/cubic_current_with_commutator_corrections.md` | yes |
| `src/nlo_current/nlo_current_skeleton.py` | yes |
| `src/nlo_current/nlo_velocity_evaluator.py` | yes |
| `src/nlo_current/synthetic_kernels.py` | yes |
| `src/nlo_current/finite_difference_scores.py` | yes |
| `src/nlo_current/su3_adjoint.py` | yes |
| `tests/nlo_current/test_nlo_current_skeleton.py` | yes |
| `tests/nlo_current/test_nlo_sector_assembly.py` | yes |
| `scripts/nlo_current/build_nlo_current_skeleton_demo.py` | yes |
| `reports/nlo_current/nlo_current_skeleton_demo_report.md` | yes |
| `reports/nlo_current/file_manifest.md` | yes |

## Baseline Tests

Before this workflow, the existing suite was run:

```text
python3 -m pytest tests/nlo_current -q
45 passed in 1.30s
```

## Current Velocity Evaluator Contract

`evaluate_velocity_from_terms(...)` accepts explicit dense coefficient-derivative arrays:

- `dK2`: shape `(D,)`, representing \(L_BK_2^{AB}\);
- `dK3_first["LC_K3_ABC"]`: shape `(D,D)`, representing \(L_CK_3^{ABC}\);
- `dK3_first["LB_K3_ABC"]`: shape `(D,D)`, representing \(L_BK_3^{ABC}\);
- `d2K3`: shape `(D,)`, representing \(L_BL_CK_3^{ABC}\).

When any are omitted, the evaluator treats them as zero and records warnings.

## Scope Statement

This workflow only adds a finite-difference diagnostic backend for tiny dense
lattices. It is not a production coefficient-derivative implementation.
