# Dipole Validation Start Status

## Mode

- no_git_mode: true
- production_evolution_scope: not modified
- score_or_hessian_training_scope: not implemented
- validation_scope: observable-side dipole action, dense small-lattice only

## Previous Workflow Files

| file | exists |
|---|---:|
| `docs/nlo_current/two_generator_sector_summary.md` | yes |
| `docs/nlo_current/three_generator_sector_summary.md` | yes |
| `docs/nlo_current/nlo_current_map_status.md` | yes |
| `docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md` | yes |
| `docs/nlo_current/Kqbarq_ordered_current.md` | yes |
| `docs/nlo_current/KJJSJ_cubic_ordered_current.md` | yes |
| `docs/nlo_current/KJJSSJ_cubic_ordered_current.md` | yes |
| `docs/nlo_current/cubic_current_with_commutator_corrections.md` | yes |
| `docs/nlo_current/coefficient_derivative_strategy.md` | yes |
| `src/nlo_current/su3_adjoint.py` | yes |
| `src/nlo_current/two_generator_terms.py` | yes |
| `src/nlo_current/three_generator_terms.py` | yes |
| `src/nlo_current/lie_word_algebra.py` | yes |
| `src/nlo_current/cubic_commutator_terms.py` | yes |
| `src/nlo_current/nlo_current_skeleton.py` | yes |
| `src/nlo_current/nlo_velocity_evaluator.py` | yes |
| `src/nlo_current/coefficient_derivatives.py` | yes |
| `src/nlo_current/synthetic_kernels.py` | yes |
| `scripts/nlo_current/validate_dipole_two_generator_terms.py` | yes |
| `scripts/nlo_current/validate_dipole_kjjsj_skeleton.py` | yes |
| `scripts/nlo_current/validate_dipole_kjjssj_skeleton.py` | yes |
| `reports/nlo_current/file_manifest.md` | yes |

## Baseline Tests

Before this workflow, the existing suite was run:

```text
python3 -m pytest tests/nlo_current -q
53 passed in 3.00s
```

## Current Dipole Validation Skeletons

The existing dipole scripts are scaffolds:

- `scripts/nlo_current/validate_dipole_two_generator_terms.py`
- `scripts/nlo_current/validate_dipole_kjjsj_skeleton.py`
- `scripts/nlo_current/validate_dipole_kjjssj_skeleton.py`

They contain notes and partial targets, but do not provide complete confirmed
Appendix A comparisons for all sectors.

## Appendix A Target Availability

- `KJSJ`: exact target available locally and explicitly stated in the workflow.
- `KJSSJ`: target notes exist, but the scaffold says precise normalization/sign
  still needs checking, so it is not treated as an exact available target.
- `Kqbarq`: target notes exist, but the scaffold says precise normalization/sign
  still needs checking, so it is not treated as an exact available target.
- `KJJSJ`: only partial TODO notes are available locally.
- `KJJSSJ`: only partial TODO notes are available locally.

Unavailable targets are marked `pending-target`, not passed.
