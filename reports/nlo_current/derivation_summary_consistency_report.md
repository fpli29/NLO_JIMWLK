# Derivation Summary Consistency Report

## Scope

This report records the Phase 4 checks for the documentation-only derivation
summary workflow. No physical kernels, production evolution code, or
score/Hessian-score model training were added in this workflow.

## File Checks

| check | result |
|---|---|
| Main summary exists: `docs/nlo_current/NLO_JIMWLK_generalized_current_derivation_summary.md` | pass |
| Executive summary exists: `docs/nlo_current/NLO_JIMWLK_generalized_current_executive_summary.md` | pass |
| Start status exists: `reports/nlo_current/derivation_summary_start_status.md` | pass |
| This consistency report exists | pass |

## Test Check

Command:

```bash
python3 -m pytest tests/nlo_current -q
```

Result:

```text
84 passed in 9.25s
```

## Appendix A Sector Status Check

`docs/nlo_current/dipole_validation_status.md` marks all five sectors as
implemented and passed:

| sector | status file result | status-file residual |
|---|---|---:|
| \(K_{JSJ}\) | passed | `1.6172982426630268e-16` |
| \(K_{JSSJ}\) | passed | `1.3877787807814457e-16` |
| \(K_{q\bar q}\) | passed | `3.236828524569469e-16` |
| \(K_{JJSJ}\) | passed | `1.1801832636420706e-15` |
| \(K_{JJSSJ}\) | passed | `1.9087382418229414e-15` |

The latest full validation report,
`reports/nlo_current/full_dipole_validation_report.md`, also marks all five
sectors passed:

| sector | latest full-validation residual |
|---|---:|
| \(K_{JSJ}\) | `1.6172982426630268e-16` |
| \(K_{JSSJ}\) | `1.3877787807814457e-16` |
| \(K_{q\bar q}\) | `5.6325111571547335e-17` |
| \(K_{JJSJ}\) | `1.1801832636420706e-15` |
| \(K_{JJSSJ}\) | `1.9087382418229414e-15` |

The only numerical difference between the status file and the latest full
validation report in this check is the recorded \(K_{q\bar q}\) residual. Both
values are near machine precision and both sources mark the sector passed. The
new derivation summaries cite the latest full-validation residuals.

## Caveat Checks

The main derivation summary explicitly states:

- dense small-lattice, non-production scope;
- no physical coordinate kernels are implemented;
- no production evolution code is included;
- no score or Hessian-score model training is included;
- coefficient derivatives are currently diagnostic finite differences only;
- physical kernels are a next step;
- score/Hessian-score strategy remains future work.

The executive summary repeats the same caveats and states that the work is not
production-ready.

## Production-Readiness Check

The summaries do not claim production readiness. They frame the current state as
validated dense synthetic algebra, direct-action tests, and non-production
diagnostics.

## Result

All Phase 4 consistency checks passed. The only noted source-status difference
is a benign historical residual mismatch for \(K_{q\bar q}\), with both sources
agreeing on pass/fail status.
