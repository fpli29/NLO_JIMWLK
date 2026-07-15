# KJJSJ Analytic Cubic Benchmark

## Scope

This benchmark is non-production and covers only the diagnostic \(K_{JJSJ}\)
analytic cubic coefficient-derivative backend. It does not claim asymptotic
production scaling.

## Nonzero Synthetic Oracle Timing

From `reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`, the
nonzero two-site synthetic KJJSJ oracle scan measured:

| seed | eps second | analytic s | FD s | observed speedup | full d2 max residual |
|---:|---:|---:|---:|---:|---:|
| `20260756` | `1e-3` | `3.020e+01` | `5.259e+01` | `1.74` | `3.985e-08` |
| `20260756` | `5e-4` | `3.020e+01` | `5.292e+01` | `1.75` | `9.884e-09` |
| `20260758` | `1e-3` | `3.134e+01` | `5.768e+01` | `1.84` | `3.605e-08` |
| `20260758` | `5e-4` | `3.134e+01` | `5.258e+01` | `1.68` | `9.071e-09` |

The analytic path still loops over ordered derivative indices in dense form.
The measured speedup is therefore modest on this tiny lattice. The point of
this diagnostic is sign/order validation, not performance optimization.

## Physical Two-Site Backend Timing

The physical two-site unbarred KJJSJ kernel is structurally zero in the current
diagnostic setup. It is still useful as a backend-routing check:

| backend | time s | FD s | observed speedup | dK2 max | LC max | LB max | d2 max | fallback used |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `analytic` | `3.334854e+00` | `6.280168e+00` | `1.88` | `0.000e+00` | `0.000e+00` | `0.000e+00` | `0.000e+00` | `False` |
| `hybrid_local_fd` | `6.130139e+00` | `6.280168e+00` | `1.02` | `0.000e+00` | `0.000e+00` | `0.000e+00` | `0.000e+00` | `True` |

## Conclusion

`KJJSJ` analytic cubic derivatives are validated for the dense diagnostic
backend and avoid global FD fallback under `backend="analytic"`. `KJJSSJ`
remains pending.

## Command Run

```text
python3 scripts/nlo_current/benchmark_coefficient_derivative_backends.py
KJJSJ analytic: time=3.3349e+00s fd=6.2802e+00s speedup=1.88 d2_max=0.000e+00 fallback=False
KJJSJ hybrid_local_fd: time=6.1301e+00s fd=6.2802e+00s speedup=1.02 d2_max=0.000e+00 fallback=True
```
