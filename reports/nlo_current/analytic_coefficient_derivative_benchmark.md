# Analytic Coefficient Derivative Benchmark

This benchmark compares implemented analytic derivatives with the
finite-difference oracle on the smallest physical diagnostic setup. It
does not claim asymptotic production scaling. The physical two-site
`KJJSJ` kernel is expected-zero; nonzero KJJSJ residuals are reported in
`reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`.

| sector | analytic s | FD s | observed speedup | dK2 max residual | dK2 relative residual |
|---|---:|---:|---:|---:|---:|
| `KJSJ` | 1.066063e-02 | 4.673000e-02 | 4.38 | 0.000e+00 | 0.000e+00 |
| `KJSSJ` | 3.264204e-01 | 3.179942e-01 | 0.97 | 3.028e-14 | 4.491e-11 |
| `Kqbarq` | 3.505325e-02 | 4.248071e-02 | 1.21 | 2.873e-09 | 2.899e-06 |

## KJJSJ Physical Two-Site Backend Benchmark

| backend | time s | FD s | observed speedup | dK2 max | LC max | LB max | d2 max | fallback used |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `analytic` | 3.334854e+00 | 6.280168e+00 | 1.88 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | False |
| `hybrid_local_fd` | 6.130139e+00 | 6.280168e+00 | 1.02 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |

`KJJSJ` is analytic-complete for the validated diagnostic backend; `KJJSSJ` remains pending.