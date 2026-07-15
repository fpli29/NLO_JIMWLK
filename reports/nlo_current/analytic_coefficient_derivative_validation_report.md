# Analytic Coefficient Derivative Validation Report

## Scope

This report compares implemented analytic/local physical coefficient derivatives
against the preserved finite-difference oracle. It is non-production and does
not claim regulator independence or physical positivity.

## Implemented Analytic Sectors

- `KJSJ`: analytic `dK2`.
- `KJSSJ`: analytic `dK2`.
- `Kqbarq`: analytic trace, subtraction, and full `dK2`.
- `KJJSJ`: analytic `dK2_comm`, `LC_K3`, `LB_K3`, and ordered `d2K3`.

## Pending Sectors

- `KJJSSJ`: cubic `LC_K3`, `LB_K3`, and `d2K3` pending.

`KJJSSJ` is not marked analytic-complete. Use `finite_difference` or
`hybrid_local_fd` explicitly for that sector.

## Residual Table

| seed | fd eps | sector | dK2 max | dK2 rel | LC max | LB max | d2 max | analytic s | FD s | speedup |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 20260753 | 2e-05 | `KJSJ` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 1.975e-02 | 5.127e-02 | 2.60 |
| 20260753 | 2e-05 | `KJSSJ` | 9.227e-14 | 5.587e-11 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 3.219e-01 | 3.146e-01 | 0.98 |
| 20260753 | 2e-05 | `Kqbarq` | 2.881e-09 | 1.875e-06 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 3.454e-02 | 4.162e-02 | 1.21 |
| 20260753 | 1e-05 | `KJSJ` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 8.667e-03 | 2.813e-02 | 3.25 |
| 20260753 | 1e-05 | `KJSSJ` | 1.200e-13 | 3.348e-11 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 2.952e-01 | 3.177e-01 | 1.08 |
| 20260753 | 1e-05 | `Kqbarq` | 4.875e-09 | 2.928e-06 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 3.441e-02 | 4.157e-02 | 1.21 |
| 20260754 | 2e-05 | `KJSJ` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 8.653e-03 | 2.932e-02 | 3.39 |
| 20260754 | 2e-05 | `KJSSJ` | 6.632e-14 | 6.312e-11 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 2.958e-01 | 3.172e-01 | 1.07 |
| 20260754 | 2e-05 | `Kqbarq` | 4.333e-09 | 2.565e-06 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 3.452e-02 | 4.162e-02 | 1.21 |
| 20260754 | 1e-05 | `KJSJ` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 8.640e-03 | 2.888e-02 | 3.34 |
| 20260754 | 1e-05 | `KJSSJ` | 9.974e-14 | 6.938e-11 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 3.081e-01 | 3.338e-01 | 1.08 |
| 20260754 | 1e-05 | `Kqbarq` | 8.590e-09 | 4.372e-06 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 3.539e-02 | 4.326e-02 | 1.22 |

## Dtype and Complex Checks

The analytic paths preserve complex intermediate values. The tested physical
outputs are real after `np.real_if_close`; no `ComplexWarning` was emitted by
the validation tests.

## KJJSJ Cubic Validation

The dedicated KJJSJ report
`reports/nlo_current/kjjsj_analytic_cubic_validation_report.md` validates
nonzero synthetic dense KJJSJ tensors against the finite-difference oracle.

Best stable-window residuals from seed `20260758`, `eps_second=5e-4`:

- `dK2_comm` max residual: `6.234e-12`;
- `LC_K3` max residual: `1.434e-11`;
- `LB_K3` max residual: `5.497e-12`;
- ordered `d2K3` max residual: `9.071e-09`.

Velocity and closure diagnostics from the same KJJSJ validation:

- analytic-vs-FD velocity max residual: `1.646e-09`;
- analytic-vs-FD velocity relative residual: `2.651e-08`;
- density-closure absolute residual: `1.612e-09`;
- density-closure relative residual: `2.797e-08`.

The physical two-site KJJSJ smoke check is expected-zero and is used only to
exercise the physical adapter/backend routing.

## Closure Status

`tests/nlo_current/test_analytic_coefficient_derivatives.py` checks analytic
two-generator velocity agreement with the FD reference and projected density
closure for `KJSSJ` plus `Kqbarq`. `tests/nlo_current/test_analytic_cubic_derivatives.py`
checks KJJSJ velocity agreement and synthetic density closure.

The explicit readout for the tested projected closure setup is:

- direct: `4.0167047558776765e-04`;
- current: `4.0167052449343070e-04 - 9.286236696325818e-13 i`;
- absolute residual: `4.8914478634410496e-11`;
- relative residual: `6.088881138827276e-08`.

## Test Commands

```text
python3 -m pytest tests/nlo_current/test_analytic_lie_derivatives.py -q
7 passed in 0.17s

python3 -m pytest tests/nlo_current/test_analytic_coefficient_derivatives.py -q
7 passed in 11.61s

python3 -m pytest tests/nlo_current/test_analytic_lie_derivatives.py tests/nlo_current/test_analytic_coefficient_derivatives.py tests/nlo_current/test_physical_nlo_current.py -q
21 passed in 36.22s

python3 -m pytest tests/nlo_current -q
143 passed in 61.23s (0:01:01)

python3 -m pytest tests/nlo_current/test_analytic_cubic_derivatives.py -q
8 passed in 269.26s (0:04:29)

python3 -m pytest tests/nlo_current/test_analytic_coefficient_derivatives.py -q
7 passed in 14.98s

python3 scripts/nlo_current/check_kjjsj_analytic_cubic_derivatives.py
KJJSJ analytic cubic validation report written; best d2K3 max residual 9.071e-09

python3 scripts/nlo_current/benchmark_coefficient_derivative_backends.py
KJJSJ analytic physical two-site backend speedup 1.88 against FD

python3 -m pytest tests/nlo_current -q
151 passed in 330.76s (0:05:30)
```

## Backend Completeness

The analytic backend is partial. It is justified for the listed two-generator
`dK2` sectors and the KJJSJ cubic diagnostic derivatives only. It is not a
full analytic NLO coefficient derivative backend until KJJSSJ cubic
`LC_K3`, `LB_K3`, `d2K3`, and commutator-correction derivatives pass the same
FD oracle and closure checks.
