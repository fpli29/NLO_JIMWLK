# KJJSJ Analytic Cubic Derivative Validation Report

## Scope

This report validates only the non-production analytic coefficient derivatives
for the `KJJSJ` cubic sector. `KJJSSJ` remains explicitly pending. The
finite-difference backend is the reference oracle.

## Implemented Formulas

- `LC_K3`: `(LC_K3)^{AB} = L_C K3^{ABC}` from the same canonicalized KJJSJ normal-form tensors.
- `LB_K3`: `(LB_K3)^{AC} = L_B K3^{ABC}` with ordered combined-index contractions.
- `d2K3`: `d2K3^A = L_B L_C K3^{ABC}` using ordered second adjoint derivatives.
- `dK2_comm`: derivative of the KJJSJ quadratic commutator correction produced by canonicalization.
- `dK1_comm`: classified from canonical linear terms; not assumed zero.

The physical adapter supplies KLM-normalized cubic coefficients. The analytic
derivative code does not apply the `(-i)` normalization again.

## Full Synthetic Nonzero FD Step Scan

| seed | eps first | eps second | quantity | max abs | relative | real max | imag max | analytic norm | FD norm | expected zero | analytic s | FD s | speedup |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 20260756 | 2e-05 | 1e-03 | `dK2_comm` | 5.766e-12 | 5.946e-11 | 5.766e-12 | 0.000e+00 | 1.768e-01 | 1.768e-01 | False | 3.020e+01 | 5.259e+01 | 1.74 |
| 20260756 | 2e-05 | 1e-03 | `LC_K3` | 1.249e-11 | 5.298e-11 | 1.249e-11 | 0.000e+00 | 6.234e-01 | 6.234e-01 | False | 3.020e+01 | 5.259e+01 | 1.74 |
| 20260756 | 2e-05 | 1e-03 | `LB_K3` | 4.464e-12 | 5.526e-11 | 4.464e-12 | 0.000e+00 | 2.202e-01 | 2.202e-01 | False | 3.020e+01 | 5.259e+01 | 1.74 |
| 20260756 | 2e-05 | 1e-03 | `d2K3` | 3.985e-08 | 3.150e-07 | 3.985e-08 | 0.000e+00 | 1.675e-01 | 1.675e-01 | False | 3.020e+01 | 5.259e+01 | 1.74 |
| 20260756 | 2e-05 | 5e-04 | `dK2_comm` | 5.766e-12 | 5.946e-11 | 5.766e-12 | 0.000e+00 | 1.768e-01 | 1.768e-01 | False | 3.020e+01 | 5.292e+01 | 1.75 |
| 20260756 | 2e-05 | 5e-04 | `LC_K3` | 1.249e-11 | 5.298e-11 | 1.249e-11 | 0.000e+00 | 6.234e-01 | 6.234e-01 | False | 3.020e+01 | 5.292e+01 | 1.75 |
| 20260756 | 2e-05 | 5e-04 | `LB_K3` | 4.464e-12 | 5.526e-11 | 4.464e-12 | 0.000e+00 | 2.202e-01 | 2.202e-01 | False | 3.020e+01 | 5.292e+01 | 1.75 |
| 20260756 | 2e-05 | 5e-04 | `d2K3` | 9.884e-09 | 7.827e-08 | 9.884e-09 | 0.000e+00 | 1.675e-01 | 1.675e-01 | False | 3.020e+01 | 5.292e+01 | 1.75 |
| 20260758 | 2e-05 | 1e-03 | `dK2_comm` | 6.234e-12 | 3.857e-11 | 6.234e-12 | 0.000e+00 | 2.983e-01 | 2.983e-01 | False | 3.134e+01 | 5.768e+01 | 1.84 |
| 20260758 | 2e-05 | 1e-03 | `LC_K3` | 1.434e-11 | 5.092e-11 | 1.434e-11 | 0.000e+00 | 8.085e-01 | 8.085e-01 | False | 3.134e+01 | 5.768e+01 | 1.84 |
| 20260758 | 2e-05 | 1e-03 | `LB_K3` | 5.497e-12 | 5.780e-11 | 5.497e-12 | 0.000e+00 | 2.575e-01 | 2.575e-01 | False | 3.134e+01 | 5.768e+01 | 1.84 |
| 20260758 | 2e-05 | 1e-03 | `d2K3` | 3.605e-08 | 2.174e-07 | 3.605e-08 | 0.000e+00 | 2.247e-01 | 2.247e-01 | False | 3.134e+01 | 5.768e+01 | 1.84 |
| 20260758 | 2e-05 | 5e-04 | `dK2_comm` | 6.234e-12 | 3.857e-11 | 6.234e-12 | 0.000e+00 | 2.983e-01 | 2.983e-01 | False | 3.134e+01 | 5.258e+01 | 1.68 |
| 20260758 | 2e-05 | 5e-04 | `LC_K3` | 1.434e-11 | 5.092e-11 | 1.434e-11 | 0.000e+00 | 8.085e-01 | 8.085e-01 | False | 3.134e+01 | 5.258e+01 | 1.68 |
| 20260758 | 2e-05 | 5e-04 | `LB_K3` | 5.497e-12 | 5.780e-11 | 5.497e-12 | 0.000e+00 | 2.575e-01 | 2.575e-01 | False | 3.134e+01 | 5.258e+01 | 1.68 |
| 20260758 | 2e-05 | 5e-04 | `d2K3` | 9.071e-09 | 5.441e-08 | 9.071e-09 | 0.000e+00 | 2.247e-01 | 2.247e-01 | False | 3.134e+01 | 5.258e+01 | 1.68 |

## Per-Block Residuals

The block scan uses seed `20260756`, `eps_first=2e-5`, and
`eps_second=5e-4`.

| block | quantity | max abs | relative | real max | imag max | analytic norm | FD norm | expected zero |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `LLR` | `dK2_comm` | 2.700e-12 | 6.683e-11 | 2.700e-12 | 0.000e+00 | 9.339e-02 | 9.339e-02 | False |
| `LLR` | `LC_K3` | 6.311e-12 | 5.336e-11 | 6.311e-12 | 0.000e+00 | 3.117e-01 | 3.117e-01 | False |
| `LLR` | `LB_K3` | 2.297e-12 | 6.198e-11 | 2.297e-12 | 0.000e+00 | 1.101e-01 | 1.101e-01 | False |
| `LLR` | `d2K3` | 4.910e-09 | 7.782e-08 | 4.910e-09 | 0.000e+00 | 8.376e-02 | 8.376e-02 | False |
| `LRR` | `dK2_comm` | 3.281e-12 | 6.111e-11 | 3.281e-12 | 0.000e+00 | 1.000e-01 | 1.000e-01 | False |
| `LRR` | `LC_K3` | 5.825e-12 | 5.406e-11 | 5.825e-12 | 0.000e+00 | 3.117e-01 | 3.117e-01 | False |
| `LRR` | `LB_K3` | 2.004e-12 | 5.649e-11 | 2.004e-12 | 0.000e+00 | 1.101e-01 | 1.101e-01 | False |
| `LRR` | `d2K3` | 4.945e-09 | 7.825e-08 | 4.945e-09 | 0.000e+00 | 8.376e-02 | 8.376e-02 | False |
| `virtual_LLL` | `dK2_comm` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `virtual_LLL` | `LC_K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `virtual_LLL` | `LB_K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `virtual_LLL` | `d2K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `virtual_RRR` | `dK2_comm` | 5.529e-13 | 1.000e+00 | 5.529e-13 | 0.000e+00 | 0.000e+00 | 7.451e-13 | False |
| `virtual_RRR` | `LC_K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `virtual_RRR` | `LB_K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `virtual_RRR` | `d2K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |

## Velocity and Closure

| diagnostic | value |
|---|---:|
| velocity max residual vs FD | 1.646e-09 |
| velocity relative residual vs FD | 2.651e-08 |
| direct density operator | `(0.02881954217407266+0j)` |
| current divergence | `(0.028819543786385085+0j)` |
| closure absolute residual | 1.612e-09 |
| closure relative residual | 2.797e-08 |
| omit-Hessian absolute residual | 1.198e-04 |

## Physical Two-Site Expected-Zero Smoke Check

For the two-site physical coordinate setup, the diagnostic unbarred physical
`KJJSJ` array is structurally zero. This smoke check exercises the physical
adapter/backend path but does not replace the synthetic nonzero oracle scan.

| quantity | max abs | relative | real max | imag max | analytic norm | FD norm | expected zero |
|---|---:|---:|---:|---:|---:|---:|---|
| `dK2_comm` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `LC_K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `LB_K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |
| `d2K3` | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | True |

## Backend Status

- `backend="analytic"` uses no global FD fallback for `KJJSJ`.
- `KJJSSJ` remains pending and must raise under `backend="analytic"`.
- No complex-to-real cast is required by the KJJSJ analytic path; tested
  KLM-normalized derivative arrays are real to numerical precision.

## Commands Run

```text
python3 scripts/nlo_current/check_kjjsj_analytic_cubic_derivatives.py
KJJSJ analytic cubic validation report written; best d2K3 max residual 9.071e-09

python3 -m pytest tests/nlo_current/test_analytic_cubic_derivatives.py -q
8 passed in 269.26s (0:04:29)

python3 -m pytest tests/nlo_current -q
151 passed in 330.76s (0:05:30)
```
