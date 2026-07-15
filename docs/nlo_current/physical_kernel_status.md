# Physical Kernel Status

This document summarizes the non-production physical KLM kernel diagnostic
implementation. It does not claim production readiness.

## Implemented Kernels

Implemented unbarred singlet pointwise kernels:

| kernel | implementation | WORKNLO.tex lines |
|---|---|---|
| \(K_{JSJ}\) | `KJSJ_unbarred_value(...)` with explicit `KJSJIntegrationPolicy` | 324--332 |
| \(K_{JSSJ}\) | `KJSSJ_unbarred_value(...)` | 313--323 |
| \(K_{q\bar q}\) | `Kqbarq_unbarred_value(...)` | 301--306 |
| \(K_{JJSJ}\) | `KJJSJ_unbarred_value(...)` | 298--300 |
| \(K_{JJSSJ}\) | `KJJSSJ_unbarred_value(...)` | 288--297 |
| \(\widetilde K\) | `tilde_K_JJSSJ_unbarred_value(...)` | 307--311 |

## Pending Kernels

| kernel | status | reason |
|---|---|---|
| none | n/a | All five unbarred pointwise diagnostic interfaces exist. |

`KJSJ_unbarred_value(...)` still raises by default unless an explicit
`KJSJIntegrationPolicy` is supplied. This prevents accidental use of an
implicit finite-grid quadrature or regulator.

## Barred Kernels

Barred nonsinglet modifications are documented but not implemented:

- \(\bar K_{JSJ}\): `WORKNLO.tex` lines 367--370.
- \(\bar K_{JSSJ}\): `WORKNLO.tex` lines 374--375.
- \(\bar K_{q\bar q}\): `WORKNLO.tex` lines 378--383.
- \(K_{JJSJ}\) and \(K_{JJSSJ}\) remain unchanged in line 385.

## Symmetry Checks

`reports/nlo_current/physical_kernel_integration_report.md` records the
following pointwise symmetry residuals:

| diagnostic | residual |
|---|---:|
| `KJSJ_xy` | `0.0000000000000000e+00` |
| `Kqbarq_xy` | `0.0000000000000000e+00` |
| `Kqbarq_zzp` | `0.0000000000000000e+00` |
| `KJSSJ_xy` | `2.7105054312137611e-20` |
| `KJSSJ_zzp` | `4.0657581468206416e-20` |
| `KJJSJ_xy_antisym` | `8.2718061255302767e-25` |
| `KJJSSJ_simultaneous_antisym` | `6.7762635780344027e-21` |

## Singularity Policy

The coordinate kernel helpers expose:

```python
singularity_policy = "raise" | "nan" | "eps"
```

- `"raise"` is the default and raises on exact zero denominators.
- `"nan"` exposes singular dense-array entries as `np.nan`.
- `"eps"` requires an explicit positive `eps` and is a diagnostic regulator
  only.

No singularity policy here is a physical UV/IR prescription.

## Calling the Diagnostic Builder

```python
from nlo_current.physical_kernel_adapter import physical_kernels_for_skeleton

kernels = physical_kernels_for_skeleton(
    coords,
    Nc=3,
    nf=2,
    alpha_s=0.3,
    singularity_policy="eps",
    eps=1e-6,
    integration_policy=kjsj_policy,
)
```

The returned dict contains all five implemented physical kernels plus metadata.
The integration policy must specify the \(z'\) weights, \(\mu\),
coincident-point exclusions, principal-value label, subtraction label, and
finite-volume boundary label.

## Physical Dipole Recheck

`scripts/nlo_current/full_dipole_validation_physical_kernels.py` was run with
an explicit diagnostic finite-grid policy.

| sector | status | residual | relative residual |
|---|---|---:|---:|
| \(K_{JSJ}\) | passed | `7.7688469830766467e-12` | `2.8320504613732165e-09` |
| \(K_{JSSJ}\) | passed | `2.1438317028805034e-13` | `1.7081190426370902e-11` |
| \(K_{q\bar q}\) | passed | `3.6313016435989580e-13` | `1.9663547428238200e-10` |
| \(K_{JJSJ}\) | passed | `6.1080414466655290e-20` | `1.4244164473524333e-15` |
| \(K_{JJSSJ}\) | passed | `1.1408362482102489e-17` | `2.2558036620974433e-15` |

The \(K_{q\bar q}\) finite-grid mismatch was traced to the local compact
reduced trace-current target used in the previous diagnostic recheck. The
current recheck uses the exact WORKNLO trace-product expression from lines
1174--1177 and matches the direct action to roundoff. This resolves the
non-production dipole recheck; it does not define a production quadrature or
regulator.

## Cubic Dtype Convention

Raw physical \(K_{JJSJ}\) and \(K_{JJSSJ}\) kernels contain the explicit
\(-i\) factors in WORKNLO. The physical adapter now applies the explicit
diagnostic convention

```text
KLM-normalized cubic coefficient = (-1j) * raw physical cubic kernel
```

before dense normal-form assembly. Raw complex cubic tensors are still
supported for diagnostics, but the adapter output used by
`assemble_nlo_current_terms(...)` is real for the tested physical kernels.

## Positivity Caveat

Physical-kernel positivity checks are a future step. The existing
non-production Pawula toy diagnostic demonstrates positivity risk for finite
third-order generalized Fokker-Planck operators, but it does not prove physical
NLO JIMWLK positivity or non-positivity.

## Remaining Work Before Production

- Define a physical UV/IR regularization and subtraction strategy.
- Replace the diagnostic \(K_{JSJ}\) finite sum with a physical
  \(\int_{z'}\widetilde K\) prescription if production use is attempted.
- Implement barred/nonsinglet kernels if needed.
- Derive analytic coefficient derivatives for physical kernels.
- Design a score/Hessian-score estimator or contracted-Hessian strategy.
- Replace dense arrays with sparse/local structures for performance.
- Run physical validation beyond synthetic coordinate tests.

## Density-Side Closure Diagnostic

The physical kernels are now exercised in a tiny-lattice density-side closure
diagnostic:

- direct density operator:
  \(-L_A(K_1^AW)+\frac12L_AL_B(K_2^{AB}W)-\frac16L_AL_BL_C(K_3^{ABC}W)\);
- generalized-current side: \(-L_A(v^AW)\);
- positive finite-difference test densities;
- explicit diagnostic finite-grid policy.

The default physical scan is projected to active outer index `0` and
`KJSSJ`/`Kqbarq` for runtime control. The best projected residual is
`3.5453411040275995e-13` absolute. Sparse cubic diagnostic tests cover the
Hessian-score and raw-cubic normalization failure modes.

This does not change the production status: no production evolution,
regulator-independence proof, score/Hessian-score training, or physical
positivity claim is made.

## Analytic Coefficient Derivative Diagnostic

The physical kernels now feed a partial analytic/local coefficient-derivative
backend:

- `KJSJ`: analytic `dK2`;
- `KJSSJ`: analytic `dK2`;
- `Kqbarq`: analytic trace, subtraction, and full `dK2`.
- `KJJSJ`: analytic `dK2_comm`, `LC_K3`, `LB_K3`, and ordered `d2K3`
  validated in the dense diagnostic backend.

The finite-difference backend remains the oracle and is still required for
unproven sectors. Cubic physical derivative contractions for `KJJSSJ` remain
pending for a true analytic backend.
