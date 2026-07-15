# Physical Density-Side Closure Report

## Scope

This report covers the dense tiny-lattice, non-production density-side closure
diagnostic. It does not implement production evolution, train score/Hessian-
score models, claim regulator independence, or claim physical positivity.

The checked identity is

\[
{\cal G}_{\rm direct}[W]
=
{\cal G}_{\rm current}[W],
\]

with

\[
{\cal G}_{\rm direct}[W]
=
-L_A(K_1^AW)
+\frac12L_AL_B(K_2^{AB}W)
-\frac16L_AL_BL_C(K_3^{ABC}W),
\]

and

\[
{\cal G}_{\rm current}[W] = -L_A(v^A W).
\]

## Coordinate and Wilson-Line Setup

- Coordinate set for the physical projected scan:
  `[[0.0, 0.0], [1.0, 0.2]]`.
- Wilson-line configurations: deterministic random `SU(3)` matrices generated
  with fixed seeds in the script and tests.
- Active outer derivative direction for the physical scan: index `0`.

The active-index restriction keeps the dense finite-difference diagnostic small.
It is a projected tiny-lattice closure check, not a production-volume run.

## Physical Finite-Grid Policy

The physical scan used an explicit `KJSJIntegrationPolicy`:

- weights: `np.ones(2) / 2`;
- `mu = 1.3`;
- excluded coincident labels: `("x", "y", "z")`;
- singularity policy: `eps`;
- `eps = 1e-6`.

This policy is diagnostic and is not a production regulator.

## Cubic Convention

The dense physical current uses the established convention:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

The raw-cubic toggle in unit tests intentionally restores an `i`-weighted
cubic block and produces the expected complex failure.

## Test Densities

Implemented positive densities:

- Density A: `single_link_trace`, \(W=\exp[\lambda_1 \mathrm{Re\,tr}(U_x)]\).
- Density B: `dipole_trace`, adding
  \(\lambda_2 \mathrm{Re\,tr}(U_x^\dagger U_y)\).
- Density C: `multilink_nonlinear`, adding a squared neighboring trace.
- Constant density: \(W=1\).

Scores and ordered Hessian-scores are finite-difference derivatives of the
same \(\log W\). No learned score or Hessian-score model is used.

## Backend

- Direct density operator: dense finite differences of the written normal-form
  operator.
- Current divergence: dense finite differences of \(v^A(U)W(U)\), with
  score, Hessian-score, physical terms, and coefficient derivatives recomputed
  under outer perturbations.
- Physical coefficient derivatives: finite-difference diagnostic backend.
- For non-cubic sector filters, the physical derivative wrapper now records
  structurally zero `K3` derivative contractions instead of running cubic
  finite-difference loops.

## Physical Projected Step Scan

Command:

```text
python3 scripts/nlo_current/check_physical_density_closure.py
```

Result: completed successfully in about 33 seconds.

| density | fd eps | direct | current | abs residual | rel residual |
|---|---:|---:|---:|---:|---:|
| `single_link_trace` | `2.0e-3` | `3.531011439833710e-05` | `3.530999012155974e-05` | `1.242767773591705e-10` | `1.7597931666067723e-06` |
| `single_link_trace` | `1.0e-3` | `3.530984073194448e-05` | `3.530981217813179e-05` | `2.8553812688857788e-11` | `4.0433238499793905e-07` |
| `single_link_trace` | `5.0e-4` | `3.530977204688982e-05` | `3.530977169235571e-05` | `3.5453411040275995e-13` | `5.020339861042377e-09` |
| `single_link_trace` | `2.5e-4` | `3.530975413044892e-05` | `3.530974044852224e-05` | `1.368192668535495e-11` | `1.937414982495373e-07` |
| `dipole_trace` | `2.0e-3` | `2.907326013790655e-05` | `2.907313233137010e-05` | `1.2780653644510042e-10` | `2.198013170165815e-06` |
| `dipole_trace` | `1.0e-3` | `2.907298209193007e-05` | `2.907295700020531e-05` | `2.509172475914992e-11` | `4.315301317842803e-07` |
| `dipole_trace` | `5.0e-4` | `2.907291271400246e-05` | `2.907292373318282e-05` | `1.1019180359839176e-11` | `1.8950936185857535e-07` |
| `dipole_trace` | `2.5e-4` | `2.907289469930574e-05` | `2.907286796336806e-05` | `2.673593768251492e-11` | `4.5980887442513227e-07` |

The best projected physical residual in this scan is
`3.5453411040275995e-13` absolute.

## Per-Sector Residuals

For `dipole_trace`, `fd_eps=1e-3`, active outer index `0`:

| sector | abs residual | rel residual |
|---|---:|---:|
| `KJSSJ` | `4.0532749140571867e-11` | `9.460754901963364e-08` |
| `Kqbarq` | `1.5347200235920483e-11` | `3.154119778504921e-08` |

A separate direct measurement with the same coordinate setup and seed
`20260731` gave:

- direct: `1.0418862199382197e-04`;
- current: `1.0418861169977465e-04`;
- absolute residual: `1.029404731215603e-11`;
- relative residual: `4.940101723057073e-08`.

## Toggle Diagnostics

Physical two-generator projected scan:

- baseline residual: `5.1535998504739083e-11`;
- omit coefficient derivatives: `3.322883299670487e-05`.

Synthetic sparse cubic closure, with nonzero ordered Hessian-score:

- baseline residual: `1.4251942364435465e-10`;
- omit Hessian-score residual: `4.976507581037093e-05`;
- remove cubic normalization residual: `6.872404894813686e-05`;
- raw-cubic current imaginary part: `4.859531230150793e-05`.

Synthetic commutator-correction toggle:

- complete sparse block residual: below `2e-9`;
- omitting the lower-order commutator block changes the current by more than
  `1e-7`.

Constant density:

- \(W=1\), \(s=0\), and \(H=0\);
- constant-coefficient synthetic closure gives zero direct and current values
  to the tested tolerance.

## Warnings and Caveats

- The physical step scan is projected to active outer index `0` and restricted
  to `KJSSJ` and `Kqbarq` for runtime control.
- A one-off full all-sector projected check was run during development with
  active outer index `0`; it gave absolute residual
  `7.041882246587583e-11` and relative residual
  `4.7325147231493296e-08`, but took about 242 seconds and is not part of the
  default unit suite.
- The cubic closure and failure toggles are validated by sparse synthetic
  tests rather than a full physical all-sector finite-difference scan.
- No physical positivity or regulator-independence claim is made.

## Test Commands

Commands run for this workflow:

```text
python3 -m pytest tests/nlo_current -q
```

Baseline result before edits:

```text
117 passed in 65.25s (0:01:05)
```

Final result:

```text
129 passed in 49.76s
```

Additional targeted checks:

```text
python3 -m pytest tests/nlo_current/test_physical_density_closure.py -q
12 passed in 13.71s

python3 -m pytest tests/nlo_current/test_physical_nlo_current.py tests/nlo_current/test_coefficient_derivatives.py tests/nlo_current/test_physical_density_closure.py -q
24 passed in 40.72s

python3 scripts/nlo_current/check_physical_density_closure.py
completed successfully in about 33 seconds
```

## Analytic-Derivative Closure Recheck

The analytic coefficient derivative backend is complete for two-generator
`dK2` sectors and for the `KJJSJ` cubic diagnostic derivatives. The projected
physical density closure recheck with `derivative_backend="analytic"` covers
`KJSSJ` plus `Kqbarq` and remains within the established closure window:

```text
tests/nlo_current/test_analytic_coefficient_derivatives.py
```

The `KJJSJ` analytic cubic closure check is synthetic and nonzero, using the
same normal-form/canonicalization path with analytic `dK2_comm`, `LC_K3`,
`LB_K3`, and ordered `d2K3`. From
`reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`:

- analytic-vs-FD velocity max residual: `1.646e-09`;
- analytic-vs-FD velocity relative residual: `2.651e-08`;
- density-closure absolute residual: `1.612e-09`;
- density-closure relative residual: `2.797e-08`;
- omit-Hessian residual: `1.198e-04`.

`KJJSSJ` closure remains on the FD or explicitly labeled `hybrid_local_fd`
path until its analytic cubic derivatives are implemented.

The full current test suite after adding KJJSJ analytic cubic derivatives
passed:

```text
python3 -m pytest tests/nlo_current -q
151 passed in 330.76s (0:05:30)
```

\[
\boxed{
\text{The non-production physical NLO density operator is numerically reproduced by the generalized score/Hessian-score current on the tested tiny lattice.}
}
\]

This is a controlled dense diagnostic closure test. It is not a production
evolution result, a regulator-independence proof, or a positivity proof.
