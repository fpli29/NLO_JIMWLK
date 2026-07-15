# Analytic Coefficient Derivative Start Status

## Scope

This workflow starts from the completed non-production physical density-side
closure validation. The goal is to add source-grounded analytic/local
coefficient-derivative diagnostics while preserving the finite-difference
backend as the reference oracle.

This is not production evolution, score/Hessian-score training, a regulator-
independence proof, or a positivity claim.

## Baseline Tests

Command:

```text
python3 -m pytest tests/nlo_current -q
```

Result:

```text
129 passed in 49.95s
```

## Current Finite-Difference Runtime

The full dense finite-difference coefficient-derivative path remains expensive
because it perturbs every color/site direction and, for cubic sectors, computes
nested finite differences. The current full suite baseline is roughly 50
seconds on this workspace. The previous density-closure report records a
one-off all-sector projected closure check taking about 242 seconds.

## Current Closure Residuals

From `reports/nlo_current/physical_density_closure_report.md`:

- best projected physical closure residual: `3.5453411040275995e-13`;
- `KJSSJ` projected residual at `fd_eps=1e-3`: `4.0532749140571867e-11`;
- `Kqbarq` projected residual at `fd_eps=1e-3`: `1.5347200235920483e-11`;
- one-off all-sector projected residual: `7.041882246587583e-11`.

## Current Dtype and Cubic Convention

The physical adapter uses:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

The dense normal-form tensors consume the KLM-normalized coefficients. Complex
intermediate diagnostics are allowed; silent complex-to-real casting is not.

## Current Physical Finite-Grid Policy

The existing physical diagnostics require an explicit `KJSJIntegrationPolicy`.
The standard tiny-lattice tests use equal weights, `mu=1.3`, excluded labels
`("x", "y", "z")`, and `singularity_policy="eps"` with `eps=1e-6`.

This remains a diagnostic finite-grid policy, not a production regulator.

## Available Derivative Contractions

The finite-difference oracle currently provides:

- `dK2`: \(L_BK_2^{AB}\);
- `LC_K3_ABC`: \(L_CK_3^{ABC}\);
- `LB_K3_ABC`: \(L_BK_3^{ABC}\);
- `d2K3`: \(L_BL_CK_3^{ABC}\).

No analytic physical coefficient-derivative backend was available at baseline.

