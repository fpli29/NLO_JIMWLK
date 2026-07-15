# Physical NLO Current Assembly Report

## Scope

This report records the first dense, small-lattice, non-production assembly of
the physical unbarred NLO generalized current. It does not implement production
evolution, train score/Hessian-score models, claim regulator independence, or
claim physical positivity.

## Implemented API

Implemented in `src/nlo_current/physical_nlo_current.py`:

- `assemble_physical_K1(...)`
- `assemble_physical_K2(...)`
- `assemble_physical_K3(...)`
- `assemble_physical_terms(...)`
- `compute_physical_coefficient_derivatives(...)`
- `evaluate_physical_nlo_velocity(...)`

## Implemented Sectors

The wrapper uses the existing `physical_kernel_adapter`, so all validated
unbarred physical sectors are included:

| sector | normal-form role |
|---|---|
| `KJSJ` | `K2` |
| `KJSSJ` | ordered `K2` |
| `Kqbarq` | ordered `K2` |
| `KJJSJ` | `K3` plus commutator corrections |
| `KJJSSJ` | `K3` plus commutator corrections |

## Cubic Convention

Raw physical cubic kernels are complex because WORKNLO includes explicit
\(-i\) factors. The adapter supplies the dense skeleton with KLM-normalized
real cubic coefficients:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

Tests verify no imaginary residual remains after adapter normalization.

## Derivative Backend

Implemented backend selector:

- `finite_difference`: dense central finite differences for tiny diagnostics.
- `diagnostic`: alias for `finite_difference`.
- `analytic`: explicit `NotImplementedError`; no physical analytic derivative
  formulas are implemented yet.

The finite-difference backend computes:

- `dK2`: \(L_B K2^{A B}\)
- `dK3_first`: first-derivative contractions of \(K3^{A B C}\)
- `d2K3`: \(L_B L_C K3^{A B C}\)

## Validation Summary

`tests/nlo_current/test_physical_nlo_current.py` covers:

- LO/K2-only score-current limit.
- zero-score pure coefficient-derivative contribution.
- zero-`K3` second-order-current limit.
- cubic normalization with no imaginary residual.
- sector-by-sector sum equals full physical assembly.
- finite-difference backend and `diagnostic` alias consistency.
- explicit `analytic` backend status.

Targeted test result:

```text
7 passed in 52.11s
```

Full `tests/nlo_current` result after adding the physical-current wrapper:

```text
117 passed in 62.52s (0:01:02)
```

Full result after adding the density-side closure diagnostic:

```text
129 passed in 49.76s
```

Full result after adding the partial analytic coefficient-derivative backend:

```text
143 passed in 61.23s (0:01:01)
```

## Density-Side Closure Extension

The physical current assembly is now used by the non-production density-side
closure diagnostic:

- `src/nlo_current/physical_density_operator.py`
- `src/nlo_current/physical_current_divergence.py`
- `src/nlo_current/physical_density_closure.py`
- `tests/nlo_current/test_physical_density_closure.py`

The diagnostic compares the direct normal-form density operator with
\(-L_A(v^AW)\) for positive test densities. The physical projected scan uses
active outer index `0` and `KJSSJ`/`Kqbarq`; the best recorded residual is
`3.5453411040275995e-13` absolute. Sparse synthetic cubic tests cover nonzero
Hessian-score dependence and raw-cubic normalization failure.

The physical coefficient-derivative wrapper now skips structurally zero cubic
derivative loops for non-cubic sector filters while preserving the existing
finite-difference backend for cubic sectors.

## Analytic Coefficient Derivative Extension

`src/nlo_current/physical_coefficient_derivatives.py` adds a structured backend
with:

- `analytic`: FD-validated for two-generator `dK2` only;
- `finite_difference`: preserved reference oracle;
- `diagnostic`: FD alias;
- `hybrid_local_fd`: explicit mixed path for pending sectors such as
  `KJJSSJ`.

Implemented analytic sectors: `KJSJ`, `KJSSJ`, and `Kqbarq`. Pending analytic
sectors: `KJJSJ` and `KJJSSJ` cubic first/second derivative contractions.

The existing `evaluate_physical_nlo_velocity(...)` path accepts
`derivative_backend="analytic"` for non-cubic sector filters and raises rather
than silently falling back for cubic sectors.

## Remaining Assumptions

- Wilson-line score and Hessian-score are supplied by the caller.
- The `KJSJ` finite \(z'\) sum is controlled by an explicit
  `KJSJIntegrationPolicy`.
- `singularity_policy='eps'` remains a diagnostic finite-array regulator.
- Dense finite differences are only for tiny coordinate sets.

## Known Limitations

- No production evolution loop.
- No trained score/Hessian-score model.
- No analytic physical coefficient derivatives.
- No barred/nonsinglet kernels.
- No physical positivity claim; Pawula caveat remains.
- No regulator independence claim beyond the current diagnostics.
