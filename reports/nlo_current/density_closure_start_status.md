# Density Closure Start Status

## Scope

This workflow is a dense, tiny-lattice, non-production density-side closure
validation for the physical NLO generalized-current normal form. It does not
implement production evolution, train score/Hessian-score models, optimize for
large lattices, claim physical positivity, or claim regulator independence.

## Baseline Audit

The requested source, test, documentation, and report files were present at the
start of this workflow, including the physical current assembly wrapper,
finite-difference coefficient derivative backend, physical kernel adapter,
cubic convention layer, and prior physical-kernel integration reports.

## Baseline Tests

Command:

```text
python3 -m pytest tests/nlo_current -q
```

Result:

```text
117 passed in 65.25s (0:01:05)
```

The expected count in the workflow was approximate; the actual clean baseline
count is 117.

## Physical Assembly Status

The current non-production physical assembly exposes:

- `assemble_physical_K1(...)`
- `assemble_physical_K2(...)`
- `assemble_physical_K3(...)`
- `assemble_physical_terms(...)`
- `compute_physical_coefficient_derivatives(...)`
- `evaluate_physical_nlo_velocity(...)`

The implemented unbarred physical sectors are:

- `KJSJ`
- `KJSSJ`
- `Kqbarq`
- `KJJSJ`
- `KJJSSJ`

## Derivative Backend Status

The available physical coefficient-derivative backends are:

- `finite_difference`: dense central finite differences for tiny diagnostics.
- `diagnostic`: alias for `finite_difference`.
- `analytic`: explicitly unavailable and raises `NotImplementedError`.

## Cubic Dtype Convention

The physical adapter applies the established convention:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

The dense current assembly consumes KLM-normalized real cubic coefficients.
Silent complex-to-real casting remains disallowed.

## Finite-Grid Policy

Physical kernel assembly requires an explicit `KJSJIntegrationPolicy` for the
finite diagnostic \(z'\) sum. The policy is a non-production quadrature and
singularity prescription. It is not a production regulator or a regulator-
independence statement.

