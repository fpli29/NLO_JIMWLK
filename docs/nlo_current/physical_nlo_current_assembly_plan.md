# Physical NLO Current Assembly Plan

## Scope

This document describes the first end-to-end non-production physical NLO
generalized current assembly. It is dense, small-lattice, and diagnostic only.
It does not implement production evolution, train score/Hessian-score models,
optimize large lattices, claim regulator independence, or claim physical
positivity.

## Data Flow

```text
physical unbarred KLM kernels
      |
      v
physical_kernel_adapter
      |
      v
KLM-normalized dense kernel arrays
      |
      v
assemble_nlo_current_terms(...)
      |
      v
K1, K2, K3 coefficients
      |
      v
generalized current / product-rule contractions
      |
      v
velocity v^A
```

The implemented API is `src/nlo_current/physical_nlo_current.py`.

## Kernel Conventions

The adapter uses the already validated physical kernel layer:

- `KJSJ`, `KJSSJ`, and `Kqbarq` are real unbarred singlet diagnostic kernels.
- raw physical `KJJSJ` and `KJJSSJ` carry explicit `-i` factors from WORKNLO.
- the physical adapter maps raw cubic kernels to KLM-normalized coefficients:

```text
KLM-normalized cubic coefficient = (-1j) * raw physical cubic kernel
```

The dense normal-form assembly uses the KLM-normalized real cubic coefficients.
Raw complex cubic tensors remain available only as explicit diagnostics; no
silent complex-to-real cast is allowed.

## Coefficient Assembly

`assemble_physical_terms(...)` builds all dense coefficient tensors:

- `K1`: lower-order commutator drift from validated cubic sectors.
- `K2`: second-order ordered current coefficients from `KJSJ`, `KJSSJ`, and
  `Kqbarq`, plus cubic commutator corrections when enabled.
- `K3`: third-order ordered current coefficients from `KJJSJ` and `KJJSSJ`.

Convenience accessors:

- `assemble_physical_K1(...)`
- `assemble_physical_K2(...)`
- `assemble_physical_K3(...)`

All outputs carry metadata with kernel origin, sector labels, dtype, maximum
imaginary component, and the cubic normalization convention.

## Derivative Backends

The velocity evaluator requires product-rule derivative contractions:

```text
L_B K2^{A B}
L_B L_C K3^{A B C}
```

The wrapper exposes `compute_physical_coefficient_derivatives(...)` with an
explicit backend selector:

- `finite_difference`: dense central finite differences on tiny lattices.
- `diagnostic`: alias for `finite_difference`.
- `analytic`: reserved for future closed-form physical derivatives; currently
  raises `NotImplementedError` because no analytic physical-kernel derivative
  implementation is available in this layer.

The finite-difference backend uses the existing diagnostic coefficient
derivative implementation and recomputes dense `K2` and `K3` callbacks under
left perturbations of the Wilson lines.

## Velocity Evaluation

`evaluate_physical_nlo_velocity(...)` computes:

```text
v^A = K1^A
      - 1/2 [L_B K2^{A B} + K2^{A B} s_B]
      + 1/6 [L_B L_C K3^{A B C}
             + L_C K3^{A B C} s_B
             + L_B K3^{A B C} s_C
             + K3^{A B C}(H_BC + s_B s_C)]
```

The score `s_B` and Hessian-score `H_BC = L_B s_C` are supplied inputs. This
layer does not train or estimate them.

## Diagnostic Controls

The assembly requires explicit physical-kernel policies:

- `KJSJIntegrationPolicy` for the \(z'\) integral in `KJSJ`.
- explicit `singularity_policy` and optional `eps` for finite dense arrays.
- optional sector filters for sector-by-sector checks and LO/K2-only tests.

The finite-grid and singularity policies are diagnostic choices only.

## Limitations

- No barred/nonsinglet kernels.
- No production UV/IR regulator or quadrature prescription.
- No score or Hessian-score model training.
- No physical positivity claim; Pawula caveats remain.
- No large-lattice sparse implementation.
- Analytic physical coefficient derivatives remain future work.
