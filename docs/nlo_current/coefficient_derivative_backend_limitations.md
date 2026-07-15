# Coefficient-Derivative Backend Limitations

The coefficient-derivative backend is finite-difference diagnostic
infrastructure for tiny dense lattices only.

## Scaling

- First derivatives require \(O(D)\) evaluations of full dense coefficient
  arrays.
- The nested \(K_3\) second derivative requires \(O(D^2)\) finite-difference
  derivative pairs.
- Dense \(K_3\) storage already scales as \(O(D^3)\).

This is acceptable for \(N_{\rm site}=1\) or similarly tiny diagnostic runs,
but it is not suitable for realistic lattice volumes.

## Numerical limitations

- Small finite-difference steps suffer from roundoff.
- Large finite-difference steps suffer from truncation error.
- Nested second derivatives are especially sensitive to the chosen
  `eps_second`.
- Product-rule tests use loose tolerances appropriate for finite differences,
  not exact symbolic identities.

## Scope limitations

- The backend differentiates coefficient callbacks after sector assembly.
- It does not implement production automatic differentiation.
- It does not exploit locality, sparse kernels, or translational structure.
- It does not train or evaluate score/Hessian-score models.

## Future production options

Possible production strategies remain open:

- analytic coefficient derivatives;
- automatic differentiation through coefficient builders;
- sparse or local kernel structure;
- stochastic trace or derivative estimators;
- combined validation against physical-kernel dipole evolution.
