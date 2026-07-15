# K_JSJ Physical Kernel Validation Report

## Scope

This report validates the non-production diagnostic implementation of the
unbarred physical \(K_{JSJ}\) kernel. It does not define a production regulator,
claim regulator independence, or claim physical positivity.

## Source

- primary formula: `references/WORKNLO.tex` lines 324--332
- tilde-K definition: `references/WORKNLO.tex` lines 307--311
- derivation note: `docs/nlo_current/KJSJ_physical_kernel_derivation_note.md`

## Implementation Decomposition

The implementation separates:

- `KJSJ_unbarred_local_value(...)`: local coordinate-dependent and
  \(\mu\)-dependent logarithmic terms;
- `KJSJ_unbarred_tilde_integral_value(...)`: explicit diagnostic finite sum
  for \(\int_{z'}\widetilde K(x,y,z,z')\);
- `KJSJ_unbarred_value(...)`: full diagnostic value
  \(K_{\rm local}-N_c\int\widetilde K/2\);
- diagonal-zero handling for \(x=y\), as used by KLM for singlet kernels.

The default behavior raises unless a `KJSJIntegrationPolicy` is supplied.

## Diagnostic Policy Used

Representative policy:

```text
KJSJIntegrationPolicy(
    quadrature_weights = equal finite-coordinate weights,
    mu = 1.3,
    exclude_coincident_labels = ("x", "y", "z"),
    principal_value = "none",
    subtraction = "diagonal_zero",
    finite_volume_boundary = "finite_coordinate_sum",
)
```

This is a finite-grid diagnostic policy only.

## Pointwise Checks

Using coordinates
`[[0.0, 0.0], [1.0, 0.2], [0.3, 1.1], [1.4, 1.3], [2.1, 0.7]]`,
`Nc=3`, `nf=2`, `alpha_s=0.3`, `mu=1.3`, and `eps=1e-6`:

- full value: `-1.7522658138122986e-03`
- local value: `-1.7519671644370793e-03`
- tilde integral value: `1.9909958347957940e-07`
- \(x\leftrightarrow y\) symmetry residual: `0.0`
- diagonal \(x=y\) value: `0.0`
- tilde-integral linearity residual under doubled quadrature weights: `0.0`
- \(\mu=1\) versus \(\mu=2\) sensitivity: `1.4960694107815898e-03`

The local sign check is consistent with the negative prefactor in
`WORKNLO.tex` line 325 for the tested positive-bracket configuration.

## Quadrature Refinement Diagnostic

For a deterministic ring quadrature around a fixed \((x,y,z)\) triple:

| ring points | K_JSJ value |
|---:|---:|
| 8 | `-3.8152465169998787e-04` |
| 16 | `-3.5631009142516115e-05` |
| 32 | `-2.0320832475300365e-04` |

The refinement deltas were:

- `|K_16 - K_8| = 3.4589364255747175e-04`
- `|K_32 - K_16| = 1.6757731561048753e-04`

This is only a small-grid convergence diagnostic. It is not a continuum
quadrature proof.

## Physical-Kernel Dipole Recheck

`scripts/nlo_current/full_dipole_validation_physical_kernels.py` was run with
an explicit diagnostic finite-grid policy.

| sector | status | residual | relative residual |
|---|---|---:|---:|
| \(K_{JSJ}\) | passed | `7.7688469830766467e-12` | `2.8320504613732165e-09` |
| \(K_{JSSJ}\) | passed | `2.1438317028805034e-13` | `1.7081190426370902e-11` |
| \(K_{q\bar q}\) | passed | `3.6313016435989580e-13` | `1.9663547428238200e-10` |
| \(K_{JJSJ}\) | passed | `6.1080414466655290e-20` | `1.4244164473524333e-15` |
| \(K_{JJSSJ}\) | passed | `1.1408362482102489e-17` | `2.2558036620974433e-15` |

The \(K_{q\bar q}\) finite-grid mismatch was resolved by using the exact
WORKNLO trace-product expression from lines 1174--1177 as the diagnostic
target. No production quadrature or regulator is fixed by this check.

## Conclusion

The diagnostic \(K_{JSJ}\) implementation satisfies the requested local
checks, policy requirements, symmetry, diagonal-zero condition, tilde-linearity,
\(\mu\) sensitivity, and small-grid quadrature refinement check.

The physical-kernel dipole recheck is fully passed under the explicit
non-production diagnostic policies used here. This report does not claim
physical NLO JIMWLK positivity or non-positivity.
