# K_JJSSJ Appendix Target Validation Report

Primary formula source: `references/WORKNLO.tex`.

This pass implements the full isolated `K_JJSSJ` Appendix A target: real plus
virtual contributions, including the pure eight-kernel combination and the
`\widetilde K` terms.

## Source Formulas

Tilde-K definition:

- `WORKNLO.tex` lines 307--311

```text
tilde K(x,y,z,z') = (i/2) [
  K(x; x,y; z,z') - K(y; x,y; z,z')
  - K(x; y,x; z,z') + K(y; y,x; z,z')
]
```

Real contribution:

- `WORKNLO.tex` lines 1304--1328
- final formula lines 1322--1326

```text
(i/2Nc) int_{z,z'} [
  K(u;u,u) - K(u;v,u) + K(u;v,v) - K(u;u,v)
  + K(v;u,v) - K(v;u,u) + K(v;v,u) - K(v;v,v)
] Nc^3 s(z,v)s(z',z)s(u,z')

+ (1/Nc) int_{z,z'} tildeK(u,v,z,z') [
  Nc^3 s(z,v)s(z',z)s(u,z')
  - tr(S(v)Sdag(z)S(z')Sdag(u)S(z)Sdag(z'))
]
```

Virtual contribution:

- `WORKNLO.tex` lines 1334--1341

```text
-((Nc^2 - 1)/3) int_{z,z'} tildeK(u,v,z,z') s(u,v)
```

## Convention

The K_JSJ sector fixes the local target convention as:

```text
Appendix target = H_sector s
```

The cubic one-`f` convention calibration fixes comparison of the current
Hermitian-generator direct action as:

```text
Appendix cubic target = (-i) * raw direct cubic action
```

Existing direct-action functions were not redefined.

## Dense Validation

Configuration:

- seed: `70001`
- `nsite = 3`
- dipoles: `(0,1)`, `(1,2)`, `(0,2)`
- kernel: synthetic KLM-antisymmetric `K_JJSSJ[w,x,y,z,z']`

| dipole | real residual | virtual residual | full residual | raw uncalibrated full residual | no-`1/3` virtual change |
|---|---:|---:|---:|---:|---:|
| `(0,1)` | `2.6990424254992210e-16` | `3.3489362443588235e-16` | `9.0205620750793970e-17` | `3.1356671734340913e-01` | `9.3806988396103060e-03` |
| `(1,2)` | `3.5529186548606245e-15` | `4.8321064166698255e-15` | `1.2827670180810443e-15` | `4.8499571852332940e-01` | `8.0453761306219870e-01` |
| `(0,2)` | `5.5250329970541970e-16` | `1.0111756560930368e-15` | `5.3604920268812600e-16` | `1.6856446621016427e-01` | `2.4358552288881750e-01` |

Summary:

- max real residual: `3.5529186548606245e-15`
- max virtual residual: `4.8321064166698255e-15`
- max full residual: `1.2827670180810443e-15`
- min raw uncalibrated full residual: `1.6856446621016427e-01`
- min no-`1/3` virtual change: `9.3806988396103060e-03`
- max `tildeK` magnitude: `1.8216929024694406e-01`

Contribution exercise:

| dipole | pure eight-kernel real norm | tilde-K real norm |
|---|---:|---:|
| `(0,1)` | `4.7703097666632770e-02` | `2.7330244211019095e-01` |
| `(1,2)` | `4.7887288802575195e-01` | `2.7651396087213270e-01` |
| `(0,2)` | `2.7284804228652476e-01` | `4.5043614137456180e-01` |

The raw direct-action failure confirms that the cubic `i` calibration remains
necessary. The no-`1/3` virtual check confirms the virtual normalization is
active. Both the pure eight-kernel contribution and the tilde-K contribution
are nonzero for the tested synthetic kernel.

## Decision

`K_JJSSJ` is a locally implemented and tested Appendix A target. With this
pass, all five Appendix A dipole sector targets in the non-production dense
validation are implemented and pass.

