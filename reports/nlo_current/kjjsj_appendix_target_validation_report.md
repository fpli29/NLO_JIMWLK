# K_JJSJ Appendix Target Validation Report

Primary formula source: `references/WORKNLO.tex`.

This report validates only the `K_JJSJ` Appendix A target. It does not
implement `K_JJSSJ`; `K_JSSJ` is handled by its separate validation report.

## Source Formulas

Real contribution:

- `WORKNLO.tex` lines 1208--1224
- compact formula body lines 1216--1223
- uses the antisymmetry stated on line 1253:
  `K_JJSJ(w,x,y;z) = -K_JJSJ(w,y,x;z)`

Implemented formula:

```text
(i/2) int_z [
  K(v;u,v;z) - K(v;v,u;z) - K(u;u,v;z) + K(u;v,u;z)
] [s(u,v) - Nc^2 s(u,z)s(z,v)]
```

Virtual contribution:

- `WORKNLO.tex` lines 1244--1253
- formula body lines 1246--1251

Implemented formula:

```text
i (Nc^2 - 1)/3 int_z [K(u,v,u,z) + K(v,u,v,z)] s(u,v)
```

## Convention

The K_JSJ sector fixes the local target convention as:

```text
Appendix target = H_sector s
```

The cubic one-`f` convention calibration in
`reports/nlo_current/cubic_i_convention_calibration_report.md` fixes the
comparison of current Hermitian-generator direct actions as:

```text
Appendix cubic target = (-i) * raw direct cubic action
```

The implementation therefore adds the explicit helper:

```text
klm_normalized_cubic_direct_action(raw_direct) = (-i) * raw_direct
```

Existing direct-action functions were not redefined.

## Dense Validation

Configuration:

- seed: `67001`
- `nsite = 3`
- dipoles: `(0,1)`, `(1,2)`, `(0,2)`
- kernel: synthetic antisymmetric `K_JJSJ[w,x,y,z]`

Residuals compare the implemented Appendix targets against
`klm_normalized_cubic_direct_action(action_KJJSJ_direct(...))`.

| dipole | real residual | virtual residual | full residual | raw uncalibrated full residual |
|---|---:|---:|---:|---:|
| `(0,1)` | `1.3608726004012153e-15` | `3.5596458096434965e-15` | `2.2211398302543545e-15` | `9.1773691740376000e-01` |
| `(1,2)` | `9.777686967753793e-16` | `2.8609792490763985e-16` | `6.855121490028353e-16` | `2.8327374324235194e-01` |
| `(0,2)` | `4.322212449724815e-16` | `1.2600542234970832e-15` | `8.2826097710246705e-16` | `1.0098287391359802e-01` |

Summary:

- max real residual: `1.3608726004012153e-15`
- max virtual residual: `3.5596458096434965e-15`
- max full residual: `2.2211398302543545e-15`
- min raw uncalibrated full residual: `1.0098287391359802e-01`

The raw direct-action failure confirms that the cubic `i`-convention
calibration is numerically meaningful.

## Decision

`K_JJSJ` is a locally implemented and tested Appendix A target, with
real/virtual separation preserved. `K_JJSSJ` remains pending.
