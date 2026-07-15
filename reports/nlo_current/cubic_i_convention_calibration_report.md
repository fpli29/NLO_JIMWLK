# Cubic `i` Convention Calibration Report

Primary source: `references/WORKNLO.tex`.

This report diagnoses the one-`f` cubic convention only. It does not implement
`K_JJSJ` or `K_JJSSJ` Appendix A targets and does not mark either sector as
Appendix A passed.

## Calibration Setup

The K_JSJ comparison fixes the target convention used locally:

```text
target = TeX H_sector s
```

not `-H_sector s`.

The local dense generator convention is Hermitian:

```text
L F(U) = d F(exp(i eps t^a) U) / d eps
R F(U) = d F(U exp(i eps t^a)) / d eps
```

while `WORKNLO.tex` writes the algebraic generator action without the explicit
finite-difference `i` in the Wilson-line variation. For ordered cubic one-`f`
Hamiltonian words, the three Hermitian Lie derivatives and the existing direct
action sign convention combine to give the tested relation

```text
TeX target = (-i) * direct_action
```

for the one-`f` cubic sectors below.

Numerical setup:

- lattice size: `nsite = 3`
- seed: `66001`
- SU(3) Hermitian generators with `Tr(t^a t^b)=delta^{ab}/2`
- dipoles checked: `(0,1)`, `(1,2)`, `(0,2)`
- `K_JJSJ`: synthetic antisymmetric kernel
- `K_JJSSJ`: synthetic kernel antisymmetric under simultaneous
  `x <-> y`, `z <-> z'`

## Minimal One-`f` Operators

The tested representatives were:

| block | operator structure |
|---|---|
| K_JJSJ real | `f^{bde}[J_L^d J_L^e S_A^{ba} J_R^a - J_L^a S_A^{ab} J_R^d J_R^e]` |
| K_JJSJ virtual | `(1/3) f^{bde}[J_L^d J_L^e J_L^b - J_R^d J_R^e J_R^b]` |
| K_JJSSJ-like real | `f^{acb}[J_L^d J_L^e S_A^{dc} S_A^{eb} J_R^a - J_L^a S_A^{cd} S_A^{be} J_R^d J_R^e]` |
| K_JJSSJ-like virtual | `(1/3) f^{acb}[J_L^c J_L^b J_L^a - J_R^c J_R^b J_R^a]` |

The K_JJSSJ-like comparison used the isolated real and virtual formulas
transcribed from `WORKNLO.tex` lines 1308--1326 and 1336--1340. No full
K_JJSSJ implementation was added.

## K_JJSJ Checks

| dipole | real `direct/target` | virtual `direct/target` | real `|(-i)direct-target|` | virtual `|(-i)direct-target|` |
|---|---:|---:|---:|---:|
| `(0,1)` | `1.38e-15 + 1.00e+00 i` | `3.93e-15 + 1.00e+00 i` | `2.09e-15` | `1.98e-15` |
| `(1,2)` | `-2.47e-16 + 1.00e+00 i` | `3.33e-15 + 1.00e+00 i` | `1.28e-16` | `4.11e-16` |
| `(0,2)` | `0.00e+00 + 1.00e+00 i` | `4.28e-15 + 1.00e+00 i` | `6.28e-16` | `6.05e-15` |

Result: the real and virtual K_JJSJ isolated targets both satisfy
`TeX target = (-i) * direct_action` at roundoff for the tested small lattice.

## K_JJSSJ-Like Checks

| dipole | real `direct/target` | virtual `direct/target` | real `|(-i)direct-target|` | virtual `|(-i)direct-target|` |
|---|---:|---:|---:|---:|
| `(0,1)` | `1.19e-15 + 1.00e+00 i` | `-7.21e-15 + 1.00e+00 i` | `3.90e-16` | `4.43e-16` |
| `(1,2)` | `-1.16e-15 + 1.00e+00 i` | `-2.48e-15 + 1.00e+00 i` | `7.77e-16` | `6.17e-16` |
| `(0,2)` | `-5.48e-16 + 1.00e+00 i` | `1.59e-16 + 1.00e+00 i` | `3.33e-16` | `1.11e-16` |

Result: the same one-`f` cubic convention factor holds for the tested
K_JJSSJ-like real and virtual isolated blocks.

## Decision

The one-`f` cubic convention factor is calibrated for the tested ordered cubic
blocks:

```text
TeX target = (-i) * direct_action
```

This is sufficient to explain the K_JJSJ `i`-type mismatch observed in
`appendix_A_convention_diagnosis_report.md`.

No cubic target was implemented in this task. K_JJSJ and K_JJSSJ remain
pending because applying this convention factor to durable Appendix targets
requires a separate target-implementation pass, and K_JJSSJ still has
`tilde K` and combined-equation bookkeeping that should not be folded into this
K_qbarq task.

