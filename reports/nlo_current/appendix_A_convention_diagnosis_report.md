# Appendix A Convention Diagnosis Report

Primary source: `references/WORKNLO.tex`.

This pass did not modify production evolution code and did not update
`src/nlo_current/dipole_appendix_targets.py`. No new sector is marked as an
Appendix A pass.

## Source Formula Ranges

The transcription notes in
`docs/nlo_current/KLM_appendix_A_dipole_targets_notes.md` now record exact TeX
line ranges for:

| formula | `WORKNLO.tex` source lines |
|---|---:|
| `tilde K` definition | 307--311; body 309--310 |
| `K_JSJ` formula | 1115--1120; body 1116--1118 |
| `K_JSSJ` subsection formula | 1156--1166; body 1158--1164 |
| `K_qbarq` subsection formula | 1172--1181; body 1174--1179 |
| `K_JJSJ` real formula | 1208--1224; compact body 1216--1223 |
| `K_JJSJ` virtual formula | 1244--1253; body 1246--1251 |
| `K_JJSSJ` real formula | 1304--1328; final body 1308--1326 |
| `K_JJSSJ` virtual formula | 1334--1341; body 1336--1340 |
| combined `ourdipole` equation | 1344--1360; body 1348--1359 |

## Numerical Setup

Diagnostics used the existing dense small-lattice utilities only:

- lattice size: `nsite = 3`
- color group: existing SU(3) Hermitian generators, `Tr(t^a t^b)=delta^{ab}/2`
- seed: `64001`
- kernels: synthetic endpoint-compatible `K_JSJ`, endpoint-compatible
  symmetric `K_qbarq`, and antisymmetric `K_JJSJ`
- dipoles checked: `(u,v) = (0,1), (1,2), (0,2)`

## Direct-Action Convention

The direct-action functions cannot currently be classified by one global
`H_sector s` versus `-H_sector s` sign.

Calibration anchor:

| check | max residual |
|---|---:|
| `K_JSJ`: `direct - Appendix target` | `7.85e-17` |
| `K_JSJ`: `direct + Appendix target` | `1.09e+00` |

Therefore `action_KJSJ_direct` is calibrated to the TeX Appendix A
`H_KJSJ s` convention, not to `-H_KJSJ s`, despite the local docstring wording.

Sector diagnosis:

| sector | diagnosis |
|---|---|
| `K_JSJ` | returns TeX `H_sector s` under the stated Appendix kernel symmetry and zero-diagonal condition |
| `K_qbarq` | trace-current component returns the subsection target; full direct action contains the extra subtraction term and should not be compared to the subsection alone |
| `K_JJSJ` | current direct action equals `i` times the transcribed TeX real+virtual target on the tested small lattices |
| `K_JJSSJ` | not convention-fit in this pass; classified only because of `tilde K` and combined-equation complications |

Conclusion: the current direct-action set is sector-dependent under the
implementation conventions. The `K_JSJ` anchor proves there is no global
`Hs` versus `-Hs` sign flip.

## `K_qbarq` Diagnosis

The `WORKNLO.tex` subsection at lines 1174--1179 isolates only

```tex
2 J_L^a(x) tr[S^\dagger(z)t^a S(z')t^b] J_R^b(y)
```

The existing direct full-sector builder uses

```text
2 Tr(U_z^\dagger t^a U_z' t^b) - S_A^{ab}(z)
```

inside the ordered `J_L A J_R` block.

Component split:

| check | max residual |
|---|---:|
| `direct_full - (trace_current + subtraction)` | `1.39e-16` |
| local `J_L J_R` order versus reversed `J_R J_L` order | `0.00e+00` |

Comparison to the subsection target, using `S = U` as in the local dipole
observable:

| dipole | `|full - target|` | `|-full - target|` | `|trace_current - target|` |
|---|---:|---:|---:|
| `(0,1)` | `3.70e-02` | `5.53e-01` | `1.14e-16` |
| `(1,2)` | `8.53e-02` | `3.50e-01` | `5.72e-17` |
| `(0,2)` | `4.45e-03` | `2.50e-01` | `4.65e-17` |

Trace orientation check:

| dipole | `|full - target(S=U^\dagger)|` | `|-full - target(S=U^\dagger)|` |
|---|---:|---:|
| `(0,1)` | `7.27e-02` | `5.50e-01` |
| `(1,2)` | `2.20e-01` | `2.85e-01` |
| `(0,2)` | `8.88e-02` | `2.34e-01` |

Cause tests:

| candidate cause | result |
|---|---|
| missing `-J_L S_A J_R` subtraction in subsection target | yes, if the subsection target is used as a full-sector target; the trace-current component alone matches the subsection at roundoff |
| Hermitian versus anti-Hermitian generator convention | not supported for this two-generator comparison; a two-generator convention flip would act as a sign flip, which fails, and `K_JSJ` calibrates the two-generator sign |
| trace orientation `U` versus TeX `S` notation | not the cause; `S=U` makes the trace-current component match exactly, while `S=U^\dagger` is worse |
| `Hs` versus `-Hs` sign | not the cause; `-full` is farther from the subsection target |
| left/right generator ordering | not the cause; the local LR/RL order check is zero within arithmetic precision |

Conclusion: `K_qbarq` subsection target is partial unless the subtraction is
included. The combined equation target is available in `WORKNLO.tex` line 1352
and should be used for a full-sector Appendix comparison.

## `K_JJSJ` Diagnosis

The real and virtual isolated targets are available in `WORKNLO.tex` lines
1208--1224 and 1244--1253. With the transcribed real+virtual target:

| dipole | `direct / target` | `|direct - target|` | `|-direct - target|` | `|(-i) direct - target|` |
|---|---:|---:|---:|---:|
| `(0,1)` | `2.57e-15 + 1.00e+00 i` | `7.43e-01` | `7.43e-01` | `2.01e-15` |
| `(1,2)` | `5.18e-16 + 1.00e+00 i` | `2.06e-01` | `2.06e-01` | `4.75e-16` |
| `(0,2)` | `5.53e-16 + 1.00e+00 i` | `1.40e-01` | `1.40e-01` | `1.39e-16` |

The same `(-i)` relation holds separately for the real and virtual pieces:

| piece | max `|(-i) direct_piece - target_piece|` |
|---|---:|
| real | `1.84e-15` |
| virtual | `6.96e-16` |

Cause tests:

| candidate cause | result |
|---|---|
| structure-constant convention `f^{abc}` sign | not supported; flipping `f` flips the direct action and has the same failing residual as `-direct` |
| Hermitian generator convention | possible but not uniquely identified; the observed error is a uniform `i` factor in a cubic sector, which is exactly where generator/i conventions matter |
| missing factor `i` in direct action | numerically consistent; multiplying the current direct result by `-i` matches the transcribed TeX target at roundoff |
| `Hs` versus `-Hs` sign | not supported; both `direct` and `-direct` have comparable large residuals |
| kernel argument ordering | not supported by the tested alternate swapped-argument target; the mismatch is a uniform `i`, not an argument-dependent residual |

Conclusion: the `K_JJSJ` mismatch is an `i`-type convention ambiguity. The
data rule out `f` sign, global Hamiltonian sign, and tested kernel-argument
ordering as the cause. They do not by themselves distinguish a missing `-i`
factor in the direct cubic action from a Hermitian/anti-Hermitian convention
translation. No conversion factor should be added until this convention is
fixed from an independent source.

## `K_JJSSJ` Classification Only

No full implementation or full convention fit was attempted.

| item | classification |
|---|---|
| real isolated target | available at lines 1308--1326; contains a pure eight-kernel piece and a `tilde K` trace/triple-dipole piece |
| virtual isolated target | available at lines 1336--1340; collapses entirely to a `tilde K` contribution |
| `tilde K` terms | definition lines 307--311; appears in the real isolated target line 1326, virtual target line 1340, and combined equation lines 1350 and 1353 |
| pure eight-kernel combined term | available in `ourdipole` lines 1356--1358 |

## Target Availability Summary

| sector | status |
|---|---|
| `K_JSJ` | full target, implemented |
| `K_JSSJ` | partial subsection; full combined-equation target required |
| `K_qbarq` | subsection target is partial unless subtraction is included; combined-equation target available |
| `K_JJSJ` | real + virtual isolated target available, convention diagnosis pending |
| `K_JJSSJ` | real + virtual target available with `tilde K` complications, convention diagnosis pending |

