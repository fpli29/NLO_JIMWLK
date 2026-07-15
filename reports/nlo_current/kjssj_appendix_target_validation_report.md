# K_JSSJ Appendix Target Validation Report

Primary formula source: `references/WORKNLO.tex`.

This pass implements only the full isolated `K_JSSJ` Appendix A target. It does
not implement `K_JJSSJ`.

## Source Formula

The subsection expression in `WORKNLO.tex` lines 1158--1164 covers only the
`f f J_L S S J_R` component. It is not the full `K_JSSJ` sector target because
the Hamiltonian also contains the `-N_c J_L S_A J_R` subtraction.

The full isolated `K_JSSJ` target is taken from the `K_JSSJ` part of the
combined equation in `WORKNLO.tex` lines 1353--1355:

```text
-(1/Nc) int_{z,z'} K_JSSJ(u,v;z,z') [
  Nc^3 s(u,z') s(z',z) s(z,v)
  - tr(S(v) S^\dagger(z) S(z') S^\dagger(u) S(z) S^\dagger(z'))
  - Nc^3 s(u,z) s(z,v)
  + Nc s(u,v)
]
```

The `\widetilde K` part of the combined factor
`K_JSSJ(u,v;z,z') - \widetilde K(u,v,z,z')` is excluded from this isolated
target because it belongs to the `K_JJSSJ`-related combination.

## Kernel Convention

The simplified combined-equation form is validated with the Appendix endpoint
condition used in the TeX simplifications:

- symmetric in endpoint coordinates `(x,y)`
- zero endpoint diagonal `K_JSSJ(x,x;z,z') = 0`
- symmetric in `(z,z')`

Under these conditions, the isolated combined-equation target matches the
direct dense `K_JSSJ` action.

## Dense Validation

Configuration:

- seed: `69001`
- `nsite = 3`
- dipoles: `(0,1)`, `(1,2)`, `(0,2)`
- kernel: synthetic Appendix-compatible `K_JSSJ`

| dipole | full target residual | subsection-as-full residual | fake-tilde contamination delta |
|---|---:|---:|---:|
| `(0,1)` | `1.4571677198205180e-16` | `2.1503160210840777e-01` | `7.1336076700438900e-01` |
| `(1,2)` | `5.0250058654113200e-17` | `3.3001781996773796e-02` | `3.1995236396439400e+00` |
| `(0,2)` | `1.4752290795525882e-16` | `4.3815424383340174e-01` | `2.1751463593252837e-01` |

Summary:

- max full target residual: `1.4752290795525882e-16`
- min subsection-as-full residual: `3.3001781996773796e-02`
- min fake-tilde contamination delta: `2.1751463593252837e-01`

The subsection formula fails as a full-sector target for the tested generic
Appendix-compatible kernel. The fake-tilde check verifies that adding a
`\widetilde K` contribution would change the target, so the implemented
isolated target contains no tilde-K contamination.

## Decision

`K_JSSJ` is now a locally implemented and tested full Appendix A target.
`K_JJSSJ` remains pending.

