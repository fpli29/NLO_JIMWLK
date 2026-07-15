# Physical KLM Kernel Integration Plan

This workflow adds non-production coordinate-space KLM kernel functions for
small dense diagnostics. It does not implement production evolution, physical
UV/IR regularization, score/Hessian-score training, or large-lattice
optimization.

## Goal

Implement diagnostic coordinate-space kernels for

\[
K_{JSJ},\quad
K_{JSSJ},\quad
K_{q\bar q},\quad
K_{JJSJ},\quad
K_{JJSSJ}.
\]

The immediate target is a small dense coordinate module that can feed
non-production skeleton metadata checks and later physical-kernel dipole
diagnostics.

## Source of Truth

`references/WORKNLO.tex` is the primary source. The PDF
`references/1405.0418v2.pdf` is only a cross-check. No formula is inferred if
the TeX source is ambiguous.

If a formula is not confidently available in a pointwise form usable by the
diagnostic interface, the corresponding function raises `NotImplementedError`
and the status file marks it pending.

## Kernel Types

The workflow distinguishes:

- unbarred singlet kernels used in KLM dipole-action formulas;
- barred/generalized kernels used for nonsinglet/full configuration-level
  Hamiltonian variants.

This first diagnostic implements unbarred singlet kernels where the pointwise
formula is explicit. Barred kernels are documented from the source but not
implemented here.

Names must keep this distinction explicit:

```python
KJJSJ_unbarred_value(...)
KJSSJ_unbarred_value(...)
```

No barred kernel is silently substituted for an unbarred one.

## Coordinate Conventions

Coordinates are two-dimensional transverse vectors:

```python
coords: np.ndarray  # shape (Nsite, 2)
```

The KLM notation from `WORKNLO.tex` lines 200--201 is:

\[
X=x-z,\qquad X'=x-z',\qquad Y=y-z,\qquad Y'=y-z',
\]

\[
W=w-z,\qquad W'=w-z'.
\]

The code uses integer site indices into `coords`, with helper functions for
differences, squared norms, dot products, and the two-dimensional scalar cross
product.

## Singularity Handling

KLM coordinate kernels contain singular denominators for coincident points and
for degenerate cross-ratio combinations. The diagnostic exposes an explicit
policy:

```python
singularity_policy = "raise" | "nan" | "eps"
```

- `"raise"` is the default and raises `KernelSingularityError` on exact zero
  denominators.
- `"nan"` returns `np.nan` in singular entries and is useful for dense array
  shape/nonfinite-count diagnostics.
- `"eps"` replaces exact zero denominators by an explicit positive `eps`
  argument. This is only a diagnostic regulator and is not physically final.

No regulator in this workflow is production-ready.

## Expected Symmetries

The source states the following symmetry information:

- `WORKNLO.tex` lines 256--258: \(K_{JSSJ}\) and \(K_{q\bar q}\) are symmetric
  under \(z\leftrightarrow z'\) and/or \(x\leftrightarrow y\).
- `WORKNLO.tex` lines 256--258 and 1304--1305: \(K_{JJSSJ}\) is antisymmetric
  under simultaneous
  \[
  x\leftrightarrow y,\qquad z\leftrightarrow z'.
  \]
- `WORKNLO.tex` line 1253: \(K_{JJSJ}(w,x,y;z)=-K_{JJSJ}(w,y,x;z)\).

Tests enforce the unambiguous pointwise symmetries for implemented kernels. If
a later formula requires integration or a barred nonsinglet modification, it
must be tested separately.

## Positivity Caveat

The existing non-production Pawula/positivity diagnostic remains a caveat. This
physical-kernel workflow does not claim physical NLO JIMWLK positivity or
non-positivity. Physical-kernel positivity checks are future work.
