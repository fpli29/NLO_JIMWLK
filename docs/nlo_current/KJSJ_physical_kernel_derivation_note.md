# Diagnostic K_JSJ Physical Kernel Derivation Note

Primary source: `references/WORKNLO.tex` lines 324--332.

This note transcribes the unbarred singlet \(K_{JSJ}\) formula and records the
finite-grid diagnostic policy used by the non-production implementation. It
does not define a production regulator or a physical quadrature prescription.

## Source Formula

`WORKNLO.tex` lines 324--332 give:

\[
K_{JSJ}(x,y;z)
=
-\frac{\alpha_s^2}{16\pi^3}
\frac{(x-y)^2}{X^2Y^2}
\left[
b\ln (x-y)^2\mu^2
-b\frac{X^2-Y^2}{(x-y)^2}\ln\frac{X^2}{Y^2}
+\left(\frac{67}{9}-\frac{\pi^2}{3}\right)N_c
-\frac{10}{9}n_f
\right]
-\frac{N_c}{2}\int_{z'}\widetilde K(x,y,z,z').
\]

The same lines define

\[
b=\frac{11}{3}N_c-\frac{2}{3}n_f,
\]

and state that \(\mu\) is the \(\overline{\rm MS}\) normalization point.

The coordinate notation from `WORKNLO.tex` lines 200--201 is:

\[
X=x-z,\qquad Y=y-z.
\]

The \(\widetilde K\) combination is defined in lines 307--311:

\[
\widetilde K(x,y,z,z')
=
\frac{i}{2}\left[
K_{JJSSJ}(x;x,y;z,z')
-K_{JJSSJ}(y;x,y;z,z')
-K_{JJSSJ}(x;y,x;z,z')
+K_{JJSSJ}(y;y,x;z,z')
\right].
\]

## Decomposition

The implementation decomposes \(K_{JSJ}\) into:

1. Local coordinate-dependent terms:

\[
-\frac{\alpha_s^2}{16\pi^3}
\frac{(x-y)^2}{X^2Y^2}
\left[
-b\frac{X^2-Y^2}{(x-y)^2}\ln\frac{X^2}{Y^2}
+\left(\frac{67}{9}-\frac{\pi^2}{3}\right)N_c
-\frac{10}{9}n_f
\right].
\]

2. Logarithmic scheme-scale term:

\[
-\frac{\alpha_s^2}{16\pi^3}
\frac{(x-y)^2}{X^2Y^2}
b\ln (x-y)^2\mu^2.
\]

3. Explicit diagnostic finite-sum version of the source integral:

\[
-\frac{N_c}{2}\int_{z'}\widetilde K(x,y,z,z')
\quad\rightarrow\quad
-\frac{N_c}{2}\sum_{z'\in{\cal Q}} w_{z'}\widetilde K(x,y,z,z').
\]

4. Diagonal-zero subtraction:

The source discussion states that kernels were assigned so they vanish at
\(x=y\), analogously to the LO dipole kernel (`WORKNLO.tex` lines 344--346),
and Appendix A uses \(K_{JSJ}(u,u;z)=0\) (`WORKNLO.tex` lines 1117--1120).
The diagnostic implementation therefore returns zero for \(x=y\) before
evaluating singular logarithms.

## Diagnostic Policy Interface

The finite sum is controlled by `KJSJIntegrationPolicy`:

```python
KJSJIntegrationPolicy(
    quadrature_weights=...,
    mu=...,
    excluded_indices=...,
    exclude_coincident_labels=...,
    principal_value="none",
    subtraction="diagonal_zero",
    finite_volume_boundary="finite_coordinate_sum",
    description="...",
)
```

The policy explicitly records:

- \(z'\) quadrature weights;
- excluded coincident points, e.g. `("x", "y", "z")`;
- whether any principal-value handling is used;
- subtraction handling;
- finite-volume boundary label;
- \(\mu\) and the scheme constants entering the formula.

Only `principal_value="none"` is implemented. Any physical principal-value or
subtraction prescription remains future work.

## Ambiguity and Scope

The KLM formula fixes the continuum expression but does not uniquely determine
a finite-grid quadrature, coincident-point exclusion set, or production
regulator. The code therefore refuses to evaluate `KJSJ_unbarred_value(...)`
unless an explicit `KJSJIntegrationPolicy` is supplied.

This resolves the diagnostic implementation without claiming regulator
independence, physical positivity, or production readiness.
