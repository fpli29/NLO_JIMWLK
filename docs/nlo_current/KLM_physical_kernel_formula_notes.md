# KLM Physical Kernel Formula Notes

Primary source: `references/WORKNLO.tex`. Line ranges refer to the local TeX
source. These notes are for non-production dense coordinate diagnostics.

## Common Coordinate Notation

`WORKNLO.tex` lines 200--201 define:

\[
X\equiv x-z,\qquad X^\prime\equiv x-z^\prime,\qquad
Y\equiv y-z,\qquad Y^\prime\equiv y-z^\prime,
\]

\[
W\equiv w-z,\qquad W^\prime\equiv w-z^\prime.
\]

All squared quantities denote two-dimensional transverse norms.

## Kernel: `K_JSJ`

Formula:

```tex
K_{JSJ}(x,y;z) =&-&\frac{\alpha_s^2}{16 \pi^3}
\frac{(x-y)^2}{X^2 Y^2}\Big[b\ln(x-y)^2\mu^2
-b\frac{X^2-Y^2}{ (x-y)^2}\ln\frac{X^2}{Y^2}+
(\frac{67}{9}-\frac{\pi^2}{ 3})N_c-\frac{10}{ 9}n_f\Big]\nonumber \\
&-& \frac{N_c}{2}\ \int_{z^\prime}\, \tilde K(x,y,z,z^\prime).
```

WORKNLO.tex line range: 324--332.

Arguments: \(x,y,z\), plus scheme scale \(\mu\), \(N_c\), \(n_f\), and
\(\alpha_s\). The full formula also contains an integral over \(z'\).

Symmetries: Appendix A uses \(K_{JSJ}(u,v;z)=K_{JSJ}(v,u;z)\) and
\(K_{JSJ}(u,u;z)=0\), see lines 1117--1120.

Singular denominators: \(X^2\), \(Y^2\), \((x-y)^2\), plus singularities from
the integrated \(\tilde K\).

Barred/unbarred status: unbarred formula above; barred nonsinglet modification
is given in lines 367--370.

Implementation status: implemented as `KJSJ_unbarred_value(...)` only when an
explicit `KJSJIntegrationPolicy` is supplied. The default behavior raises rather
than silently choosing a quadrature measure or regulator.

Open questions: physical integration measure for \(z'\), UV/IR subtraction,
principal-value prescription, regulator independence, and production
scheme-scale handling.

## Kernel: `K_JJSSJ`

Formula:

```tex
K_{JJSSJ}(w;x,y;z,z^\prime)=-i
\frac{\alpha_s^2}{ 2\,\pi^4}
\left(\frac{X_iY^\prime_j}{ X^2Y^{\prime 2}}
\right)
\times \Big(\frac{\delta_{ij}}{2 (z-z^\prime)^2}
+\frac{(z^\prime-z)_i W^\prime_j}{ (z^\prime-z)^2 W^{\prime 2}}+
\frac{(z-z^\prime)_j W_i}{ (z-z^\prime)^2 W^{ 2}}
-\frac{W_i W^\prime_j}{ W^2 W^{\prime 2}}
\Big)\ln\frac{W^2}{ {W'}^2}
```

WORKNLO.tex line range: 288--297.

Arguments: \(w;x,y;z,z'\).

Symmetries: antisymmetric under simultaneous \(x\leftrightarrow y\) and
\(z\leftrightarrow z'\), stated in lines 256--258 and used explicitly in
lines 1304--1305.

Singular denominators: \(X^2\), \(Y'^2\), \(W^2\), \(W'^2\),
\((z-z')^2\), and the logarithm ratio \(W^2/W'^2\).

Barred/unbarred status: unbarred singlet kernel. Lines 385 state
\(K_{JJSJ}\) and \(K_{JJSSJ}\) remain unchanged in the nonsinglet
modification.

Implementation status: implemented as `KJJSSJ_unbarred_value`.

Open questions: physical singularity subtraction and integration policy.

## Kernel: `K_JJSJ`

Formula:

```tex
K_{JJSJ}(w;x,y;z)\,=\,-\,i\,\frac{\alpha_s^2}{ 4\, \pi^3 }\,
\Big[ \frac{X\cdot W}{ X^2\,W^2}\,-\,
\frac{Y\cdot W}{ Y^2\,W^2}    \Big]
\ln\frac{Y^2}{ (x-y)^2}\,\ln\frac{X^2}{ (x-y)^2},
```

WORKNLO.tex line range: 298--300.

Arguments: \(w;x,y;z\).

Symmetries: antisymmetric in \(x,y\), line 1253:
\(K_{JJSJ}(w,x,y;z)=-K_{JJSJ}(w,y,x;z)\).

Singular denominators: \(X^2\), \(Y^2\), \(W^2\), and \((x-y)^2\) in the
logarithms.

Barred/unbarred status: unbarred singlet kernel; lines 385 state it remains
unchanged in the nonsinglet modification.

Implementation status: implemented as `KJJSJ_unbarred_value`.

Open questions: physical singularity subtraction and integration policy.

## Kernel: `K_qbarq`

Formula:

```tex
K_{q\bar q}(x,y;z,z^\prime) =
-\frac{\alpha_s^2\,n_f}{ 8\,\pi^4}
\Big\{
\frac{{X'}^2Y^2+{Y'}^2X^2-(x-y)^2(z-z')^2}
{ (z-z')^4(X^2{Y'}^2-{X'}^2Y^2)}
\ln\frac{X^2{Y'}^2}{ {X'}^2Y^2}
-\frac{2}{(z-z^\prime)^4}\Big\}
```

WORKNLO.tex line range: 301--306.

Arguments: \(x,y;z,z'\).

Symmetries: \(K_{q\bar q}\) is symmetric under the appropriate
\(z\leftrightarrow z'\) and/or \(x\leftrightarrow y\) exchanges, lines
256--258.

Singular denominators: \((z-z')^4\), \(X^2Y'^2-X'^2Y^2\), and the logarithm
ratio.

Barred/unbarred status: unbarred formula above; barred nonsinglet modification
is given in lines 378--383 using \(I_f\).

Implementation status: implemented as `Kqbarq_unbarred_value`.

Open questions: barred nonsinglet extension and physical singularity
subtraction.

## Kernel: `tilde_K`

Formula:

```tex
\tilde K(x,y,z,z^\prime)\,=\frac{i}{2}\,
\Big[K_{JJSSJ}(x;x,y;z,z^\prime)
-K_{JJSSJ}(y;x,y;z,z^\prime)
-K_{JJSSJ}(x;y,x;z,z^\prime)
+K_{JJSSJ}(y;y,x;z,z^\prime)\Big]\ ,
```

WORKNLO.tex line range: 307--311.

Arguments: \(x,y,z,z'\).

Symmetries: Appendix A comments record a symmetric dipole combination in lines
1331--1333; pointwise helper is tested directly against the definition.

Singular denominators: inherited from `K_JJSSJ`.

Barred/unbarred status: unbarred singlet combination.

Implementation status: implemented as `tilde_K_JJSSJ_unbarred_value`.

Open questions: physical integration/subtraction when \(\tilde K\) is
integrated into `K_JSJ`.

## Kernel: `K_JSSJ`

Formula:

```tex
K_{JSSJ}(x,y;z,z^\prime) = \frac{\alpha_s^2}{16\,\pi^4}
\Bigg[\,-\,\frac{4}{ (z-z^\prime)^4}\,+\,
\Big\{2\frac{X^2{Y'}^2+{X'}^2Y^2-4(x-y)^2(z-z')^2}
{ (z-z')^4[X^2{Y'}^2-{X'}^2Y^2]}
+
~\frac{(x-y)^4}{ X^2{Y'}^2-{X'}^2Y^2}\Big[
\frac{1}{X^2{Y'}^2}+\frac{1}{ Y^2{X'}^2}\Big]
+\frac{(x-y)^2}{(z-z')^2}\Big[
\frac{1}{X^2{Y'}^2}-\frac{1}{ {X'}^2Y^2}\Big]\Big\}
\ln\frac{X^2{Y'}^2}{ {X'}^2Y^2}\Bigg]
+\,\tilde K(x,y,z,z^\prime)
```

WORKNLO.tex line range: 313--323.

Arguments: \(x,y;z,z'\), with `tilde_K` from lines 307--311.

Symmetries: \(K_{JSSJ}\) is symmetric under the appropriate
\(z\leftrightarrow z'\) and/or \(x\leftrightarrow y\) exchanges, lines
256--258.

Singular denominators: \((z-z')^4\), \(X^2Y'^2-X'^2Y^2\), \(X^2Y'^2\),
\(Y^2X'^2\), \(X'^2Y^2\), and inherited `tilde_K` denominators.

Barred/unbarred status: unbarred formula above; barred nonsinglet modification
is given in lines 374--383 using \(I\).

Implementation status: implemented as `KJSSJ_unbarred_value`.

Open questions: barred nonsinglet extension and physical singularity
subtraction.
