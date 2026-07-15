# KLM Appendix A Dipole Target Notes

Primary source: `references/WORKNLO.tex`.

PDF cross-check: `references/1405.0418v2.pdf`; `pdftotext` was unavailable, so
`PyPDF2` text extraction was used only to confirm the Appendix A heading and
the presence of the corresponding equations around PDF pages 26--32. The TeX
source remains the source of record.

## Scope

The Appendix A observable is the quark dipole \(s(u,v)\). The TeX states that
the Hamiltonian action is obtained by applying \(H^{NLO\,JIMWLK}\) to
\(s(u,v)\), and then lists the action of each kernel-labeled term.

These notes preserve the source signs, \(N_c\) factors, \(i\) factors, trace
ordering, kernel arguments, and exchange notation. They do not infer missing
formulas.

## \(\widetilde K\) Definition

Source range: `references/WORKNLO.tex` lines 307--311; formula body lines
309--310.

The Appendix A \(K_{JJSSJ}\) formulas refer to \(\widetilde K\), defined earlier
in the TeX as

```tex
\tilde K(x,y,z,z^\prime)\,=\frac{i}{2}\,\Big[K_{JJSSJ}(x;x,y;z,z^\prime)&-&K_{JJSSJ}(y;x,y;z,z^\prime)-K_{JJSSJ}(x;y,x;z,z^\prime)\nonumber \\
&+&K_{JJSSJ}(y;y,x;z,z^\prime)\Big]\ ,
```

## \(K_{JSJ}\)

Source range: `references/WORKNLO.tex` lines 1115--1120; formula body lines
1116--1118.

Classification: isolated sector target after applying the stated kernel
conditions \(K_{JSJ}(u,v;z)=K_{JSJ}(v,u;z)\) and \(K_{JSJ}(u,u;z)=0\).

```tex
\begin{equation}\begin{split}&\int_{x,y,z}K_{JSJ}(x,y;z)\left[J_{L}^{a}(x)J_{L}^{a}(y)+J_{R}^{a}(x)J_{R}^{a}(y)-2J_{L}^{a}(x)S_A^{ab}(z)J_{R}^{b}(y)\right]s(u,v)\\
&=N_{c}\int_{z}\left[K_{JSJ}(v,v;z)+K_{JSJ}(u,u;z)-K_{JSJ}(u,v;z)-K_{JSJ}(v,u;z)\right]\\
&\times \left[s(u,v)-s(u,z)s(z,v)\right]=2N_{c}\int_{z}K_{JSJ}(u,v;z)\left[s(u,z)s(z,v)-s(u,v)\right]\,,
\end{split}\end{equation}
```

## \(K_{JSSJ}\)

Subsection source range: `references/WORKNLO.tex` lines 1156--1166; formula
body lines 1158--1164.

Combined-equation source range: `references/WORKNLO.tex` lines 1353--1355.

Classification: the subsection gives the isolated action of the
\(f f J_L S S J_R\) part. The full \(K_{JSSJ}\) Hamiltonian also contains the
\(-N_c J_L S_A J_R\) subtraction; the fully combined dipole equation below
contains the complete \(K_{JSSJ}\) coefficient in the combination
\(K_{JSSJ}-\widetilde K\).

Subsection formula:

```tex
\begin{equation}\begin{split}
&\int_{x\, y\, z\, z^{\prime}}K_{JSSJ}(x,y;z,z^{\prime})\left[f^{abc}f^{def}J_{L}^{a}(x)S_A^{be}(z)S_A^{cf}(z^{\prime})J_{R}^{d}(y)\right]s(u,v)=\\
&=\frac{1}{4}\int_{x\, y\, z\, z^{\prime}}\,\left(K_{JSSJ}(v,v;z,z^{\prime})-K_{JSSJ}(u,v;z,z^{\prime})-K_{JSSJ}(v,u;z,z^{\prime})+K_{JSSJ}(u,u;z,z^{\prime})\right)\\
&\times \left[N_{c}^{2}s(u,z^{\prime})s(z^{\prime},z)s(z,v)+N_{c}^{2}s(u,z)s(z,z^{\prime})s(z^{\prime},v)-\frac{1}{N_{c}}tr\left[S^{\dagger}(u)S(z^{\prime})S^{\dagger}(z)S(v)S^{\dagger}(z^{\prime})S(z)\right]\right.\\
&\left.-\frac{1}{N_{c}}tr\left[S^{\dagger}(u)S(z)S^{\dagger}(z^{\prime})S(v)S^{\dagger}(z)S(z^{\prime})\right]\right]\\
&=\frac{1}{2}\int_{z\, z^{\prime}}\,\left(K_{JSSJ}(u,u;z,z^{\prime})+K_{JSSJ}(v,v;z,z^{\prime})-K_{JSSJ}(v,u;z,z^{\prime})-K_{JSSJ}(u,v;z,z^{\prime})\right)\\
&\times \left[N_{c}^{2}s(u,z^{\prime})s(z^{\prime},z)s(z,v)-\frac{1}{N_{c}}tr\left[S^{\dagger}(u)S(z)S^{\dagger}(z^{\prime})S(v)S^{\dagger}(z)S(z^{\prime})\right]\right]\\
&=-\int_{z\, z^{\prime}}\, K_{JSSJ}(u,v;z,z^{\prime})\left[N_{c}^{2}s(u,z^{\prime})s(z^{\prime},z)s(z,v)-\frac{1}{N_{c}}tr\left[S^{\dagger}(u)S(z)S^{\dagger}(z^{\prime})S(v)S^{\dagger}(z)S(z^{\prime})\right]\right]\\
\end{split}\end{equation}
```

Complete \(K_{JSSJ}\) occurrence in the combined equation:

```tex
&&-\frac{1}{N_{c}}\int_{z,z^{\prime}}\,\left[K_{JSSJ}(u,v;z,z^{\prime})-\widetilde{K}(u,v,z,z^{\prime})\right]
\left[N_{c}^{3}s(u,z^{\prime})s(z^{\prime},z)s(z,v)\right.\nonumber \\
&&\left. -tr\left[S(v)S^{\dagger}(z)S(z^{\prime})S^{\dagger}(u)S(z)S^{\dagger}(z^{\prime})\right]-N_{c}^{3}s(u,z)s(z,v)+N_{c}s(u,v)\right]\nonumber \\
```

## \(K_{q\bar q}\)

Subsection source range: `references/WORKNLO.tex` lines 1172--1181; formula
body lines 1174--1179.

Combined-equation source range: `references/WORKNLO.tex` line 1352.

Classification: the subsection isolates the quark trace-current term. The
same sector also appears in the combined equation in trace-product form with
explicit \((z\rightarrow z^\prime)\) exchange notation.

Implementation note: the finite-grid diagnostic target uses the exact
generator trace-product expression in lines 1174--1177. The later compact
reduction is retained as source text, but it is not used as the diagnostic
implementation because it introduced a finite-grid mismatch for four or more
sampled sites in `reports/nlo_current/kqbarq_physical_finite_grid_diagnosis.md`.

Subsection formula:

```tex
\begin{equation}\begin{split}
&\int_{x,y,z,z^{\prime}}K_{q\bar{q}}(x,y;z,z^{\prime})\left[2J_{L}^{a}(x)tr\left[S^{\dagger}(z)t^{a}S(z^{\prime})t^{b}\right]J_{R}^{b}(y)\right]s(u,v)\\
&=\frac{1}{N_{c}}\int_{z\, z^{\prime}}\,\left(K_{q\bar{q}}(u,u;z,z^{\prime})+K_{q\bar{q}}(v,v;z,z^{\prime})-K_{q\bar{q}}(v,u;z,z^{\prime})-K_{q\bar{q}}(u,v;z,z^{\prime})\right)\\
&\times 2tr\left[S^{\dagger}(u)t^{a}S(v)t^{b}\right]tr\left[S^{\dagger}(z)t^{a}S(z^{\prime})t^{b}\right]\\
&=\frac{1}{2N_{c}}\int_{z\, z^{\prime}}\,\left(K_{q\bar{q}}(u,u;z,z^{\prime})+K_{q\bar{q}}(v,v;z,z^{\prime})-K_{q\bar{q}}(v,u;z,z^{\prime})-K_{q\bar{q}}(u,v;z,z^{\prime})\right)\\
&\times \Big(N_{c}s(u,z^{\prime})s(z,v)-\frac{1}{N_{c}^{2}}tr\left[S^{\dagger}(u)S(v)S^{\dagger}(z)S(z^{\prime})\right]-\frac{1}{N_{c}^{2}}tr\Big[S^{\dagger}(u)S(v)S^{\dagger}(z^{\prime})S(z)\Big]\\ &+\frac{1}{N_{c}}s(u,v)s(z,z^{\prime})\Big)= -\int_{z\, z^{\prime}}\, K_{q\bar{q}}(u,v;z,z^{\prime})\times\\
&\times\Big(N_{c}s(u,z^{\prime})s(z,v)-\frac{1}{N_{c}^{2}}tr\left[S^{\dagger}(u)S(v)S^{\dagger}(z)S(z^{\prime})\Big]-\frac{1}{N_{c}^{2}}tr\Big[S^{\dagger}(u)S(v)S^{\dagger}(z^{\prime})S(z)\right]\\ & +\frac{1}{N_{c}}s(u,v)s(z,z^{\prime})\Big)\,.
\end{split}\end{equation}
```

Combined-equation form:

```tex
&&-\frac{4}{N_{c}}\int_{z\, z^{\prime}}\, K_{q\bar{q}}(u,v;z,z^{\prime})\left[tr\left[S^{\dagger}(u)t^{a}S(v)t^{b}\right]tr\left[S^{\dagger}(z)t^{a}S(z^{\prime})t^{b}\right]-(z\rightarrow z^{\prime})\right]\nonumber \\
```

## \(K_{JJSJ}\)

Real source range: `references/WORKNLO.tex` lines 1208--1224; compact real
formula lines 1216--1223.

Virtual source range: `references/WORKNLO.tex` lines 1244--1253; formula body
lines 1246--1251.

Combined-equation source range: `references/WORKNLO.tex` lines 1350--1351.

Classification: isolated real and virtual sector contributions are available
in the subsection. The real contribution uses the antisymmetry of
\(K_{JJSJ}(w,x,y;z)=-K_{JJSJ}(w,y,x;z)\). The combined equation uses the
\(K_{JJSJ}\) contribution together with \(K_{JSJ}\) and \(\widetilde K\) in the
first line.

Real \(LLR-LRR\) contribution:

```tex
\begin{eqnarray}
&&\int_{w,x,y,z}K_{JJSJ}(w;x,y;z)f^{bde}\left[J_{L}^{d}(x)J_{L}^{e}(y)S_A^{ba}(z)J_{R}^{a}(w)-J_{L}^{a}(w)S_A^{ab}(z)J_{R}^{d}(x)J_{R}^{e}(y)\right]s(u,v)\nonumber \\
&&=\frac{i}{2}\int_{z}\left.[
%K_{JJSJ}(v;v,v;z)+
K_{JJSJ}(v;u,v;z)-K_{JJSJ}(v;v,u;z)
%-K_{JJSJ}(v;u,u;z)-K_{JJSJ}(u;v,v;z)
-K_{JJSJ}(u;u,v;z)+K_{JJSJ}(u;v,u;z)
%+K_{JJSJ}(u;u,u;z)
\right]\nonumber \\
&&\times \left[s(u,v)-N_{c}^{2}s(u,z)s(z,v)\right]\label{rJJSJ}\,.
 \end{eqnarray}
```

Virtual contribution:

```tex
\begin{equation}\begin{split}
&\frac{1}{3}\int_{w,x,y,z}K_{JJSJ}(w;x,y;z)f^{bde}\left[J_{L}^{d}(x)J_{L}^{e}(y)J_{L}^{b}(w)-J_{R}^{d}(x)J_{R}^{e}(y)J_{R}^{b}(w)\right]s(u,v)\\
&=i\frac{N_{c}^{2}-1}{3}\int_{z}[K_{JJSJ}(u,v,u,z)+K_{JJSJ}(v,u,v,z) ]s(u,v)\,. \\
 \end{split}\end{equation}
```

Combined-equation occurrence:

```tex
&&=\int_{z}\left[2N_{c}K_{JSJ}(u,v;z)-iN_{c}^{2}\left(K_{JJSJ}(v,u,v,z)+K_{JJSJ}(u,v,u,z)\right)+N_{c}^{2}\int_{z^{\prime}}\widetilde{K}(u,v,z,z^{\prime})\right]\nonumber\\
&&\times \left[s(u,z)s(z,v)-s(u,v)\right]\nonumber\\
```

## \(K_{JJSSJ}\)

Real source range: `references/WORKNLO.tex` lines 1304--1328; final formula
body lines 1308--1326. The preliminary separate LLR/LRR expressions are lines
1288--1303.

Virtual source range: `references/WORKNLO.tex` lines 1334--1341; formula body
lines 1336--1340.

Combined-equation source range: `references/WORKNLO.tex` lines 1344--1360;
`ourdipole` formula body lines 1348--1359.

Classification: isolated real and virtual sector contributions are available
in the subsection, but the real contribution is partly expressed using
\(\widetilde K\). The combined equation also mixes \(K_{JJSSJ}\) through
\(\widetilde K\) with the \(K_{JSJ}\), \(K_{JJSJ}\), and \(K_{JSSJ}\)
structures, and includes a separate pure \(K_{JJSSJ}\) eight-kernel
combination.

Real \(LLR-LRR\) contribution:

```tex
\begin{eqnarray}
&&\int_{w,x,y,z,z^{\prime}}K_{JJSSJ}(w;x,y;z,z^{\prime})f^{acb}\nonumber \\
&&\ \ \ \times \left[J_{L}^{d}(x)J_{L}^{e}(y)S_A^{dc}(z)S_A^{eb}(z^{\prime})J_{R}^{a}(w)-J_{L}^{a}(w)S_A^{cd}(z)S_A^{be}(z^{\prime})J_{R}^{d}(x)J_{R}^{e}(y)\right]s(u,v)\nonumber 
\\
&&=\frac{i}{2N_{c}}\int_{z,z^{\prime}}\left[K_{JJSSJ}(u;u,u;z,z^{\prime})-K_{JJSSJ}(u;v,u;z,z^{\prime})+K_{JJSSJ}(u;v,v;z,z^{\prime})\right. \nonumber \\
&& \ \ \ -K_{JJSSJ}(u;u,v;z,z^{\prime})+K_{JJSSJ}(v;u,v;z,z^{\prime})-K_{JJSSJ}(v;u,u;z,z^{\prime})\nonumber \\
&&\left. \ \ \ +K_{JJSSJ}(v;v,u,z,z^{\prime})-K_{JJSSJ}(v;v,v;z,z^{\prime})\right]
N_{c}^{3}s(z,v)s(z^{\prime},z)s(u,z^{\prime})\nonumber \\
&&\ \ \ +\frac{1}{N_{c}}\int_{z,z^{\prime}}\widetilde{K}(u,v,z,z^{\prime})\times\left[N_{c}^{3}s(z,v)s(z^{\prime},z)s(u,z^{\prime})-tr\left[S(v)S^{\dagger}(z)S(z^{\prime})S^{\dagger}(u)S(z)S^{\dagger}(z^{\prime})\right]\right]\, , \nonumber \\
\end{eqnarray}
```

Virtual contribution:

```tex
\begin{eqnarray}
&&\frac{1}{3}\int_{w,x,y,z,z^{\prime}}K_{JJSSJ}(w;x,y;z,z^{\prime})f^{acb}\left[J_{L}^{c}(x)J_{L}^{b}(y)J_{L}^{a}(w)-J_{R}^{c}(x)J_{R}^{b}(y)J_{R}^{a}(w)\right]s(u,v)\nonumber \\
&&=-\int_{z,z^{\prime}}\left[-K_{JJSSJ}(u;v,v;z,z^{\prime})+K_{JJSSJ}(u;u,v;z,z^{\prime})+K_{JJSSJ}(v;v,v;z,z^{\prime})\right.\nonumber \\
&&-K_{JJSSJ}(v;u,v;z,z^{\prime})+K_{JJSSJ}(v;v,u;z,z^{\prime})-K_{JJSSJ}(v;u,u;z,z^{\prime})\nonumber \\
&&\left.+K_{JJSSJ}(u;u,u;z,z^{\prime})-K_{JJSSJ}(u;v,u;z,z^{\prime})\right]\frac{N_{c}^{2}-1}{6}is(u,v)\nonumber\\
&&=-\frac{(N_{c}^{2}-1)}{3}\,\int_{z,z^{\prime}}\widetilde{K}(u,v,z,z^{\prime})s(u,v)\,.
  \end{eqnarray}
```

Pure \(K_{JJSSJ}\) occurrence in the combined equation:

```tex
&&+\frac{i}{2N_{c}}\int_{z,z^{\prime}}\left[K_{JJSSJ}(u,u,u,z,z^{\prime})-K_{JJSSJ}(u,v,u,z,z^{\prime})+K_{JJSSJ}(u,v,v,z,z^{\prime})\right. \nonumber \\
&&-K_{JJSSJ}(u,u,v,z,z^{\prime})+K_{JJSSJ}(v,u,v,z,z^{\prime})-K_{JJSSJ}(v,u,u,z,z^{\prime})\nonumber \\
&&\left.+K_{JJSSJ}(v,v,u,z,z^{\prime})-K_{JJSSJ}(v,v,v,z,z^{\prime})\right] \,N_{c}^{3}s(z,v)s(z^{\prime},z)s(u,z^{\prime})\,.
```

## Target Availability Summary

| sector | availability status | notes |
|---|---:|---|
| \(K_{JSJ}\) | full target, implemented | Requires symmetric endpoint kernel and zero diagonal condition stated in TeX. |
| \(K_{JSSJ}\) | partial subsection; full combined-equation target required | The \(f f J_L S S J_R\) part is isolated; the full subtraction-complete sector appears in the combined equation through \(K_{JSSJ}-\widetilde K\). |
| \(K_{q\bar q}\) | subsection target is partial unless subtraction is included; combined-equation target available | The subsection isolates the quark trace-current term. The direct local full-sector current also includes the \(-J_L S_A J_R\) subtraction. |
| \(K_{JJSJ}\) | real + virtual isolated target available, convention diagnosis pending | Real and virtual pieces are isolated; combined equation also mixes the endpoint contribution with \(K_{JSJ}\) and \(\widetilde K\). |
| \(K_{JJSSJ}\) | real + virtual target available with \(\widetilde K\) complications, convention diagnosis pending | Real and virtual pieces are isolated but use \(\widetilde K\); combined equation contains both \(\widetilde K\) and a pure eight-kernel combination. |
