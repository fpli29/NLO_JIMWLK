# Dipole Validation Status

This document summarizes the non-production observable-side dipole validation.

## Direct generator-action validation

The analytic dipole generator actions use the local convention

\[
L_x^aF(U)=\frac{d}{d\epsilon}F(e^{i\epsilon t^a}U_x)\bigg|_{\epsilon=0},
\qquad
R_x^aF(U)=\frac{d}{d\epsilon}F(U_xe^{i\epsilon t^a})\bigg|_{\epsilon=0}.
\]

Tests compare analytic single-, two-, and three-generator words against nested
finite differences for representative `LL`, `LR`, `RL`, `RR`, `LLR`, `LRR`,
`LLL`, and `RRR` words.

## Appendix A target status

| sector | local Appendix A target implemented | status | max residual |
|---|---:|---|---:|
| \(K_{JSJ}\) | yes | passed | `1.6172982426630268e-16` |
| \(K_{JSSJ}\) | yes | passed | `1.3877787807814457e-16` |
| \(K_{q\bar q}\) | yes | passed | `3.236828524569469e-16` |
| \(K_{JJSJ}\) | yes | passed | `1.1801832636420706e-15` |
| \(K_{JJSSJ}\) | yes | passed | `1.9087382418229414e-15` |

The \(K_{JSJ}\) comparison uses an Appendix-compatible symmetric synthetic
kernel with zero diagonal in the endpoint coordinates. For a raw arbitrary
`K(x,y;z)`, direct action also receives endpoint-diagonal contributions; those
are not the closed Appendix A dipole target \(K(u,v;z)\).

The \(K_{q\bar q}\) full target uses the exact Hamiltonian bracket in
`WORKNLO.tex` lines 268--269: the trace-current term plus the
\(-J_L S_A J_R\) subtraction. The subsection formula at lines 1174--1179 is a
partial target for the trace-current component only. Tests verify that treating
the subsection expression as a full-sector target fails.

The \(K_{JSSJ}\) full target uses only the \(K_{JSSJ}\) part of the combined
Appendix equation in `WORKNLO.tex` lines 1353--1355. The subsection expression
at lines 1158--1164 is partial because it omits the \(-N_c J_L S_A J_R\)
subtraction. The \(\widetilde K\) contribution in the combined factor is
excluded from the isolated \(K_{JSSJ}\) target and remains part of the
\(K_{JJSSJ}\)-related bookkeeping.

The \(K_{JJSJ}\) target uses the exact isolated real contribution from
`WORKNLO.tex` lines 1208--1224 and the isolated virtual contribution from lines
1244--1253. Comparisons use the calibrated cubic convention
`Appendix target = (-i) * raw direct cubic action`; existing direct-action
functions were not redefined.

The \(K_{JJSSJ}\) target uses the exact \(\widetilde K\) definition from
`WORKNLO.tex` lines 307--311, the isolated real contribution from lines
1304--1328, and the isolated virtual contribution from lines 1334--1341.
Comparisons use the same calibrated cubic convention. Tests exercise both the
pure eight-kernel real term and the \(\widetilde K\) real/virtual terms.

## Internal checks

- Direct Hamiltonian actions give zero for zero kernels.
- Direct Hamiltonian actions are linear in each sector kernel.
- The \(q\bar q\) \(z'=z\) subtraction consistency test passes.
- The cubic virtual \(1/3\) sensitivity test passes for both \(K_{JJSJ}\) and
  \(K_{JJSSJ}\).

## Sign and normalization issues

No sign mismatch was found for the implemented \(K_{JSJ}\), full \(K_{JSSJ}\),
full \(K_{q\bar q}\), calibrated \(K_{JJSJ}\), or calibrated \(K_{JJSSJ}\)
Appendix targets under the K_JSJ-calibrated convention `target = H_sector s`.

The cubic one-\(f\) convention diagnostic in
`reports/nlo_current/cubic_i_convention_calibration_report.md` verifies on the
tested small lattices that `TeX target = (-i) * direct_action` for K_JJSJ real,
K_JJSJ virtual, and K_JJSSJ-like real/virtual isolated blocks. K_JJSJ and
K_JJSSJ are now implemented and tested with this calibration.

## Remaining work

- Integrate all validated sector targets into broader non-production dipole
  diagnostics as needed.
- Integrate physical kernels after sector targets pass.
- Keep this separate from production evolution and score/Hessian-score model
  design.
