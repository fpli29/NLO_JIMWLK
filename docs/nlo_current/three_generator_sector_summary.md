# Three-Generator NLO Current Sector Summary

This summary covers the distinct-site small-lattice validation of the cubic
sector. It does not implement the full NLO flow and does not solve coincident
site commutators.

## \(K_{JJSJ}\) status

The \(K_{JJSJ}\) workflow validated the distinct-site ordered cubic structure:

- LLR sign passed.
- LRR sign passed.
- Virtual LLL/RRR sign smoke tests passed.
- The virtual \(1/3\) factor is included and tested.
- Hessian-score contributions are nonzero.

The tested synthetic kernels included symmetric and antisymmetric
\(x\leftrightarrow y\) conventions. The antisymmetric synthetic kernel gave the
expected combined \((x,d)\leftrightarrow(y,e)\) behavior.

Coincident-site commutators remain unresolved.

## \(K_{JJSSJ}\) status

The \(K_{JJSSJ}\) distinct-site workflow added the coefficient blocks:

\[
A_{LLR}^{dea}(x,y,w)
=
\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}S_z^{dc}S_{z'}^{eb},
\]

\[
B_{LRR}^{ade}(w,x,y)
=
-\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')
f^{acb}S_z^{cd}S_{z'}^{be},
\]

\[
V^{cba}(x,y,w)
=
\frac13\int_{z,z'}
K_{JJSSJ}(w;x,y;z,z')f^{acb}.
\]

The primary synthetic kernel imposed the KLM-like simultaneous antisymmetry
\[
K(w;x,y;z,z')=-K(w;y,x;z',z).
\]

The kernel residual for this identity was `0.0`. The induced LLR coefficient
behavior was measured rather than assumed; the report found that \(A_{LLR}\)
is symmetric under combined \((x,d)\leftrightarrow(y,e)\) exchange for this
synthetic kernel, with minus residual `5.5511151231257827e-17`.

Distinct-site sign checks passed:

- LLR max residual in diagnostic: `2.3060939822235754e-10` for the KLM-like kernel.
- LRR residual in diagnostic: `0.0`.
- virtual LLL residual: `0.0`.
- virtual RRR residual: `1.9006860972267670e-11` for the KLM-like kernel.

The virtual \(1/3\) factor was included and tested by a factor-three norm ratio.

Hessian-score terms were nonzero. For the KLM-like kernel, the report measured:

- LLR Hessian-score component: `3.3624265535648586e-05`.
- LRR Hessian-score component: `1.7493785691600776e-05`.
- virtual Hessian-score component: `5.1891223835296681e-05`.

Coincident-site commutators remain unresolved:
\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

## Shared conclusion

The NLO three-generator sector requires score plus Hessian-score:

\[
s_A=L_A\log W,
\qquad
H_{AB}=L_As_B=L_AL_B\log W.
\]

The distinct-site ordered signs are validated for \(K_{JJSJ}\) and
\(K_{JJSSJ}\), but this is not a production current implementation.

## Coincident-site commutator status

The canonical Lie-word ordering rule is lexicographic in the combined index
\[
A=(x,a).
\]

For adjacent same-site derivatives, the validated swap rule is
\[
L_x^aL_x^b
=
L_x^bL_x^a
+
f^{abc}L_x^c.
\]

For different sites, derivatives commute:
\[
L_x^aL_y^b=L_y^bL_x^a,\qquad x\neq y.
\]

Symbolic tests passed for:

- different-site commutation;
- same-site two-word commutator;
- three-letter path independence;
- canonical words left unchanged;
- order splitting into cubic, quadratic, linear, and scalar pieces.

Finite-difference tests passed for same-site and mixed-site canonicalization.
The requested \(10^{-5}\) finite-difference entries are reported but are
roundoff dominated for nested third derivatives; stable residuals are reported
from the \(10^{-2}\) to \(10^{-3}\) range.

End-to-end coincident patterns passed:

- \(x=y\neq w\): stable max residual `4.4200754167889045e-09`.
- \(x=w\neq y\): stable max residual `4.4200754167889045e-09`.
- \(y=w\neq x\): stable max residual `4.4200754167889045e-09`.
- \(x=y=w\): stable max residual `2.6429552990592242e-08`.

Coefficient-level diagnostics found commutator-induced lower-order terms:

- \(K_{JJSJ}\): quadratic correction norm `6.1541551762714866e+00`;
  linear correction norm `0.0`.
- \(K_{JJSSJ}\): quadratic correction norm `1.0455495097967299e+01`;
  linear correction norm `3.6795897670206137e+00`.

Thus commutator-induced \(K_{2,\rm comm}\) corrections are nonzero for both
cubic kernels in the synthetic diagnostics, and \(K_{1,\rm comm}\) is nonzero
for \(K_{JJSSJ}\).

Production flow still requires assembling these corrections into the final NLO
current implementation with the correct global Hamiltonian signs and coefficient
normalizations.

## Remaining work

- Dipole validation beyond skeleton formulas.
- Non-production NLO current skeleton that includes \(K_{2,\rm comm}\) and
  \(K_{1,\rm comm}\).
- Eventual production-flow design after all ordering and commutator terms are
  controlled.
