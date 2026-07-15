# K_JSSJ Coefficient Symmetry Report

This is a dense small-lattice diagnostic, not a proof. It tests whether the
left-basis coefficient C^{(x,a)(y,h)} built from the Appendix B ordered
A_JSSJ block is symmetric under combined index exchange (x,a)<->(y,h).

| seed | n_site | kernel symmetry | ||C|| | ||C-C^T|| | r_asym |
|---:|---:|---|---:|---:|---:|
| 20260704 | 3 | K(x,y;z,z')=K(y,x;z,z') | 8.1042804532210475e+00 | 1.0895716987442919e+01 | 1.3444397747997991e+00 |
| 20260705 | 3 | K(x,y;z,z')=K(y,x;z,z') and K(x,y;z,z')=K(x,y;z',z) | 7.6955603611530101e+00 | 9.9681317280075081e+00 | 1.2953094070090567e+00 |

## Conclusion

max_r_asym: 1.3444397747997991e+00

C_JSSJ is not symmetric in these tests. Do not simplify it to a pure symmetric score-current; keep the antisymmetric/commutator drift.
