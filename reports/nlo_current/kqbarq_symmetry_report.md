# K_qbarq Coefficient Symmetry Report

This dense small-lattice diagnostic tests whether the left-basis
coefficient C^{(x,a)(y,h)} built from the ordered K_qbarq A block is
symmetric under combined index exchange (x,a)<->(y,h). It is not a proof.

| seed | n_site | kernel symmetry | ||C|| | ||C-C^T|| | r_asym | max Im A | max Im C |
|---:|---:|---|---:|---:|---:|---:|---:|
| 20260714 | 3 | K(x,y;z,z')=K(y,x;z,z') | 2.9878690619841626e+00 | 3.2654161716689138e+00 | 1.0928913228548376e+00 | 2.8161583942624385e-01 | 2.4819880778710335e-01 |
| 20260715 | 3 | K(x,y;z,z')=K(y,x;z,z') and K(x,y;z,z')=K(x,y;z',z) | 3.5522679567587927e+00 | 2.8103659277747437e+00 | 7.9114694104861816e-01 | 0.0000000000000000e+00 | 0.0000000000000000e+00 |

## Conclusion

max_r_asym: 1.0928913228548376e+00

C_qbarq has an order-one antisymmetric component in these tests. Keep the generic ordered-current / commutator-drift representation.
