# K_JJSJ Cubic Requirements Report

Dense small-lattice diagnostic for the first three-generator NLO block.
The tested current signs are distinct-site ordered-block checks, not a
production treatment of all coincident-site commutators.

## Kernel symmetric

seed: 20260724
nsite: 3

Built coefficient blocks: A_LLR, B_LRR, V_virtual.

| block | norm |
|---|---:|
| A_LLR | 4.9006773840965767e+00 |
| B_LRR | 4.9006773840965945e+00 |
| V_virtual | 1.6468871288320945e+00 |

| component | second derivative | score-linear | Hessian-score | score-product |
|---|---:|---:|---:|---:|
| LLR | 0.0000000000000000e+00 | 1.5777727038483931e-04 | 6.8662788190954074e-05 | 8.4394265875573757e-06 |
| LRR | 0.0000000000000000e+00 | 1.5777727038483931e-04 | 6.8662788190954074e-05 | 8.4394265875573757e-06 |
| virtual | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 4.2431222464970789e-04 | 5.2152730241241978e-05 |

hessian_probe_abs: 5.8478200094802535e-03

| sign check | relative residual |
|---|---:|
| llr | 3.0319018985962981e-10 |
| lrr | 0.0000000000000000e+00 |
| virtual_lll | 0.0000000000000000e+00 |
| virtual_rrr | 6.7563524631768235e-10 |

## Kernel antisymmetric

seed: 20260725
nsite: 3

Built coefficient blocks: A_LLR, B_LRR, V_virtual.

| block | norm |
|---|---:|
| A_LLR | 4.8895255827129933e+00 |
| B_LRR | 4.8895255827129924e+00 |
| V_virtual | 1.2582937403350736e+00 |

| component | second derivative | score-linear | Hessian-score | score-product |
|---|---:|---:|---:|---:|
| LLR | 3.4694469519536142e-12 | 1.3328328185997611e-03 | 1.1358480138403666e-03 | 8.3367118337033969e-05 |
| LRR | 3.4694469519536142e-12 | 1.3328328185997611e-03 | 1.1358480138403666e-03 | 8.3367118337033969e-05 |
| virtual | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 8.4625985261482057e-04 | 6.2112399209368055e-05 |

hessian_probe_abs: 1.7848910623813330e-02

| sign check | relative residual |
|---|---:|
| llr | 3.4975497433149816e-11 |
| lrr | 0.0000000000000000e+00 |
| virtual_lll | 0.0000000000000000e+00 |
| virtual_rrr | 7.0905992665890039e-11 |

## Conclusion

max_sign_residual: 6.7563524631768235e-10
max_hessian_score_component: 1.1358480138403666e-03

K_JJSJ current can be represented as score + Hessian-score with the tested distinct-site ordered signs. Coincident-site production commutator handling remains outside this workflow.
