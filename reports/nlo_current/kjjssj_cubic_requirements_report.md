# K_JJSSJ Cubic Requirements Report

Distinct-site dense small-lattice diagnostic for the final cubic NLO block.
Coincident-site commutators are explicitly not resolved here.

## Kernel klm_antisym

seed: 20260801
nsite: 3
K(w;x,y;z,z') + K(w;y,x;z',z) max residual: 0.0000000000000000e+00

Built coefficient blocks: A_LLR, B_LRR, V_virtual.

| block | norm |
|---|---:|
| A_LLR | 4.8705215401232493e+00 |
| B_LRR | 4.8705215401232511e+00 |
| V_virtual | 1.0718504778712232e+00 |

| component | second derivative | score-linear | Hessian-score | score-product |
|---|---:|---:|---:|---:|
| LLR | 4.4749041139269252e-03 | 1.2807954745390762e-03 | 3.3624265535648586e-05 | 4.0603259195243788e-04 |
| LRR | 3.4937662329176389e-03 | 2.4705241275422406e-03 | 1.7493785691600776e-05 | 2.1124765208300100e-04 |
| virtual | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 5.1891223835296681e-05 | 6.2661675363857795e-04 |

hessian_probe_abs: 1.0534461081235236e-03

| sign check | relative residual |
|---|---:|
| llr | 2.3060939822235754e-10 |
| lrr | 0.0000000000000000e+00 |
| virtual_lll | 0.0000000000000000e+00 |
| virtual_rrr | 1.9006860972267670e-11 |

## Kernel unconstrained

seed: 20260802
nsite: 3

Built coefficient blocks: A_LLR, B_LRR, V_virtual.

| block | norm |
|---|---:|
| A_LLR | 4.8523084331530777e+00 |
| B_LRR | 4.8523084331530715e+00 |
| V_virtual | 1.2648123357399943e+00 |

| component | second derivative | score-linear | Hessian-score | score-product |
|---|---:|---:|---:|---:|
| LLR | 1.4235249159999519e-03 | 3.3035866452428781e-04 | 1.4228378644586717e-05 | 2.6203575906337660e-05 |
| LRR | 1.1110351527676698e-04 | 7.0993373907178364e-04 | 1.9430741793427697e-05 | 3.5784464992027862e-05 |
| virtual | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 5.1709531249785783e-05 | 9.5230430749072785e-05 |

hessian_probe_abs: 5.7898814909140839e-04

| sign check | relative residual |
|---|---:|
| llr | 4.5529412514427881e-09 |
| lrr | 0.0000000000000000e+00 |
| virtual_lll | 0.0000000000000000e+00 |
| virtual_rrr | 1.9441355128416507e-10 |

## Conclusion

max_sign_residual: 4.5529412514427881e-09
max_hessian_score_component: 5.1891223835296681e-05

K_JJSSJ distinct-site current can be represented as score + Hessian-score with the tested signs.

Coincident-site commutators remain unresolved and must be handled in a separate workflow before production use.
