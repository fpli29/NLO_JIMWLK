# Ordered LR Current Finite-Difference Check

Checks -J_R^b J_L^a(A^{ab}W) against -L_y^h{S_y^{hb} L_x^a[A^{ab}W]} on a two-site SU(3) toy problem. The inner L_x derivative is evaluated analytically for the toy A and W; the outer left/right derivatives are forward finite differences.

seed: 7890

| eps | direct | divergence | relative residual |
|---:|---:|---:|---:|
| 1e-04 | -9.6106647624825941e-03 | -9.6096690029457560e-03 | 9.9575953683810059e-07 |
| 1e-05 | -9.6097725847590198e-03 | -9.6096730117616030e-03 | 9.9572997416808029e-08 |
| 1e-06 | -9.6096833876306055e-03 | -9.6096734363351755e-03 | 9.9512954300495782e-09 |
