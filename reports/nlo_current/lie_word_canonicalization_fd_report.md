# Lie Word Canonicalization Finite-Difference Report

Residuals at eps=1e-5 are expected to be roundoff dominated for nested third derivatives.

| pattern | eps | direct | canonical | relative residual |
|---|---:|---:|---:|---:|
| same_site_second | 1e-03 | -1.3205585425168920e-01 | -1.3205586489073440e-01 | 1.0639045200377950e-08 |
| same_site_second | 1e-04 | -1.3205585513986762e-01 | -1.3205585291498068e-01 | 2.2248869413488137e-09 |
| same_site_second | 1e-05 | -1.3205658788706387e-01 | -1.3205292201945440e-01 | 3.6658676094702969e-06 |
| same_site_third | 1e-03 | 4.7828885296752333e-01 | 4.7828946470040989e-01 | 6.1173288656846125e-07 |
| same_site_third | 1e-04 | 4.7839510131097995e-01 | 4.7796098190389102e-01 | 4.3411940708892871e-04 |
| same_site_third | 1e-05 | 2.2204460492503128e-01 | 6.2350791196763555e-01 | 4.0146330704260424e-01 |

stable_max_same_site_second: 1.1004630540156768e-06
stable_max_same_site_third: 3.9857184219727060e-06
