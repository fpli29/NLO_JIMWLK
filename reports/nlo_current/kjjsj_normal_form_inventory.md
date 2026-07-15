# KJJSJ Normal-Form Inventory

## Scope

This inventory records every current \(K_{JJSJ}\) contribution to the dense
normal-form tensors before implementing analytic derivatives.

## 1. Distinct-Site LLR K3

name: `KJJSJ LLR`

source function: `three_generator_terms.kjjsj_A_LLR_from_kernel`, followed by
`cubic_commutator_terms.llr_right_to_left` and
`cubic_block_terms_from_LLR`.

tensor order: cubic raw word `((x,d),(y,e),(w,h))`.

index ordering: `A_LLR[x,y,w,d,e,a]`, converted to
`C_left[x,y,w,d,e,h]`.

Wilson-line dependence: one adjoint Wilson line inside `A_LLR` and one adjoint
right-to-left conversion factor at site `w`.

coordinate-kernel dependence: normalized `KJJSJ[w,x,y,z]`.

raw/normalized dtype: physical adapter supplies KLM-normalized real
coefficients for the tested physical kernels.

commutator origin: same-site word canonicalization can generate lower-order
quadratic/linear terms.

expected derivative contribution: contributes to `LC_K3`, `LB_K3`, `d2K3`,
and to `dK2_comm` through canonicalized quadratic terms.

## 2. Distinct-Site LRR K3

name: `KJJSJ LRR`

source function: `three_generator_terms.kjjsj_B_LRR_from_kernel`, followed by
`cubic_commutator_terms.lrr_right_to_left` and
`cubic_block_terms_from_LRR`.

tensor order: cubic raw word `((w,a),(x,p),(y,q))`.

index ordering: `B_LRR[w,x,y,a,d,e]`, converted to
`C_left[w,x,y,a,p,q]`.

Wilson-line dependence: one adjoint Wilson line inside `B_LRR` and two adjoint
right-to-left conversion factors at sites `x` and `y`.

coordinate-kernel dependence: normalized `KJJSJ[w,x,y,z]`.

raw/normalized dtype: KLM-normalized physical coefficient.

commutator origin: same-site canonicalization of the left-basis cubic word.

expected derivative contribution: contributes to `LC_K3`, `LB_K3`, `d2K3`,
and `dK2_comm`.

## 3. Virtual Cubic K3

name: `KJJSJ virtual LLL/RRR`

source function: `three_generator_terms.kjjsj_V_virtual_from_kernel`, followed
by `cubic_block_terms_from_virtual_LLL` and
`virtual_rrr_right_to_left`/`cubic_block_terms_from_virtual_RRR`.

tensor order: LLL raw word `((x,d),(y,e),(w,b))`; RRR converted raw word
`((x,p),(y,q),(w,h))`.

index ordering: `V[x,y,w,d,e,b]`; RRR conversion uses three adjoint factors.

Wilson-line dependence: LLL virtual block has no Wilson-line dependence. RRR
virtual block has three adjoint conversion factors.

coordinate-kernel dependence: normalized `KJJSJ[w,x,y,z]` summed over `z`
with the existing `1/3` factor.

raw/normalized dtype: KLM-normalized physical coefficient.

commutator origin: canonicalization of same-site words in both LLL and RRR
virtual terms.

expected derivative contribution: LLL derivative is structurally zero. The RRR
right-to-left product is differentiated, but in the current canonicalized
diagnostic its `K3` derivative contractions are structurally zero; tiny
roundoff-level `dK2_comm` entries are treated by absolute residuals rather
than relative-error claims.

## 4. Quadratic Commutator Correction K2_comm

name: `KJJSJ quadratic commutator correction`

source function: `cubic_commutator_terms.canonicalize_cubic_block_terms`.

tensor order: quadratic word after canonicalization.

index ordering: flattened through `_word_terms_to_arrays`.

Wilson-line dependence: inherited from differentiated raw block coefficients.

coordinate-kernel dependence: inherited normalized `KJJSJ`.

raw/normalized dtype: follows block coefficients.

commutator origin: same-site left-derivative canonicalization using
`[L_x^a,L_x^b]=f^{abc}L_x^c`.

expected derivative contribution: contributes to `dK2=L_BK2^{AB}`.

## 5. Linear Commutator Correction K1_comm

name: `KJJSJ linear commutator correction`

source function: `cubic_commutator_terms.canonicalize_cubic_block_terms`.

tensor order: linear word after canonicalization.

index ordering: flattened through `_word_terms_to_arrays`.

Wilson-line dependence: inherited from raw blocks if nonzero.

coordinate-kernel dependence: inherited normalized `KJJSJ`.

raw/normalized dtype: follows block coefficients.

commutator origin: repeated same-site canonicalization can reduce cubic words
to linear words.

expected derivative contribution: not used by the velocity coefficient-
derivative contraction, but status must be classified rather than assumed.

## 6. Metadata or Structurally Zero Blocks

metadata-only path: `assemble_kjjsj_terms(..., metadata_only=True)` returns
zero arrays and records `normal_form="K3+K2_comm+K1_comm"`.

structurally zero block: virtual LLL derivative is zero because it has no
Wilson-line coefficient factor.
