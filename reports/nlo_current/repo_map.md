# Repository Map for NLO Current Checks

No Git metadata is present in this workspace, so this workflow is running in
local no-git mode.

## Relevant files found

- `NLO_JIMWLK_current_worklog_KJSJ.md`: prior worklog with KLM sign, left/right
  relation, K_JSJ derivation, and K_JSSJ planning.
- `NLO_JIMWLK_current_KJSJ_signfix_KJSSJ_derivation.md`: updated derivation note
  containing the corrected K_JSJ sign and preliminary K_JSSJ coefficient.
- `WORKNLO.tex`: KLM paper source text including Appendix B color identities and
  the left/right generator relation.
- `codex_nlo_current_ordered_lr_kjssj_workflow.md`: controlling workflow.

No existing `src`, `tests`, `docs`, `scripts`, or `reports` implementation tree
was present before this workflow run.

## Functions and classes to reuse

No reusable Python functions or classes were found. The validation code added in
`src/nlo_current/` is self-contained and dense, intended only for tiny-lattice
tests.

## Left Lie derivative convention

The existing notes state the paper relation
\[
J_L^a(x)=S_A^{ab}(x)J_R^b(x),
\qquad
J_R^a(x)=S_A^{ba}(x)J_L^b(x).
\]

No code convention existed. The new small-lattice utilities define
\[
L^aF(U)=dF(e^{i\epsilon t^a}U)/d\epsilon|_{\epsilon=0},
\]
and
\[
R^aF(U)=dF(Ue^{i\epsilon t^a})/d\epsilon|_{\epsilon=0}.
\]

The tests verify \(J_R^a=S_A^{ba}J_L^b\) for this convention.

## Adjoint Wilson line convention

The existing notes and `WORKNLO.tex` use the KLM adjoint convention
\[
S_A^{ab}=2\,{\rm tr}(t^a S t^b S^\dagger).
\]

The new utilities use the real numerical form
\[
S_A^{ab}=2\,{\rm Re}\,{\rm tr}(t^a U t^b U^\dagger),
\]
with Hermitian fundamental generators normalized by
\[
{\rm tr}(t^a t^b)=\delta^{ab}/2.
\]

