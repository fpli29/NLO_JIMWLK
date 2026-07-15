# NLO Current Map Status

The current non-production map uses
\[
\partial_YW=-L_AJ^A_{\rm NLO},
\]
with
\[
J^A_{\rm NLO}
=
K_1^AW
-
\frac12L_B(K_2^{AB}W)
+
\frac16L_BL_C(K_3^{ABC}W).
\]

## Dependency table

| sector | normal-form contribution | density derivatives needed |
|---|---|---|
| \(K_{JSJ}\) | \(K_2\) | score |
| \(K_{JSSJ}\) | \(K_2\) | score + coefficient drift |
| \(K_{q\bar q}\) | \(K_2\) | score + coefficient drift |
| \(K_{JJSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |
| \(K_{JJSSJ}\) | \(K_3+K_{2,\rm comm}+K_{1,\rm comm}\) | score + Hessian-score |

## Coefficient derivative status

- A dense finite-difference coefficient-derivative backend exists for tiny
  diagnostics.
- It computes explicit arrays for \(dK2\), \(LC\_K3\), \(LB\_K3\), and
  \(d2K3\) from coefficient callbacks \(K_2(U)\) and \(K_3(U)\).
- The velocity evaluator can now be run with explicit derivative arrays.
- The omitted-derivative path remains available only for smoke tests and
  records warnings when used.
- The production coefficient-derivative strategy remains unresolved.

## Dipole validation status

- Observable-side analytic dipole generator actions match finite differences
  for representative one-, two-, and three-generator words.
- Direct sector actions pass zero-kernel and linearity tests.
- \(K_{JSJ}\) passes the available KLM Appendix A target with max residual
  `1.6172982426630268e-16` using an Appendix-compatible symmetric zero-diagonal
  synthetic kernel.
- \(K_{JSSJ}\), \(K_{q\bar q}\), \(K_{JJSJ}\), and \(K_{JJSSJ}\) also pass
  the local Appendix A targets in the current dipole validation status.
- The \(q\bar q\) \(z'=z\) subtraction test passes.
- The cubic virtual \(1/3\) sensitivity test passes for \(K_{JJSJ}\) and
  \(K_{JJSSJ}\).

## Important caveats

- This is not a production implementation.
- Physical kernels are integrated only in dense non-production diagnostics.
  Barred/nonsinglet choices remain outside the current implementation.
- Coefficient derivatives still need a robust production implementation strategy.
- Score/Hessian-score estimation is not implemented.
- Dense tensors are for tiny lattice diagnostics only.

## Physical kernel integration status

- Non-production unbarred physical coordinate kernels are now implemented for
  \(K_{JSJ}\), \(K_{JSSJ}\), \(K_{q\bar q}\), \(K_{JJSJ}\), \(K_{JJSSJ}\),
  and \(\widetilde K\).
- \(K_{JSJ}\) requires an explicit `KJSJIntegrationPolicy` because
  `WORKNLO.tex` lines 324--332 contain \(\int_{z'}\widetilde K\) and a scheme
  scale. The diagnostic does not choose a quadrature measure or regulator by
  default.
- Barred/nonsinglet kernel modifications are documented but not implemented.
- Dense physical-kernel arrays expose singular coordinate entries with an
  explicit `singularity_policy`; no policy is a production regulator.
- The physical-kernel adapter can pass implemented kernels into
  `assemble_nlo_current_terms(...)` in metadata-only mode and in a smallest
  dense diagnostic assembly. The full assembly uses KLM-normalized cubic
  coefficients and emits no complex-cast warning.
- The physical-kernel dipole recheck passes \(K_{JSJ}\), \(K_{JSSJ}\),
  \(K_{q\bar q}\), \(K_{JJSJ}\), and \(K_{JJSSJ}\). The earlier
  \(K_{q\bar q}\) mismatch was traced to a compact reduced target diagnostic;
  the recheck now uses the exact WORKNLO trace-product expression.
- Physical-kernel positivity checks are future work. The Pawula toy diagnostic
  does not prove physical NLO JIMWLK positivity or non-positivity.

## Physical density-side closure status

- `src/nlo_current/physical_density_closure.py` compares the direct
  normal-form density operator with \(-L_A(v^A W)\) on tiny lattices.
- `src/nlo_current/test_densities.py` supplies positive differentiable test
  densities and finite-difference score/Hessian-score data.
- The physical projected scan uses `KJSSJ` and `Kqbarq`, active outer index
  `0`, and an explicit diagnostic finite-grid policy.
- Best projected physical residual: `3.5453411040275995e-13` absolute.
- Sparse synthetic cubic checks exercise nonzero Hessian-score dependence,
  raw-cubic normalization failure, and commutator-correction failure modes.
- This remains a dense diagnostic closure check, not a production evolution
  algorithm.

## Analytic coefficient derivative status

- Local analytic Lie-derivative primitives are implemented and calibrated
  against finite differences.
- Analytic \(dK2^A=L_BK_2^{AB}\) is implemented for `KJSJ`, `KJSSJ`, and
  `Kqbarq`.
- `Kqbarq` is decomposed into trace and subtraction derivative contributions
  and the full sum.
- The structured backend `src/nlo_current/physical_coefficient_derivatives.py`
  supports `analytic`, `finite_difference`, `diagnostic`, and explicit
  `hybrid_local_fd`.
- `backend="analytic"` does not silently fall back to FD.
- `KJJSJ` analytic cubic `LC_K3`, `LB_K3`, ordered `d2K3`, and `dK2_comm`
  are implemented and validated against the FD oracle in
  `reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`.
- `KJJSSJ` analytic cubic `LC_K3`, `LB_K3`, and `d2K3` remain pending and
  still raise under `backend="analytic"`.
