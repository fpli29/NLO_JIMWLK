# KJJSJ Analytic Cubic Failure Modes

## Scope

This report lists failure modes checked or guarded for the \(K_{JJSJ}\)
analytic cubic coefficient derivatives. It is specific to the dense
non-production diagnostic backend.

## Wrong Adjoint-Index Orientation

The implementation depends on the calibrated rule
\[
L^hS_A^{ab}=f^{hac}S_A^{cb}.
\]
Acting on the wrong adjoint index changes the LLR/LRR derivative tensors and
produces FD residuals in `LC_K3`, `LB_K3`, and `d2K3`.

## Wrong Structure-Constant Sign

Both the KJJSJ coefficient blocks and same-site commutator corrections carry
explicit \(f^{abc}\). Flipping the sign changes the LLR/LRR block residuals and
the `dK2_comm` comparison.

## LLR/LRR Index Swap

The validated blocks use
\[
A_{LLR}^{dea}(x,y,w)=\int_zK_{JJSJ}(w;x,y;z)f^{bde}S_z^{ba},
\]
\[
B_{LRR}^{ade}(w,x,y)=-\int_zK_{JJSJ}(w;x,y;z)f^{bde}S_z^{ab}.
\]
Swapping \((x,y,w)\) or the adjoint indices can pass shape checks while failing
the per-block FD residuals.

## Reversed Derivative Order

The ordered derivative \(L_BL_CK_3^{ABC}\) is not symmetrized. Same-site tests
verify that reversing the derivative order changes the result, while
distinct-site derivatives commute when local support makes the commutator
zero.

## Duplicate Cubic Normalization

The physical adapter already applies

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized coefficient
```

The analytic derivative code differentiates that normalized coefficient. A
second \((-i)\) factor, or omitting the established normalization upstream,
changes the complex character and fails the cubic convention tests.

## Omitted Virtual Block

The virtual LLL coefficient is Wilson-line independent. In the current
canonicalized diagnostic, the virtual K3 derivative blocks are structurally
zero, while the virtual path still exercises canonicalization and numerical
zero handling. Treating structural zeros through relative residual alone can
mislead; absolute residuals are reported.

## Omitted Commutator Correction

The KJJSJ normal form contains nonzero quadratic commutator corrections.
`dK2_comm` is differentiated analytically from the canonicalized tensors and
compared against FD. The current `K1_comm` status is classified rather than
assumed.

## Hidden FD Fallback

`backend="analytic"` is tested with the global FD function monkeypatched to
raise. KJJSJ succeeds without fallback. Requests including `KJJSSJ` still raise
`NotImplementedError`.

## Complex-to-Real Cast

The KJJSJ analytic path preserves complex dtype through intermediate arrays
and only reports real-if-close values after validation. Tests assert that no
`ComplexWarning` is emitted and that raw-cubic diagnostics do not silently
discard imaginary information.

## Structurally Zero Relative-Error Instability

The virtual and physical expected-zero checks can have tiny FD norms. Absolute
residuals, expected-zero tags, and block labels are therefore reported
alongside relative residuals.
