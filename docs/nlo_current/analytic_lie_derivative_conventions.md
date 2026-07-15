# Analytic Lie-Derivative Conventions

## Left Generator

The code convention is

\[
L_x^aF(U)
=
\left.
\frac{d}{d\epsilon}
F(e^{i\epsilon t^a}U_x)
\right|_{\epsilon=0}.
\]

Therefore

\[
L_x^aU_y
=
i\,\delta_{xy}\,t^aU_y,
\]

and

\[
L_x^aU_y^\dagger
=
-i\,\delta_{xy}\,U_y^\dagger t^a.
\]

The fundamental generators are \(t^a=\lambda^a/2\), normalized by
\(\mathrm{Tr}(t^a t^b)=\delta^{ab}/2\), and
\([t^a,t^b]=i f^{abc}t^c\).

## Adjoint Wilson Line

The code defines

\[
S_A^{ab}(U)=2\,\mathrm{Re\,Tr}(t^aUt^bU^\dagger).
\]

Using the left perturbation above,

\[
L^h S_A^{ab}
=
f^{hac}S_A^{cb}.
\]

Thus the left derivative acts on the first adjoint index in the project
convention. This sign and index orientation are checked against finite
differences in `tests/nlo_current/test_analytic_lie_derivatives.py`; they are
not imported from memory.

The ordered second derivative is

\[
L^gL^hS_A^{ab}
=
f^{hac}f^{gcd}S_A^{db},
\]

when both derivatives act at the target site, and zero otherwise.

## Left/Right Conversion

The existing convention tests verify

\[
J_R^a = S_A^{ba}J_L^b.
\]

The analytic coefficient derivative modules differentiate the left-basis
normal-form tensors actually used by `nlo_current_skeleton.py`.

## Commutators and Ordering

Same-site left derivatives obey

\[
[L_x^a,L_x^b]=f^{abc}L_x^c.
\]

Distinct-site left derivatives commute. Ordered second derivatives are
preserved exactly; no Hessian or coefficient derivative tensor is silently
symmetrized.

## Cubic Normalization

Raw physical \(K_{JJSJ}\) and \(K_{JJSSJ}\) kernels contain the explicit
WORKNLO complex convention. The physical adapter supplies the normal-form
skeleton with

```text
KLM-normalized cubic coefficient = (-1j) * raw physical cubic kernel
```

Analytic derivatives must differentiate the KLM-normalized coefficient used by
the skeleton and must not apply a duplicate \((-i)\) factor.

