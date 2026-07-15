# Full Dipole Validation Plan

This workflow validates the observable-side action of each current sector on
the fundamental dipole. It is not a production evolution implementation.

## Observable and target

The validation observable is

\[
s(u,v)=\frac{1}{N_c}\mathrm{tr}[U^\dagger(u)U(v)],
\qquad N_c=3.
\]

The observable-side evolution target is

\[
\frac{d}{dY}s(u,v)=-H_{\rm NLO}s(u,v).
\]

## Kernel choice

For singlet dipole validation the target kernels are unbarred:

\[
K_{JSJ},K_{JSSJ},K_{q\bar q},K_{JJSJ},K_{JJSSJ}.
\]

Barred or nonsinglet kernels are outside this validation workflow.

## Two independent paths

Path A is direct observable-side generator action. The Hamiltonian coefficients
are evaluated on the Wilson-line configuration, and the ordered generators act
on the dipole observable.

Path B is a closed-form KLM Appendix A expression. A sector is fully validated
only when Path A and Path B agree within numerical tolerance.

## Internal fallback

If the exact Appendix A target is unavailable locally, the sector is checked
only for internal consistency:

- direct analytic generator action versus finite-difference action;
- zero-kernel behavior;
- linearity in the kernel;
- symmetry or antisymmetry stress tests;
- known subtraction identities such as the \(q\bar q\) \(z'=z\) vanishing
  identity.

Such sectors are marked

```text
internal-consistency-only, Appendix A target missing
```

and must not be reported as Appendix A passed.
