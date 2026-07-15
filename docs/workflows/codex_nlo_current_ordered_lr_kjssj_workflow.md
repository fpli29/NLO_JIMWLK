# Codex Workflow: NLO JIMWLK Ordered \(J_L A J_R\) Current and \(K_{JSSJ}\) Validation

## Purpose

This workflow is for Codex to continue the NLO generalized probability-current project.

The immediate goal is **not** to implement the full NLO flow.  
The immediate goal is to make the second two-generator term \(K_{JSSJ}\) mathematically and numerically safe.

We need Codex to:

1. Fix the \(K_{JSJ}\) sign convention in the internal notes.
2. Derive and document an ordered \(J_L A[U] J_R\) density-side current lemma.
3. Apply that lemma to \(K_{JSSJ}\).
4. Implement small-lattice algebraic and finite-difference checks.
5. Produce a report deciding whether the \(K_{JSSJ}\) coefficient is symmetric or whether an antisymmetric/commutator drift must be kept.

---

## Background

Kovner--Lublinsky--Mulian use the observable-side convention

\[
\frac{d}{dY}{\cal O}=-H^{\rm JIMWLK}{\cal O}.
\]

Therefore the density-side evolution is

\[
\partial_YW=-H^\dagger W.
\]

For a two-generator observable-side term

\[
H_2=C_H^{AB}[U]L_A L_B,
\]

the density-side diffusion-like tensor satisfies

\[
D^{AB}=-2C_H^{(AB)}.
\]

For an LO-like divergence kernel,

\[
\chi^{AB}=-C_H^{AB}.
\]

This sign has already been fixed for \(K_{JSJ}\).

The corrected \(K_{JSJ}\) contribution is

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
\left(S_x^{bd}-S_z^{bd}\right)
\left(S_y^{cd}-S_z^{cd}\right),
\]

and

\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
\]

For \(K_{JSSJ}\), use the KLM Appendix B identity:

\[
f^{abc}f^{def}J_L^a(x)S_z^{be}S_{z'}^{cf}J_R^d(y)
-
N_cJ_L^a(x)S_z^{ab}J_R^b(y)
\]

\[
=
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right]
J_L^a(x)J_R^b(y).
\]

Then the \(K_{JSSJ}\) observable-side coefficient before converting \(J_R\) is

\[
A_{JSSJ}^{ab}(x,y;U)
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
\left[
S_{z'}^{cf}-S_z^{cf}
\right].
\]

Using

\[
J_R^b(y)=S_y^{hb}J_L^h(y),
\]

the left-basis coefficient is

\[
C_{JSSJ}^{(x,a)(y,h)}
=
A_{JSSJ}^{ab}(x,y;U)S_y^{hb}.
\]

The \(K_{JSSJ}\) term is still second order in charge generators, so it should need the score \(s_A=L_A\log W\), but not Hessian-score \(L_As_B\).  
However, unlike \(K_{JSJ}\), it is not manifestly a symmetric positive square kernel. The antisymmetric part must be checked.

---

## Important Warnings

Do **not** assume

\[
v_{JSSJ}=-\chi_{JSSJ}s.
\]

That may miss coefficient-derivative and commutator drift terms.

Do **not** assume

\[
C_{JSSJ}^{AB}=C_{JSSJ}^{BA}.
\]

Measure it.

Do **not** use barred kernels for dipole-only singlet validation unless explicitly testing the nonsinglet/generalized Hamiltonian.

Use:

- unbarred kernels for dipole/NLO BK validation;
- barred kernels for configuration-level current.

---

## Branch

Create a working branch:

```bash
git checkout -b nlo-current-ordered-lr-kjssj
```

If the repo is not clean, stop and report the dirty files before editing.

---

## Phase 0: Locate Existing Files

Search for existing files related to:

- LO JIMWLK kernel action;
- adjoint Wilson line construction;
- SU(3) generators;
- structure constants \(f^{abc}\);
- score model output convention;
- left Lie derivative convention;
- current derivation notes.

Use:

```bash
find . -maxdepth 4 -type f | sed 's#^\./##' | sort | grep -Ei 'jimwlk|su3|score|current|kernel|nlo|lie|adjoint'
```

Then inspect likely files.

Create a short note:

```text
reports/nlo_current/repo_map.md
```

containing:

- relevant files found;
- functions/classes to reuse;
- whether left Lie derivative is implemented as \(U\to e^{i\epsilon t^a}U\) or another convention;
- whether adjoint Wilson line convention is \(S_A^{ab}=2\mathrm{tr}(t^a S t^b S^\dagger)\).

---

## Phase 1: Add/Update Derivation Notes

Create or update:

```text
docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
```

This document must include the following sections.

### 1. KLM sign convention

Include:

\[
\frac{d}{dY}{\cal O}=-H{\cal O},
\qquad
\partial_YW=-H^\dagger W.
\]

Then state:

\[
D^{AB}=-2C_H^{(AB)},
\qquad
\chi^{AB}=-C_H^{AB}
\]

for a symmetric LO-like two-generator piece.

### 2. Corrected \(K_{JSJ}\) result

Include:

\[
\chi_{JSJ}^{(x,b)(y,c)}
=
-\int_z
\bar K_{JSJ}(x,y;z)
(S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}),
\]

and

\[
v_{JSJ}^A=-\chi_{JSJ}^{AB}s_B.
\]

Also state that using \(K_{JSJ}\to -M/2\) recovers the LO sign.

### 3. Ordered \(J_L A J_R\) current lemma

Derive the density-side current for

\[
H_{LR}[A]=\int_{x,y}A^{ab}(x,y;U)J_L^a(x)J_R^b(y),
\]

with observable evolution

\[
\partial_Y{\cal O}=-H_{LR}{\cal O}.
\]

The exact density-side expression is

\[
(\partial_YW)_{LR}
=
-J_R^b(y)J_L^a(x)
\left[
A^{ab}(x,y;U)W
\right].
\]

Using \(J_R^b(y)=S_y^{hb}L_y^h\), choose a canonical left-divergence form.

A safe form to document is

\[
(\partial_YW)_{LR}
=
-L_y^h
\left\{
S_y^{hb}
L_x^a
\left[
A^{ab}(x,y;U)W
\right]
\right\}
+
\Delta_{\rm div}[A,W],
\]

where \(\Delta_{\rm div}\) is zero if the chosen left/right conversion has no extra divergence term, or otherwise must be explicitly computed.

Then test whether

\[
L_y^hS_y^{hb}=0
\]

under the repo convention.

If this identity is verified, the current component is

\[
J_{LR}^{(y,h)}
=
S_y^{hb}
L_x^a[A^{ab}W],
\]

and the velocity is

\[
v_{LR}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA^{ab}
+
A^{ab}s_x^a
\right].
\]

If the identity is **not** verified, keep the extra term and document it.

### 4. Apply to \(K_{JSSJ}\)

Define

\[
A_{JSSJ}^{ab}(x,y;U)
=
\int_{z,z'}
\bar K_{JSSJ}(x,y;z,z')
f^{adc}f^{bef}
S_z^{de}
(S_{z'}^{cf}-S_z^{cf}).
\]

Then

\[
C_{JSSJ}^{(x,a)(y,h)}
=
A_{JSSJ}^{ab}S_y^{hb}.
\]

If the ordered-block identity is verified, write

\[
v_{JSSJ}^{(y,h)}
=
S_y^{hb}
\left[
L_x^aA_{JSSJ}^{ab}
+
A_{JSSJ}^{ab}s_x^a
\right].
\]

Also include the generic second-order decomposition:

\[
C^{AB}=C^{(AB)}+C^{[AB]},
\]

\[
D^{AB}=-2C^{(AB)}.
\]

State clearly:

\[
K_{JSSJ}\text{ needs }s_A\text{ but not }L_As_B.
\]

---

## Phase 2: Implement Algebra Utilities if Missing

Create or reuse a module such as:

```text
src/nlo_current/su3_adjoint.py
```

or place under the existing package layout.

Required functions:

```python
def su3_generators_fundamental():
    """Return normalized SU(3) fundamental generators t^a with Tr(t^a t^b)=delta^{ab}/2."""

def structure_constants(gens):
    """Return f[a,b,c] from [t^a,t^b]=i f^{abc} t^c."""

def adjoint_from_fundamental(U, gens):
    """Return S_A^{ab}=2 Re Tr(t^a U t^b U^\dagger), or the repo's convention."""

def random_su3(rng):
    """Generate a random SU(3) matrix for tests."""

def left_perturb(U, a, eps, gens):
    """Return exp(i eps t^a) U."""

def finite_diff_left_derivative(func, U_list, site, a, eps=1e-6):
    """Finite-difference left Lie derivative of func with respect to U_site."""
```

If equivalents already exist, do not duplicate them. Add thin wrappers only.

---

## Phase 3: Tests for Conventions and Identities

Create tests in a location consistent with the repo, e.g.

```text
tests/nlo_current/test_su3_adjoint_conventions.py
```

### Test 1: adjoint orthogonality

For random SU(3) matrices \(U\), verify:

\[
S_A S_A^T=I
\]

within tolerance.

### Test 2: adjoint structure-constant invariance

Verify numerically:

\[
f^{abc}S^{bd}S^{ce}=S^{af}f^{fde}.
\]

This is needed for the KLM Appendix B identity.

### Test 3: \(J_R\) to \(J_L\) convention

Numerically verify the repo convention:

\[
J_R^a=S_A^{ba}J_L^b.
\]

Use a simple scalar test function of \(U\), for example

\[
F(U)=\mathrm{Re}\,\mathrm{Tr}(A U)
\]

with random fixed \(A\).

### Test 4: \(L_y^hS_y^{hb}\)

Numerically evaluate

\[
\sum_h L_y^h S_y^{hb}
\]

for random \(U_y\), all \(b\), and report whether it is zero.

This test is critical for the ordered-block current lemma.

Do not silently force this to pass. If it is nonzero, write the measured values to a report and update the derivation accordingly.

---

## Phase 4: Test KLM Appendix B Identity

Create:

```text
tests/nlo_current/test_klm_B5_identity.py
```

Numerically test:

\[
f^{abc}f^{def}S_z^{be}S_{z'}^{cf}
-
N_cS_z^{ad}
=
f^{adc}f^{bef}S_z^{de}
(S_{z'}^{cf}-S_z^{cf})
\]

for random adjoint Wilson lines \(S_z,S_{z'}\).

Acceptance:

```text
relative_error < 1e-10
```

or a tolerance appropriate for double precision.

Also test the subtraction:

\[
z'=z\Rightarrow \text{coefficient}=0.
\]

---

## Phase 5: Implement \(K_{JSJ}\) and \(K_{JSSJ}\) Symbolic Coefficient Builders

Create or update:

```text
src/nlo_current/two_generator_terms.py
```

Add:

```python
def kjsj_chi_from_kernel(S_adj, KJSJ, use_barred=True):
    """
    Build or apply the corrected K_JSJ chi:
    chi^{(x,b)(y,c)} = - sum_z K_JSJ(x,y;z)
        (S_x^{bd}-S_z^{bd})(S_y^{cd}-S_z^{cd}).
    """

def kjssj_A_from_kernel(S_adj, KJSSJ, f):
    """
    Build A_JSSJ^{ab}(x,y):
    A^{ab}(x,y) = sum_{z,z'} K_JSSJ(x,y;z,z')
        f^{adc} f^{bef} S_z^{de} (S_zp^{cf}-S_z^{cf}).
    """

def kjssj_C_left_from_A(A, S_adj):
    """
    Convert A^{ab}(x,y) J_L^a(x) J_R^b(y)
    to left-basis C^{(x,a)(y,h)} = A^{ab}(x,y) S_y^{hb}.
    """

def sym_asym_parts(C):
    """Return C_sym, C_asym under combined index exchange."""
```

For first implementation, small dense tensors are fine. We only need tiny-lattice validation.

Use tiny lattice sizes, e.g. \(N_{\rm site}=3\) or \(4\), not production volumes.

---

## Phase 6: Symmetry Report for \(C_{JSSJ}\)

Create script:

```text
scripts/nlo_current/check_kjssj_symmetry.py
```

It should:

1. Generate random SU(3) Wilson lines for \(N_{\rm site}=3\) or \(4\).
2. Generate a synthetic kernel satisfying KLM-like symmetries:
   \[
   K(x,y;z,z')=K(y,x;z,z')
   \]
   and optionally
   \[
   K(x,y;z,z')=K(x,y;z',z).
   \]
3. Build \(A_{JSSJ}^{ab}\).
4. Convert to \(C^{(x,a)(y,h)}\).
5. Compute
   \[
   r_{\rm asym}=\frac{\|C-C^T\|}{\|C\|}.
   \]
6. Save a report:

```text
reports/nlo_current/kjssj_symmetry_report.md
```

The report must include:

- kernel symmetry used;
- random seed;
- \(r_{\rm asym}\);
- conclusion:
  - if \(r_{\rm asym}\approx 0\), the coefficient is symmetric under tested assumptions;
  - otherwise, keep antisymmetric/commutator drift.

This report is diagnostic, not a proof. It is used to avoid making a wrong simplification.

---

## Phase 7: Ordered \(J_L A J_R\) Current Finite-Difference Check

Create:

```text
tests/nlo_current/test_ordered_lr_current.py
```

Use a tiny finite-dimensional approximation with random SU(3) Wilson lines and a simple coefficient \(A^{ab}(U)\). Test the identity:

\[
(\partial_YW)_{LR}
=
-J_R^bJ_L^a(A^{ab}W)
\]

against the proposed divergence-current form.

Suggested test object:

\[
W_{\rm toy}(U)=\exp(\lambda\,\mathrm{ReTr}(Q U_0))
\]

so that the score can be computed by finite differences.

Do not need physical \(W\). This is a convention/order test.

Acceptance:

- finite-difference residual decreases with \(\epsilon\);
- report residuals for \(\epsilon=10^{-4},10^{-5},10^{-6}\);
- assert only a loose tolerance initially, e.g. \(10^{-4}\) to \(10^{-5}\).

If the test fails, do not patch around it. Update the lemma and report the extra term.

---

## Phase 8: Dipole Validation Skeleton

Create a non-production script:

```text
scripts/nlo_current/validate_dipole_two_generator_terms.py
```

This may initially be a scaffold with TODOs, but it must clearly define the target checks.

Targets:

### \(K_{JSJ}\)

Verify KLM Appendix A structure:

\[
-H_{JSJ}s(u,v)
=
2N_c\int_z K_{JSJ}(u,v;z)
[s(u,z)s(z,v)-s(u,v)].
\]

With LO replacement

\[
K_{JSJ}\to -\frac12M,
\]

recover LO BK sign.

### \(K_{JSSJ}\)

Target KLM Appendix A expression:

\[
-H_{JSSJ}s(u,v)
=
-\int_{z,z'}
K_{JSSJ}(u,v;z,z')
\left[
N_c^2s(u,z')s(z',z)s(z,v)
-
\frac1{N_c}
\mathrm{tr}
\left(
S^\dagger(u)S(z)S^\dagger(z')S(v)S^\dagger(z)S(z')
\right)
\right],
\]

with the precise normalization checked against the KLM appendix.

At this phase it is acceptable to mark this as TODO if the dipole observable utilities are missing. But the script should contain the formula and implementation hooks.

---

## Phase 9: Acceptance Criteria

Codex should stop when all of the following are done:

1. `docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md` exists and contains:
   - sign convention;
   - corrected \(K_{JSJ}\);
   - ordered \(J_L A J_R\) lemma;
   - \(K_{JSSJ}\) coefficient;
   - generic second-order current;
   - statement that \(K_{JSSJ}\) needs score but not Hessian-score.

2. Tests pass:
   - adjoint orthogonality;
   - structure-constant invariance;
   - \(J_R\leftrightarrow J_L\) relation;
   - Appendix B identity;
   - \(z'=z\) subtraction.

3. A report exists:
   ```text
   reports/nlo_current/kjssj_symmetry_report.md
   ```

4. The report states clearly whether \(C_{JSSJ}\) appears symmetric in the tested setting.

5. No production NLO flow is added yet.

---

## Phase 10: Suggested Commit Messages

Use small commits:

```bash
git add docs/nlo_current/KJSJ_signfix_KJSSJ_ordered_current.md
git commit -m "docs: derive ordered NLO two-generator current"
```

```bash
git add src/nlo_current tests/nlo_current
git commit -m "test: add SU3 adjoint and KLM identity checks"
```

```bash
git add scripts/nlo_current reports/nlo_current
git commit -m "analysis: report KJSSJ coefficient symmetry"
```

---

## Final Codex Response Required

At the end, Codex should summarize:

1. Which files were added/modified.
2. Which tests were run and their results.
3. Whether \(L_y^hS_y^{hb}=0\) holds under the repo convention.
4. The measured \(K_{JSSJ}\) asymmetry ratio.
5. Whether the next step should be:
   - simplify \(K_{JSSJ}\) to a pure score-current, or
   - keep a generic second-order current with commutator drift.
6. Any unresolved sign/convention issues.

Do not claim the NLO current is complete. This workflow only completes the ordered two-generator block and the \(K_{JSSJ}\) safety checks.
