# KJJSJ Analytic Cubic Start Status

## Scope

This workflow covers only the \(K_{JJSJ}\) analytic cubic coefficient
derivatives. It does not implement \(K_{JJSSJ}\), production evolution,
score/Hessian-score training, regulator-independence claims, or positivity
claims.

## Baseline Tests

Command:

```text
python3 -m pytest tests/nlo_current -q
```

Result:

```text
143 passed in 61.13s (0:01:01)
```

## Current Analytic Two-Generator Status

The analytic backend is already FD-validated for:

- `KJSJ`: analytic `dK2`;
- `KJSSJ`: analytic `dK2`;
- `Kqbarq`: analytic trace, subtraction, and full `dK2`.

## Current Pending Cubic Status

Before this workflow, both cubic sectors were pending:

- `KJJSJ`: `LC_K3`, `LB_K3`, `d2K3`;
- `KJJSSJ`: `LC_K3`, `LB_K3`, `d2K3`.

This task may only complete `KJJSJ`. `KJJSSJ` must remain pending.

## FD Reference Steps

The prior analytic coefficient derivative validation used FD reference steps
`2e-5` and `1e-5` for two-generator `dK2`. The current KJJSJ cubic workflow
will scan larger second-derivative-stable steps where needed, with the
finite-difference backend preserved as the oracle.

## Cubic Normalization

The physical adapter applies:

```text
raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient
```

Analytic derivatives must act on the KLM-normalized coefficient used by the
normal-form skeleton. The \((-i)\) factor must not be applied again.

## Current Closure Residuals

From the previous analytic validation report:

- two-generator analytic projected closure absolute residual:
  `4.8914478634410496e-11`;
- relative residual: `6.088881138827276e-08`.

KJJSJ analytic closure is not validated at workflow start.

