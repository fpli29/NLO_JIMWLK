"""Direct density-operator diagnostics for tiny NLO current lattices."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .nlo_current_skeleton import NLOCurrentTerms, unflatten_index
from .su3_adjoint import left_perturb


@dataclass(frozen=True)
class DensityOperatorValue:
    """Scalar density-operator value with bookkeeping metadata."""

    value: complex
    by_order: dict
    metadata: dict


def _copy_U(U_fund) -> np.ndarray:
    return np.array(U_fund, copy=True)


def _as_complex_scalar(value) -> complex:
    arr = np.asarray(value)
    if arr.shape != ():
        raise ValueError("density operator callbacks must return scalars")
    return complex(arr)


def _left_perturbed(U_fund, flat_index: int, gens, eps: float) -> np.ndarray:
    n_color = int(gens.shape[0])
    site, color = unflatten_index(flat_index, n_color)
    out = _copy_U(U_fund)
    out[site] = left_perturb(out[site], color, eps, gens)
    return out


def fd_ordered_left_derivative_scalar(
    callback: Callable[[np.ndarray], complex],
    U_fund,
    flat_indices: tuple[int, ...],
    gens,
    eps: float,
) -> complex:
    """Central finite-difference ordered derivative for scalar callbacks."""

    if not flat_indices:
        return _as_complex_scalar(callback(_copy_U(U_fund)))
    outer = int(flat_indices[0])

    def inner(V):
        return fd_ordered_left_derivative_scalar(callback, V, tuple(flat_indices[1:]), gens, eps)

    plus = _left_perturbed(U_fund, outer, gens, eps)
    minus = _left_perturbed(U_fund, outer, gens, -eps)
    return (inner(plus) - inner(minus)) / (2.0 * eps)


def _default_terms_builder(terms: NLOCurrentTerms):
    return lambda _U: terms


def _nonzero_indices(arr: np.ndarray, cutoff: float) -> np.ndarray:
    return np.argwhere(np.abs(arr) > cutoff)


def evaluate_direct_density_operator(
    U,
    W_fn,
    physical_terms: NLOCurrentTerms,
    *,
    gens,
    derivative_backend: str = "finite_difference",
    fd_eps: float = 1.0e-4,
    sector_mask=None,
    terms_builder=None,
    active_outer_indices=None,
    coefficient_cutoff: float = 1.0e-14,
) -> DensityOperatorValue:
    """Evaluate the direct normal-form density operator at ``U``.

    The implemented expression is

    ``-L_A(K1^A W) + 1/2 L_A L_B(K2^{AB} W)
      - 1/6 L_A L_B L_C(K3^{ABC} W)``.

    This is dense and intended only for tiny diagnostic lattices. Optional
    ``active_outer_indices`` restricts the outer \(A\) sum for projected
    convergence checks; inner \(B,C\) sums are left intact.
    """

    if derivative_backend != "finite_difference":
        raise ValueError("only derivative_backend='finite_difference' is implemented")
    if sector_mask is not None:
        raise NotImplementedError("sector_mask is handled by the supplied terms_builder in this layer")

    U0 = _copy_U(U)
    base_terms = physical_terms
    base_terms.validate_shapes()
    dim = base_terms.dim
    terms_at = _default_terms_builder(base_terms) if terms_builder is None else terms_builder
    outer = tuple(range(dim)) if active_outer_indices is None else tuple(int(i) for i in active_outer_indices)

    def W_at(V):
        value = W_fn(V)
        return _as_complex_scalar(value.weight if hasattr(value, "weight") else value)

    first = 0.0 + 0.0j
    for a in outer:
        if abs(base_terms.K1[a]) <= coefficient_cutoff:
            continue

        def term(V, a=a):
            return terms_at(V).K1[a] * W_at(V)

        first -= fd_ordered_left_derivative_scalar(term, U0, (a,), gens, fd_eps)

    second = 0.0 + 0.0j
    for a, b in _nonzero_indices(base_terms.K2, coefficient_cutoff):
        if int(a) not in outer:
            continue

        def term(V, a=int(a), b=int(b)):
            return terms_at(V).K2[a, b] * W_at(V)

        second += 0.5 * fd_ordered_left_derivative_scalar(term, U0, (int(a), int(b)), gens, fd_eps)

    third = 0.0 + 0.0j
    for a, b, c in _nonzero_indices(base_terms.K3, coefficient_cutoff):
        if int(a) not in outer:
            continue

        def term(V, a=int(a), b=int(b), c=int(c)):
            return terms_at(V).K3[a, b, c] * W_at(V)

        third -= (1.0 / 6.0) * fd_ordered_left_derivative_scalar(
            term,
            U0,
            (int(a), int(b), int(c)),
            gens,
            fd_eps,
        )

    value = first + second + third
    return DensityOperatorValue(
        value=value,
        by_order={"K1": first, "K2": second, "K3": third},
        metadata={
            "backend": derivative_backend,
            "fd_eps": float(fd_eps),
            "active_outer_indices": list(outer),
            "coefficient_cutoff": float(coefficient_cutoff),
            "nonproduction_only": True,
            "normalized_by_weight": False,
        },
    )


def evaluate_direct_density_operator_by_sector(
    U,
    W_fn,
    sector_builders: dict,
    *,
    gens,
    fd_eps: float = 1.0e-4,
    active_outer_indices=None,
) -> dict:
    """Return direct density-operator contributions for supplied sectors."""

    out = {}
    for name, builder in sector_builders.items():
        terms = builder(U)
        out[name] = evaluate_direct_density_operator(
            U,
            W_fn,
            terms,
            gens=gens,
            fd_eps=fd_eps,
            terms_builder=builder,
            active_outer_indices=active_outer_indices,
        )
    return out

