"""Current-divergence diagnostics for physical NLO density closure."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nlo_velocity_evaluator import evaluate_velocity_from_terms
from .physical_density_operator import fd_ordered_left_derivative_scalar


@dataclass(frozen=True)
class CurrentDivergenceValue:
    """Scalar current-divergence value with diagnostic metadata."""

    value: complex
    metadata: dict


def _weight(density_result) -> complex:
    return complex(density_result.weight if hasattr(density_result, "weight") else density_result)


def _zero_derivatives(dim: int) -> dict:
    return {
        "dK2": np.zeros(dim),
        "dK3_first": {"LC_K3_ABC": np.zeros((dim, dim)), "LB_K3_ABC": np.zeros((dim, dim))},
        "d2K3": np.zeros(dim),
    }


def evaluate_current_divergence(
    U,
    W_fn,
    score,
    hessian_score,
    physical_terms,
    coefficient_derivatives,
    *,
    gens,
    derivative_backend: str = "finite_difference",
    fd_eps: float = 1.0e-4,
    sector_mask=None,
    terms_builder=None,
    derivatives_builder=None,
    density_builder=None,
    active_outer_indices=None,
    omit_coefficient_derivatives: bool = False,
    omit_hessian_score: bool = False,
    terms_transform=None,
) -> CurrentDivergenceValue:
    """Evaluate ``-L_A(v^A W)`` with recomputation under perturbations."""

    if derivative_backend != "finite_difference":
        raise ValueError("only derivative_backend='finite_difference' is implemented")
    if sector_mask is not None:
        raise NotImplementedError("sector_mask is handled by the supplied builders in this layer")

    physical_terms.validate_shapes()
    dim = physical_terms.dim
    outer = tuple(range(dim)) if active_outer_indices is None else tuple(int(i) for i in active_outer_indices)

    def terms_at(V):
        terms = physical_terms if terms_builder is None else terms_builder(V)
        return terms if terms_transform is None else terms_transform(terms)

    def derivs_at(V, terms):
        if omit_coefficient_derivatives:
            return _zero_derivatives(terms.dim)
        if derivatives_builder is None:
            return coefficient_derivatives
        return derivatives_builder(V)

    def density_at(V):
        if density_builder is not None:
            return density_builder(V)
        return W_fn(V)

    def current_component(V, a):
        density = density_at(V)
        terms = terms_at(V)
        derivs = derivs_at(V, terms)
        local_score = np.asarray(density.score if hasattr(density, "score") else score)
        if omit_hessian_score:
            local_hessian = np.zeros((terms.dim, terms.dim), dtype=np.result_type(local_score, float))
        else:
            local_hessian = np.asarray(
                density.hessian_score if hasattr(density, "hessian_score") else hessian_score
            )
        velocity = evaluate_velocity_from_terms(
            terms,
            local_score,
            local_hessian,
            dK2=derivs.get("dK2"),
            dK3_first=derivs.get("dK3_first"),
            d2K3=derivs.get("d2K3"),
        )["velocity"]
        return velocity[a] * _weight(density)

    value = 0.0 + 0.0j
    for a in outer:
        value -= fd_ordered_left_derivative_scalar(
            lambda V, a=a: current_component(V, a),
            U,
            (a,),
            gens,
            fd_eps,
        )

    return CurrentDivergenceValue(
        value=value,
        metadata={
            "backend": derivative_backend,
            "fd_eps": float(fd_eps),
            "active_outer_indices": list(outer),
            "omit_coefficient_derivatives": bool(omit_coefficient_derivatives),
            "omit_hessian_score": bool(omit_hessian_score),
            "terms_transform": getattr(terms_transform, "__name__", None),
            "nonproduction_only": True,
            "velocity_recomputed_under_outer_derivative": True,
        },
    )


def evaluate_current_divergence_by_sector(
    U,
    W_fn,
    sector_builders: dict,
    derivative_builders: dict,
    *,
    gens,
    density_builder,
    fd_eps: float = 1.0e-4,
    active_outer_indices=None,
) -> dict:
    """Return current-divergence contributions for supplied sectors."""

    out = {}
    for name, builder in sector_builders.items():
        terms = builder(U)
        density = density_builder(U)
        derivs = derivative_builders[name](U)
        out[name] = evaluate_current_divergence(
            U,
            W_fn,
            density.score,
            density.hessian_score,
            terms,
            derivs,
            gens=gens,
            terms_builder=builder,
            derivatives_builder=derivative_builders[name],
            density_builder=density_builder,
            fd_eps=fd_eps,
            active_outer_indices=active_outer_indices,
        )
    return out

