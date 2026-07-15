"""Comparator for tiny-lattice physical NLO density-side closure."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .physical_current_divergence import evaluate_current_divergence
from .physical_density_operator import evaluate_direct_density_operator
from .physical_nlo_current import (
    PHYSICAL_SECTORS,
    PhysicalNLOCurrentConfig,
    assemble_physical_terms,
    compute_physical_coefficient_derivatives,
)
from .test_densities import TestDensityResult, evaluate_test_density


@dataclass(frozen=True)
class ClosureResult:
    """Direct/current closure comparison result."""

    direct_value: complex
    current_value: complex
    abs_residual: float
    rel_residual: float
    direct_by_sector: dict
    current_by_sector: dict
    sector_residuals: dict
    metadata: dict


def _density_builder(density, gens, score_eps, hessian_eps):
    if isinstance(density, str):
        return lambda U: evaluate_test_density(
            U,
            gens,
            density,
            score_eps=score_eps,
            hessian_eps=hessian_eps,
        )
    if callable(density):
        return density
    if isinstance(density, TestDensityResult):
        return lambda _U: density
    raise TypeError("density must be a name, callable, or TestDensityResult")


def _relative_residual(direct, current, floor=1.0e-30) -> float:
    return float(abs(direct - current) / (abs(direct) + abs(current) + floor))


def raw_cubic_terms_transform(terms):
    """Diagnostic transform that removes the established cubic normalization."""

    from .nlo_current_skeleton import NLOCurrentTerms

    return NLOCurrentTerms(
        K1=np.array(terms.K1, copy=True),
        K2=np.array(terms.K2, copy=True),
        K3=1.0j * np.array(terms.K3, copy=True),
        metadata={**terms.metadata, "diagnostic_transform": "raw_cubic_like_i_times_K3"},
    )


def compare_physical_density_closure(
    U,
    density,
    *,
    coords,
    gens,
    f,
    physical_policy,
    config: PhysicalNLOCurrentConfig | None = None,
    derivative_backend: str = "finite_difference",
    fd_eps: float = 1.0e-4,
    score_eps: float = 2.0e-5,
    hessian_eps: float = 5.0e-4,
    sector_mask=None,
    active_outer_indices=None,
    include_by_sector: bool = False,
    omit_coefficient_derivatives: bool = False,
    omit_hessian_score: bool = False,
    omit_commutators: bool = False,
    remove_cubic_normalization: bool = False,
) -> ClosureResult:
    """Compare direct density operator and generalized-current divergence."""

    cfg = PhysicalNLOCurrentConfig() if config is None else config
    current_cfg = cfg
    if omit_commutators:
        current_cfg = PhysicalNLOCurrentConfig(
            Nc=cfg.Nc,
            nf=cfg.nf,
            alpha_s=cfg.alpha_s,
            singularity_policy=cfg.singularity_policy,
            eps=cfg.eps,
            include_commutators=False,
            derivative_backend=cfg.derivative_backend,
            fd_eps_first=cfg.fd_eps_first,
            fd_eps_second=cfg.fd_eps_second,
            real_atol=cfg.real_atol,
        )

    selected = tuple(PHYSICAL_SECTORS if sector_mask is None else sector_mask)
    density_at = _density_builder(density, gens, score_eps, hessian_eps)

    def direct_terms_builder(V, sectors=selected, cfg=cfg):
        return assemble_physical_terms(
            V,
            coords,
            gens,
            f,
            integration_policy=physical_policy,
            config=cfg,
            sector_filter=sectors,
        )

    def current_terms_builder(V, sectors=selected, cfg=current_cfg):
        return assemble_physical_terms(
            V,
            coords,
            gens,
            f,
            integration_policy=physical_policy,
            config=cfg,
            sector_filter=sectors,
        )

    def current_derivatives_builder(V, sectors=selected, cfg=current_cfg):
        return compute_physical_coefficient_derivatives(
            V,
            coords,
            gens,
            f,
            integration_policy=physical_policy,
            config=cfg,
            sector_filter=sectors,
            backend=derivative_backend,
        )

    U0 = np.asarray(U)
    base_terms = direct_terms_builder(U0)
    current_terms = current_terms_builder(U0)
    base_derivatives = current_derivatives_builder(U0)
    base_density = density_at(U0)

    transform = raw_cubic_terms_transform if remove_cubic_normalization else None

    direct = evaluate_direct_density_operator(
        U0,
        density_at,
        base_terms,
        gens=gens,
        derivative_backend="finite_difference",
        fd_eps=fd_eps,
        terms_builder=direct_terms_builder,
        active_outer_indices=active_outer_indices,
    )
    current = evaluate_current_divergence(
        U0,
        density_at,
        base_density.score,
        base_density.hessian_score,
        current_terms,
        base_derivatives,
        gens=gens,
        derivative_backend="finite_difference",
        fd_eps=fd_eps,
        terms_builder=current_terms_builder,
        derivatives_builder=current_derivatives_builder,
        density_builder=density_at,
        active_outer_indices=active_outer_indices,
        omit_coefficient_derivatives=omit_coefficient_derivatives,
        omit_hessian_score=omit_hessian_score,
        terms_transform=transform,
    )

    direct_by_sector = {}
    current_by_sector = {}
    sector_residuals = {}
    if include_by_sector:
        for sector in selected:
            sector_result = compare_physical_density_closure(
                U0,
                density_at,
                coords=coords,
                gens=gens,
                f=f,
                physical_policy=physical_policy,
                config=cfg,
                derivative_backend=derivative_backend,
                fd_eps=fd_eps,
                score_eps=score_eps,
                hessian_eps=hessian_eps,
                sector_mask=(sector,),
                active_outer_indices=active_outer_indices,
                include_by_sector=False,
            )
            direct_by_sector[sector] = sector_result.direct_value
            current_by_sector[sector] = sector_result.current_value
            sector_residuals[sector] = {
                "abs": sector_result.abs_residual,
                "rel": sector_result.rel_residual,
            }

    abs_residual = float(abs(direct.value - current.value))
    return ClosureResult(
        direct_value=direct.value,
        current_value=current.value,
        abs_residual=abs_residual,
        rel_residual=_relative_residual(direct.value, current.value),
        direct_by_sector=direct_by_sector,
        current_by_sector=current_by_sector,
        sector_residuals=sector_residuals,
        metadata={
            "density": getattr(base_density, "metadata", {}),
            "direct": direct.metadata,
            "current": current.metadata,
            "direct_by_order": direct.by_order,
            "physical_policy": getattr(physical_policy, "description", None),
            "sector_mask": list(selected),
            "active_outer_indices": (
                None if active_outer_indices is None else [int(i) for i in active_outer_indices]
            ),
            "omit_coefficient_derivatives": bool(omit_coefficient_derivatives),
            "omit_hessian_score": bool(omit_hessian_score),
            "omit_commutators": bool(omit_commutators),
            "remove_cubic_normalization": bool(remove_cubic_normalization),
            "nonproduction_only": True,
            "warning": (
                "Dense finite-difference diagnostic. Active outer indices "
                "produce projected closure checks when not None."
            ),
        },
    )
