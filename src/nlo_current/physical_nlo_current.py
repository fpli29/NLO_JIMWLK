"""End-to-end non-production physical NLO generalized-current assembly."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import numpy as np

from .coefficient_derivatives import compute_all_coefficient_derivatives_fd, compute_dK2_fd
from .nlo_current_skeleton import NLOCurrentTerms, assemble_nlo_current_terms
from .nlo_velocity_evaluator import evaluate_velocity_from_terms
from .physical_kernel_adapter import physical_kernels_for_skeleton
from .su3_adjoint import adjoint_from_fundamental


PHYSICAL_SECTORS = ("KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ")


@dataclass(frozen=True)
class PhysicalNLOCurrentConfig:
    """Configuration for tiny dense physical-current diagnostics."""

    Nc: int = 3
    nf: int = 0
    alpha_s: float = 1.0
    singularity_policy: str = "raise"
    eps: float | None = None
    include_commutators: bool = True
    derivative_backend: str = "finite_difference"
    fd_eps_first: float = 1.0e-5
    fd_eps_second: float = 1.0e-4
    real_atol: float = 1.0e-10


def _S_builder(gens):
    return lambda U_fund: np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])


def _validate_sector_filter(sector_filter):
    if sector_filter is None:
        return tuple(PHYSICAL_SECTORS)
    selected = tuple(sector_filter)
    unknown = sorted(set(selected) - set(PHYSICAL_SECTORS))
    if unknown:
        raise ValueError(f"unknown physical sectors: {unknown}")
    return selected


def _filter_kernels(kernels: dict, sector_filter) -> dict:
    selected = _validate_sector_filter(sector_filter)
    out = {key: kernels[key] for key in selected}
    metadata = deepcopy(kernels.get("metadata", {}))
    metadata["sector_filter"] = list(selected)
    out["metadata"] = metadata
    return out


def _max_imag(arr) -> float:
    arr = np.asarray(arr)
    return float(np.max(np.abs(np.imag(arr)))) if arr.size else 0.0


def _assembly_metadata(terms: NLOCurrentTerms, kernels: dict, config: PhysicalNLOCurrentConfig) -> dict:
    metadata = terms.metadata
    metadata["physical_nlo_current"] = {
        "nonproduction_only": True,
        "kernel_origin": "unbarred physical KLM kernels via physical_kernel_adapter",
        "kernel_type": "unbarred singlet",
        "adapter_metadata": deepcopy(kernels.get("metadata", {})),
        "sector_labels": list(metadata.get("sectors", {}).keys()),
        "cubic_convention": (
            "raw physical cubic kernel -> (-1j) -> KLM-normalized real coefficient"
        ),
        "include_commutators": bool(config.include_commutators),
        "dtypes": {
            "K1": str(terms.K1.dtype),
            "K2": str(terms.K2.dtype),
            "K3": str(terms.K3.dtype),
        },
        "max_imag": {
            "K1": _max_imag(terms.K1),
            "K2": _max_imag(terms.K2),
            "K3": _max_imag(terms.K3),
        },
    }
    return metadata


def build_physical_kernels(
    coords,
    *,
    integration_policy,
    config: PhysicalNLOCurrentConfig | None = None,
    sector_filter=None,
) -> dict:
    """Build KLM-normalized physical kernels for the dense skeleton."""

    cfg = PhysicalNLOCurrentConfig() if config is None else config
    kernels = physical_kernels_for_skeleton(
        coords,
        Nc=cfg.Nc,
        nf=cfg.nf,
        alpha_s=cfg.alpha_s,
        singularity_policy=cfg.singularity_policy,
        eps=cfg.eps,
        integration_policy=integration_policy,
    )
    return _filter_kernels(kernels, sector_filter)


def assemble_physical_terms(
    U_fund,
    coords,
    gens,
    f,
    *,
    integration_policy,
    config: PhysicalNLOCurrentConfig | None = None,
    sector_filter=None,
    S_adj=None,
    metadata_only: bool = False,
) -> NLOCurrentTerms:
    """Assemble dense physical K1/K2/K3 terms for tiny diagnostics."""

    cfg = PhysicalNLOCurrentConfig() if config is None else config
    U_fund = np.asarray(U_fund)
    S_adj = _S_builder(gens)(U_fund) if S_adj is None else np.asarray(S_adj)
    kernels = build_physical_kernels(
        coords,
        integration_policy=integration_policy,
        config=cfg,
        sector_filter=sector_filter,
    )
    terms = assemble_nlo_current_terms(
        U_fund,
        S_adj,
        kernels,
        gens,
        f,
        include_commutators=cfg.include_commutators,
        metadata_only=metadata_only,
    )
    terms.metadata = _assembly_metadata(terms, kernels, cfg)
    return terms


def assemble_physical_K1(*args, **kwargs):
    """Return the assembled physical K1 vector."""

    return assemble_physical_terms(*args, **kwargs).K1


def assemble_physical_K2(*args, **kwargs):
    """Return the assembled physical K2 matrix."""

    return assemble_physical_terms(*args, **kwargs).K2


def assemble_physical_K3(*args, **kwargs):
    """Return the assembled physical K3 tensor."""

    return assemble_physical_terms(*args, **kwargs).K3


def physical_coefficient_callbacks(
    coords,
    gens,
    f,
    *,
    integration_policy,
    config: PhysicalNLOCurrentConfig | None = None,
    sector_filter=None,
):
    """Return dense K2 and K3 callbacks for coefficient-derivative diagnostics."""

    cfg = PhysicalNLOCurrentConfig() if config is None else config

    def _terms(U_fund, S_adj):
        return assemble_physical_terms(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=integration_policy,
            config=cfg,
            sector_filter=sector_filter,
            S_adj=S_adj,
        )

    return (
        lambda U_fund, S_adj: _terms(U_fund, S_adj).K2,
        lambda U_fund, S_adj: _terms(U_fund, S_adj).K3,
    )


def compute_physical_coefficient_derivatives(
    U_fund,
    coords,
    gens,
    f,
    *,
    integration_policy,
    config: PhysicalNLOCurrentConfig | None = None,
    sector_filter=None,
    backend: str | None = None,
) -> dict:
    """Compute coefficient derivatives for the physical dense current.

    Supported backends:
    - ``finite_difference``: dense central finite differences.
    - ``diagnostic``: alias for ``finite_difference``.
    - ``analytic``: proven local analytic dK2 for completed two-generator
      sectors and analytic KJJSJ cubic derivative contractions; raises for
      pending KJJSSJ.
    - ``hybrid_local_fd``: explicit mixed diagnostic path.
    """

    cfg = PhysicalNLOCurrentConfig() if config is None else config
    selected_backend = cfg.derivative_backend if backend is None else backend
    if selected_backend in {"analytic", "hybrid_local_fd"}:
        from .physical_coefficient_derivatives import (
            compute_physical_coefficient_derivatives as compute_structured_derivatives,
        )

        return compute_structured_derivatives(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=integration_policy,
            config=cfg,
            backend=selected_backend,
            sector_filter=sector_filter,
        ).as_legacy_dict()
    if selected_backend not in {"finite_difference", "diagnostic"}:
        raise ValueError(
            "derivative backend must be 'analytic', 'hybrid_local_fd', 'finite_difference', or 'diagnostic'"
        )

    selected_sectors = _validate_sector_filter(sector_filter)
    K2_callback, K3_callback = physical_coefficient_callbacks(
        coords,
        gens,
        f,
        integration_policy=integration_policy,
        config=cfg,
        sector_filter=sector_filter,
    )
    U_arr = np.asarray(U_fund)
    S_builder = _S_builder(gens)
    dim = int(len(U_arr)) * int(gens.shape[0])
    if not ({"KJJSJ", "KJJSSJ"} & set(selected_sectors)):
        derivatives = {
            "dK2": compute_dK2_fd(
                K2_callback,
                U_arr,
                S_builder,
                gens,
                eps=cfg.fd_eps_first,
            ),
            "dK3_first": {
                "LC_K3_ABC": np.zeros((dim, dim)),
                "LB_K3_ABC": np.zeros((dim, dim)),
            },
            "d2K3": np.zeros(dim),
        }
        computed = ["dK2", "zero structural dK3_first", "zero structural d2K3"]
    else:
        derivatives = compute_all_coefficient_derivatives_fd(
            K2_callback,
            K3_callback,
            U_arr,
            S_builder,
            gens,
            eps_first=cfg.fd_eps_first,
            eps_second=cfg.fd_eps_second,
        )
        computed = ["dK2", "dK3_first", "d2K3"]
    derivatives["metadata"] = {
        "backend": selected_backend,
        "effective_backend": "finite_difference",
        "nonproduction_only": True,
        "sector_filter": list(selected_sectors),
        "computed": computed,
        "fd_eps_first": cfg.fd_eps_first,
        "fd_eps_second": cfg.fd_eps_second,
    }
    return derivatives


def evaluate_physical_nlo_velocity(
    U_fund,
    coords,
    gens,
    f,
    score,
    hessian_score,
    *,
    integration_policy,
    config: PhysicalNLOCurrentConfig | None = None,
    sector_filter=None,
    derivative_backend: str | None = None,
) -> dict:
    """Evaluate the dense non-production physical generalized-current velocity."""

    cfg = PhysicalNLOCurrentConfig() if config is None else config
    terms = assemble_physical_terms(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=integration_policy,
        config=cfg,
        sector_filter=sector_filter,
    )
    derivatives = compute_physical_coefficient_derivatives(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=integration_policy,
        config=cfg,
        sector_filter=sector_filter,
        backend=derivative_backend,
    )
    result = evaluate_velocity_from_terms(
        terms,
        score,
        hessian_score,
        dK2=derivatives["dK2"],
        dK3_first=derivatives["dK3_first"],
        d2K3=derivatives["d2K3"],
    )
    result["terms"] = terms
    result["derivatives"] = derivatives
    result["diagnostics"]["physical_nlo_current"] = terms.metadata["physical_nlo_current"]
    result["diagnostics"]["derivative_backend"] = derivatives["metadata"]
    return result
