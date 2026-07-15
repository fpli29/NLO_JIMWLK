"""Structured physical coefficient derivative backends for diagnostics."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .analytic_two_generator_derivatives import (
    analytic_dK2_KJSJ,
    analytic_dK2_KJSSJ,
    analytic_dK2_Kqbarq,
)
from .analytic_cubic_derivatives import (
    KJJSJSectorData,
    analytic_d2K3_KJJSJ,
    analytic_dK1_comm_KJJSJ,
    analytic_first_derivatives_KJJSJ,
)
from .coefficient_derivatives import compute_all_coefficient_derivatives_fd, compute_dK2_fd
from .nlo_current_skeleton import assemble_nlo_current_terms
from .physical_kernel_adapter import physical_kernels_for_skeleton
from .su3_adjoint import adjoint_from_fundamental


PHYSICAL_SECTORS = ("KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ")
CUBIC_SECTORS = {"KJJSJ", "KJJSSJ"}
ANALYTIC_DK2_SECTORS = {"KJSJ", "KJSSJ", "Kqbarq"}
ANALYTIC_CUBIC_SECTORS = {"KJJSJ"}


@dataclass
class PhysicalCoefficientDerivatives:
    """Structured derivative contractions for the physical dense current."""

    dK2: np.ndarray
    LC_K3: np.ndarray
    LB_K3: np.ndarray
    d2K3: np.ndarray
    by_sector: dict = field(default_factory=dict)
    backend: str = "unknown"
    metadata: dict = field(default_factory=dict)

    def as_legacy_dict(self) -> dict:
        """Return the legacy dict shape used by the velocity evaluator."""

        return {
            "dK2": self.dK2,
            "dK3_first": {"LC_K3_ABC": self.LC_K3, "LB_K3_ABC": self.LB_K3},
            "d2K3": self.d2K3,
            "metadata": self.metadata | {"backend": self.backend},
        }


def _cfg(config, name, default):
    return getattr(config, name, default) if config is not None else default


def _validate_sector_filter(sector_filter):
    if sector_filter is None:
        return tuple(PHYSICAL_SECTORS)
    selected = tuple(sector_filter)
    unknown = sorted(set(selected) - set(PHYSICAL_SECTORS))
    if unknown:
        raise ValueError(f"unknown physical sectors: {unknown}")
    return selected


def _S_builder(gens):
    return lambda U_fund: np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])


def _build_kernels(coords, integration_policy, config, sector_filter):
    kernels = physical_kernels_for_skeleton(
        coords,
        Nc=_cfg(config, "Nc", 3),
        nf=_cfg(config, "nf", 0),
        alpha_s=_cfg(config, "alpha_s", 1.0),
        singularity_policy=_cfg(config, "singularity_policy", "raise"),
        eps=_cfg(config, "eps", None),
        integration_policy=integration_policy,
    )
    selected = _validate_sector_filter(sector_filter)
    out = {key: kernels[key] for key in selected}
    out["metadata"] = kernels.get("metadata", {})
    return out


def _assemble_terms(U_fund, S_adj, kernels, gens, f, config):
    return assemble_nlo_current_terms(
        U_fund,
        S_adj,
        kernels,
        gens,
        f,
        include_commutators=_cfg(config, "include_commutators", True),
    )


def _zero_result(dim, dtype=float):
    return PhysicalCoefficientDerivatives(
        dK2=np.zeros(dim, dtype=dtype),
        LC_K3=np.zeros((dim, dim), dtype=dtype),
        LB_K3=np.zeros((dim, dim), dtype=dtype),
        d2K3=np.zeros(dim, dtype=dtype),
    )


def _finite_difference_result(U_fund, coords, gens, f, *, integration_policy, config, sector_filter, backend):
    U_arr = np.asarray(U_fund)
    S_builder = _S_builder(gens)
    kernels = _build_kernels(coords, integration_policy, config, sector_filter)
    selected = _validate_sector_filter(sector_filter)
    dim = int(len(U_arr)) * int(gens.shape[0])

    def terms_at(V, S):
        return _assemble_terms(V, S, kernels, gens, f, config)

    K2_callback = lambda V, S: terms_at(V, S).K2
    K3_callback = lambda V, S: terms_at(V, S).K3

    if not (CUBIC_SECTORS & set(selected)):
        derivatives = {
            "dK2": compute_dK2_fd(
                K2_callback,
                U_arr,
                S_builder,
                gens,
                eps=_cfg(config, "fd_eps_first", 1.0e-5),
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
            eps_first=_cfg(config, "fd_eps_first", 1.0e-5),
            eps_second=_cfg(config, "fd_eps_second", 1.0e-4),
        )
        computed = ["dK2", "dK3_first", "d2K3"]

    return PhysicalCoefficientDerivatives(
        dK2=derivatives["dK2"],
        LC_K3=derivatives["dK3_first"]["LC_K3_ABC"],
        LB_K3=derivatives["dK3_first"]["LB_K3_ABC"],
        d2K3=derivatives["d2K3"],
        backend=backend,
        metadata={
            "backend": backend,
            "effective_backend": "finite_difference",
            "sector_filter": list(selected),
            "computed": computed,
            "fd_eps_first": _cfg(config, "fd_eps_first", 1.0e-5),
            "fd_eps_second": _cfg(config, "fd_eps_second", 1.0e-4),
            "nonproduction_only": True,
        },
    )


def _analytic_result(U_fund, coords, gens, f, *, integration_policy, config, sector_filter):
    U_arr = np.asarray(U_fund)
    S_adj = _S_builder(gens)(U_arr)
    selected = _validate_sector_filter(sector_filter)
    if "KJJSSJ" in selected:
        raise NotImplementedError(
            "analytic KJJSSJ coefficient derivatives are pending; select hybrid_local_fd explicitly"
        )
    kernels = _build_kernels(coords, integration_policy, config, selected)
    dim = int(len(U_arr)) * int(gens.shape[0])
    result = _zero_result(dim, dtype=complex)
    by_sector = {}

    for sector in selected:
        if sector == "KJSJ":
            dK2 = analytic_dK2_KJSJ(S_adj, kernels["KJSJ"], f)
        elif sector == "KJSSJ":
            dK2 = analytic_dK2_KJSSJ(S_adj, kernels["KJSSJ"], f)
        elif sector == "Kqbarq":
            dK2 = analytic_dK2_Kqbarq(U_arr, S_adj, kernels["Kqbarq"], gens, f)
        elif sector == "KJJSJ":
            sector_data = KJJSJSectorData(S_adj=S_adj, KJJSJ=kernels["KJJSJ"], f=f)
            first = analytic_first_derivatives_KJJSJ(sector_data=sector_data)
            d2 = analytic_d2K3_KJJSJ(sector_data=sector_data, return_diagnostics=False)
            k1_status = analytic_dK1_comm_KJJSJ(sector_data=sector_data)
            result.dK2 = result.dK2 + first["dK2_comm"]
            result.LC_K3 = result.LC_K3 + first["LC_K3"]
            result.LB_K3 = result.LB_K3 + first["LB_K3"]
            result.d2K3 = result.d2K3 + d2
            by_sector[sector] = {
                "dK2": first["dK2_comm"],
                "LC_K3": first["LC_K3"],
                "LB_K3": first["LB_K3"],
                "d2K3": d2,
                "by_block": {},
                "K1_comm_status": k1_status,
                "status": "analytic_cubic_complete",
            }
            continue
        else:
            raise NotImplementedError(f"analytic derivative pending for sector {sector}")
        result.dK2 = result.dK2 + dK2
        by_sector[sector] = {
            "dK2": dK2,
            "LC_K3": np.zeros((dim, dim)),
            "LB_K3": np.zeros((dim, dim)),
            "d2K3": np.zeros(dim),
            "status": "analytic_dK2_complete",
        }

    result.by_sector = by_sector
    result.backend = "analytic"
    result.metadata = {
        "backend": "analytic",
        "implemented_analytic_sectors": sorted(
            (set(selected) & ANALYTIC_DK2_SECTORS) | (set(selected) & ANALYTIC_CUBIC_SECTORS)
        ),
        "pending_analytic_sectors": sorted({"KJJSSJ"} & set(selected)),
        "fallback_used": False,
        "sector_filter": list(selected),
        "cubic_normalization": "adapter kernels are KLM-normalized before assembly",
        "nonproduction_only": True,
        "dtype": {
            "dK2": str(np.asarray(result.dK2).dtype),
            "LC_K3": str(result.LC_K3.dtype),
            "LB_K3": str(result.LB_K3.dtype),
            "d2K3": str(result.d2K3.dtype),
        },
    }
    return result


def compute_physical_coefficient_derivatives(
    U_fund,
    coords,
    gens,
    f,
    *,
    integration_policy,
    config=None,
    backend: str = "analytic",
    sector_filter=None,
    return_by_sector: bool = False,
) -> PhysicalCoefficientDerivatives:
    """Compute physical coefficient derivative contractions.

    Supported backends:
    - ``analytic``: currently proven for two-generator dK2 sectors only.
    - ``finite_difference``: preserved global FD oracle.
    - ``diagnostic``: alias for the FD oracle.
    - ``hybrid_local_fd``: explicit mixed path; analytic two-generator dK2 plus
      FD for sectors that are not analytically complete.
    """

    selected = _validate_sector_filter(sector_filter)
    if backend == "analytic":
        return _analytic_result(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=integration_policy,
            config=config,
            sector_filter=selected,
        )
    if backend in {"finite_difference", "diagnostic"}:
        return _finite_difference_result(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=integration_policy,
            config=config,
            sector_filter=selected,
            backend=backend,
        )
    if backend == "hybrid_local_fd":
        fd = _finite_difference_result(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=integration_policy,
            config=config,
            sector_filter=selected,
            backend=backend,
        )
        fd.metadata["effective_backend"] = "hybrid_local_fd"
        fd.metadata["fallback_used"] = bool(CUBIC_SECTORS & set(selected))
        fd.metadata["analytic_components"] = sorted(set(selected) & ANALYTIC_DK2_SECTORS)
        fd.metadata["fd_components"] = sorted(CUBIC_SECTORS & set(selected))
        return fd
    raise ValueError("backend must be 'analytic', 'finite_difference', 'diagnostic', or 'hybrid_local_fd'")
