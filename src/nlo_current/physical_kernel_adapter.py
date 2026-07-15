"""Adapter from physical coordinate kernels to the non-production skeleton API."""

from __future__ import annotations

from copy import deepcopy

import numpy as np

from .coordinate_kernels import validate_coords
from .physical_cubic_conventions import cubic_kernel_convention_diagnostics, klm_normalized_cubic_kernel
from .physical_kernels import (
    IMPLEMENTED_UNBARRED_KERNELS,
    PENDING_UNBARRED_KERNELS,
    KJSJIntegrationPolicy,
    build_all_unbarred_physical_kernels,
)


def physical_kernel_metadata(coords, params):
    """Return metadata for non-production physical kernel diagnostics."""

    arr = validate_coords(coords)
    p = deepcopy(dict(params))
    return {
        "coordinate_count": int(arr.shape[0]),
        "coordinate_shape": tuple(arr.shape),
        "singularity_policy": p.get("singularity_policy", "raise"),
        "eps": p.get("eps"),
        "kjsj_integration_policy": _policy_metadata(p.get("integration_policy")),
        "implemented_kernels": list(IMPLEMENTED_UNBARRED_KERNELS),
        "pending_kernels": list(PENDING_UNBARRED_KERNELS),
        "parameter_values": {
            "Nc": p.get("Nc", 3),
            "nf": p.get("nf", 0),
            "alpha_s": p.get("alpha_s", 1.0),
        },
        "kernel_type": "unbarred singlet",
        "nonproduction_only": True,
        "positivity_note": (
            "Physical-kernel positivity checks are future work; the Pawula toy "
            "diagnostic does not prove physical NLO JIMWLK positivity or non-positivity."
        ),
    }


def physical_kernels_for_skeleton(
    coords,
    *,
    Nc=3,
    nf=0,
    alpha_s=1.0,
    singularity_policy="raise",
    eps=None,
    integration_policy: KJSJIntegrationPolicy | None = None,
):
    """Build physical unbarred kernels for non-production skeleton diagnostics."""

    if integration_policy is None:
        raise ValueError("physical_kernels_for_skeleton requires KJSJIntegrationPolicy")
    params = {
        "Nc": Nc,
        "nf": nf,
        "alpha_s": alpha_s,
        "singularity_policy": singularity_policy,
        "eps": eps,
        "integration_policy": integration_policy,
    }
    raw_kernels = build_all_unbarred_physical_kernels(coords, **params)
    kernels = dict(raw_kernels)
    kernels["KJJSJ"] = klm_normalized_cubic_kernel(raw_kernels["KJJSJ"])
    kernels["KJJSSJ"] = klm_normalized_cubic_kernel(raw_kernels["KJJSSJ"])
    metadata = physical_kernel_metadata(coords, params)
    metadata["cubic_kernel_convention"] = {
        "adapter_output": "KLM-normalized cubic coefficients",
        "normalization": "KJJSJ,KJJSSJ adapter arrays = (-1j) * raw physical kernels",
        "raw_physical_kernel_note": "build_all_unbarred_physical_kernels returns the WORKNLO raw cubic formulas.",
        "KJJSJ": cubic_kernel_convention_diagnostics(raw_kernels["KJJSJ"], kernels["KJJSJ"]),
        "KJJSSJ": cubic_kernel_convention_diagnostics(raw_kernels["KJJSSJ"], kernels["KJJSSJ"]),
    }
    metadata["array_keys"] = sorted(key for key in kernels if key != "metadata")
    kernels["metadata"] = metadata
    return kernels


def _policy_metadata(policy):
    if policy is None:
        return None
    weights = np.asarray(policy.quadrature_weights, dtype=float)
    return {
        "mu": float(policy.mu),
        "quadrature_weight_count": int(weights.size),
        "quadrature_weight_sum": float(np.sum(weights)),
        "excluded_indices": list(policy.excluded_indices),
        "exclude_coincident_labels": list(policy.exclude_coincident_labels),
        "principal_value": policy.principal_value,
        "subtraction": policy.subtraction,
        "finite_volume_boundary": policy.finite_volume_boundary,
        "description": policy.description,
    }


def finite_kernel_stats(kernels: dict) -> dict[str, dict[str, float | int | tuple]]:
    """Return finite/nonfinite counts and finite-value norms for kernel arrays."""

    stats = {}
    for key, value in kernels.items():
        if key == "metadata":
            continue
        arr = np.asarray(value)
        finite = np.isfinite(arr)
        finite_values = arr[finite]
        stats[key] = {
            "shape": tuple(arr.shape),
            "finite_count": int(np.sum(finite)),
            "nonfinite_count": int(arr.size - np.sum(finite)),
            "finite_norm": float(np.linalg.norm(finite_values)) if finite_values.size else 0.0,
        }
    return stats
