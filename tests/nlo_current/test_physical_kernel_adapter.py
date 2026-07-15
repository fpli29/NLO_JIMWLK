from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.nlo_current_skeleton import assemble_nlo_current_terms  # noqa: E402
from nlo_current.physical_kernel_adapter import (  # noqa: E402
    finite_kernel_stats,
    physical_kernel_metadata,
    physical_kernels_for_skeleton,
)
from nlo_current.physical_cubic_conventions import klm_normalized_cubic_kernel  # noqa: E402
from nlo_current.physical_kernels import KJSJIntegrationPolicy, build_all_unbarred_physical_kernels  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


def _coords4() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0],
            [1.0, 0.2],
            [0.3, 1.1],
            [1.4, 1.3],
        ],
        dtype=float,
    )


def _policy(ncoords: int) -> KJSJIntegrationPolicy:
    return KJSJIntegrationPolicy(
        quadrature_weights=np.ones(ncoords) / ncoords,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
        description="adapter unit-test policy",
    )


def test_adapter_returns_expected_implemented_keys_and_metadata() -> None:
    kernels = physical_kernels_for_skeleton(
        _coords4(),
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_policy(4),
    )

    assert {"KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ", "metadata"}.issubset(kernels)
    assert kernels["metadata"]["implemented_kernels"] == ["KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ"]
    assert kernels["metadata"]["pending_kernels"] == []
    assert "positivity checks are future work" in kernels["metadata"]["positivity_note"]
    assert kernels["metadata"]["kjsj_integration_policy"]["mu"] == 1.3


def test_adapter_requires_explicit_kjsj_policy() -> None:
    try:
        physical_kernels_for_skeleton(_coords4(), Nc=3, nf=2, alpha_s=0.3)
    except ValueError as exc:
        assert "requires KJSJIntegrationPolicy" in str(exc)
    else:
        raise AssertionError("physical_kernels_for_skeleton should require KJSJIntegrationPolicy")


def test_physical_kernel_metadata_records_pending_kernels() -> None:
    metadata = physical_kernel_metadata(
        _coords4(),
        {
            "Nc": 3,
            "nf": 2,
            "alpha_s": 0.3,
            "singularity_policy": "eps",
            "eps": 1e-6,
            "integration_policy": _policy(4),
        },
    )

    assert metadata["coordinate_count"] == 4
    assert metadata["pending_kernels"] == []
    assert metadata["singularity_policy"] == "eps"
    assert metadata["kjsj_integration_policy"]["exclude_coincident_labels"] == ["x", "y", "z"]


def test_adapter_output_can_pass_to_skeleton_metadata_only_mode() -> None:
    rng = np.random.default_rng(20260711)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(4)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    kernels = physical_kernels_for_skeleton(
        _coords4(),
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_policy(4),
    )

    terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f, metadata_only=True)

    for key in ("KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ"):
        assert terms.metadata["sectors"][key]["metadata_only"] is True
    assert terms.metadata["sectors"]["KJSJ"]["metadata_only"] is True
    assert terms.metadata["nonproduction_only"] is True


def test_finite_kernel_stats_reports_shapes_and_counts() -> None:
    kernels = physical_kernels_for_skeleton(
        _coords4(),
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_policy(4),
    )
    stats = finite_kernel_stats(kernels)

    assert stats["KJSJ"]["shape"] == (4, 4, 4)
    assert stats["KJSSJ"]["shape"] == (4, 4, 4, 4)
    assert stats["Kqbarq"]["shape"] == (4, 4, 4, 4)
    assert stats["KJJSJ"]["shape"] == (4, 4, 4, 4)
    assert stats["KJJSSJ"]["shape"] == (4, 4, 4, 4, 4)
    assert all(item["finite_count"] > 0 for item in stats.values())


def test_physical_adapter_klm_normalizes_raw_complex_cubic_kernels() -> None:
    coords = _coords4()[:3]
    policy = _policy(coords.shape[0])
    raw = build_all_unbarred_physical_kernels(
        coords,
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=policy,
    )
    adapted = physical_kernels_for_skeleton(
        coords,
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=policy,
    )

    assert np.max(np.abs(np.imag(raw["KJJSJ"]))) > 0.0
    assert np.max(np.abs(np.imag(raw["KJJSSJ"]))) > 0.0
    np.testing.assert_allclose(adapted["KJJSJ"], klm_normalized_cubic_kernel(raw["KJJSJ"]))
    np.testing.assert_allclose(adapted["KJJSSJ"], klm_normalized_cubic_kernel(raw["KJJSSJ"]))
    assert not np.iscomplexobj(adapted["KJJSJ"])
    assert not np.iscomplexobj(adapted["KJJSSJ"])


def test_physical_full_assembly_emits_no_complex_cast_warning() -> None:
    coords = _coords4()[:3]
    rng = np.random.default_rng(20260712)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(coords.shape[0])])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    kernels = physical_kernels_for_skeleton(
        coords,
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_policy(coords.shape[0]),
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f, metadata_only=False)

    assert not caught
    assert not np.iscomplexobj(terms.K1)
    assert not np.iscomplexobj(terms.K2)
    assert not np.iscomplexobj(terms.K3)
    assert np.linalg.norm(terms.K3) > 0.0


def test_raw_complex_cubic_skeleton_path_preserves_imaginary_diagnostics() -> None:
    coords = _coords4()[:3]
    rng = np.random.default_rng(20260713)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(coords.shape[0])])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    raw = build_all_unbarred_physical_kernels(
        coords,
        Nc=3,
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_policy(coords.shape[0]),
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        terms = assemble_nlo_current_terms(U_fund, S_adj, raw, gens, f, metadata_only=False)

    assert not caught
    assert np.iscomplexobj(terms.K3)
    assert np.max(np.abs(np.imag(terms.K3))) > 0.0
