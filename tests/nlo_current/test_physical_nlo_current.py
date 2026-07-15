from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.nlo_velocity_evaluator import evaluate_velocity_from_terms  # noqa: E402
from nlo_current.physical_kernels import KJSJIntegrationPolicy  # noqa: E402
from nlo_current.physical_nlo_current import (  # noqa: E402
    PhysicalNLOCurrentConfig,
    assemble_physical_K1,
    assemble_physical_K2,
    assemble_physical_K3,
    assemble_physical_terms,
    compute_physical_coefficient_derivatives,
    evaluate_physical_nlo_velocity,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


def _coords2():
    return np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)


def _policy(nsite):
    return KJSJIntegrationPolicy(
        quadrature_weights=np.ones(nsite) / nsite,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
        description="physical current unit-test policy",
    )


def _config(**overrides):
    params = {
        "Nc": 3,
        "nf": 2,
        "alpha_s": 0.3,
        "singularity_policy": "eps",
        "eps": 1e-6,
        "fd_eps_first": 2e-5,
        "fd_eps_second": 5e-4,
    }
    params.update(overrides)
    return PhysicalNLOCurrentConfig(**params)


def _setup(seed=20260716, nsite=2):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(nsite)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    return rng, gens, f, U_fund, S_adj


def test_physical_K_arrays_have_expected_shapes_and_metadata() -> None:
    _, gens, f, U_fund, _ = _setup()
    coords = _coords2()
    terms = assemble_physical_terms(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=_policy(2),
        config=_config(),
    )

    assert assemble_physical_K1(U_fund, coords, gens, f, integration_policy=_policy(2), config=_config()).shape == (
        16,
    )
    assert assemble_physical_K2(U_fund, coords, gens, f, integration_policy=_policy(2), config=_config()).shape == (
        16,
        16,
    )
    assert assemble_physical_K3(U_fund, coords, gens, f, integration_policy=_policy(2), config=_config()).shape == (
        16,
        16,
        16,
    )
    assert terms.metadata["physical_nlo_current"]["kernel_origin"].startswith("unbarred physical")
    assert terms.metadata["physical_nlo_current"]["sector_labels"] == [
        "KJSJ",
        "KJSSJ",
        "Kqbarq",
        "KJJSJ",
        "KJJSSJ",
    ]


def test_LO_limit_K2_only_recovers_score_current_form() -> None:
    rng, gens, f, U_fund, _ = _setup(seed=20260717)
    coords = _coords2()
    terms = assemble_physical_terms(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=_policy(2),
        config=_config(),
        sector_filter=("KJSJ",),
    )
    score = rng.normal(size=terms.dim)
    hessian = rng.normal(size=(terms.dim, terms.dim))
    zeros = {
        "dK2": np.zeros(terms.dim),
        "dK3_first": {"LC_K3_ABC": np.zeros((terms.dim, terms.dim)), "LB_K3_ABC": np.zeros((terms.dim, terms.dim))},
        "d2K3": np.zeros(terms.dim),
    }

    result = evaluate_velocity_from_terms(
        terms,
        score,
        hessian,
        dK2=zeros["dK2"],
        dK3_first=zeros["dK3_first"],
        d2K3=zeros["d2K3"],
    )

    np.testing.assert_allclose(result["velocity"], -0.5 * terms.K2 @ score, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(terms.K1, 0.0, atol=1e-12)
    np.testing.assert_allclose(terms.K3, 0.0, atol=1e-12)


def test_zero_score_keeps_pure_coefficient_derivative_piece() -> None:
    _, gens, f, U_fund, _ = _setup(seed=20260718)
    coords = _coords2()
    score = np.zeros(16)
    hessian = np.zeros((16, 16))

    result = evaluate_physical_nlo_velocity(
        U_fund,
        coords,
        gens,
        f,
        score,
        hessian,
        integration_policy=_policy(2),
        config=_config(),
        sector_filter=("KJSJ",),
        derivative_backend="finite_difference",
    )
    derivatives = result["derivatives"]
    expected = result["terms"].K1 - 0.5 * derivatives["dK2"] + (1.0 / 6.0) * derivatives["d2K3"]

    np.testing.assert_allclose(result["velocity"], expected, atol=1e-10, rtol=1e-10)
    assert result["diagnostics"]["coefficient_derivatives_omitted"] is False


def test_zero_K3_sector_filter_recovers_second_order_current() -> None:
    rng, gens, f, U_fund, _ = _setup(seed=20260719)
    coords = _coords2()
    terms = assemble_physical_terms(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=_policy(2),
        config=_config(),
        sector_filter=("KJSJ", "KJSSJ", "Kqbarq"),
    )
    derivatives = compute_physical_coefficient_derivatives(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=_policy(2),
        config=_config(),
        sector_filter=("KJSJ", "KJSSJ", "Kqbarq"),
    )
    score = rng.normal(size=terms.dim)
    hessian = rng.normal(size=(terms.dim, terms.dim))
    result = evaluate_velocity_from_terms(
        terms,
        score,
        hessian,
        dK2=derivatives["dK2"],
        dK3_first=derivatives["dK3_first"],
        d2K3=derivatives["d2K3"],
    )
    expected = terms.K1 - 0.5 * (derivatives["dK2"] + terms.K2 @ score)

    np.testing.assert_allclose(terms.K3, 0.0, atol=1e-12)
    np.testing.assert_allclose(derivatives["d2K3"], 0.0, atol=1e-12)
    np.testing.assert_allclose(result["velocity"], expected, atol=1e-10, rtol=1e-10)


def test_cubic_normalization_leaves_no_imaginary_residual() -> None:
    _, gens, f, U_fund, _ = _setup(seed=20260720)
    terms = assemble_physical_terms(
        U_fund,
        _coords2(),
        gens,
        f,
        integration_policy=_policy(2),
        config=_config(),
    )

    assert not np.iscomplexobj(terms.K3)
    assert terms.metadata["physical_nlo_current"]["max_imag"]["K3"] == 0.0
    assert np.linalg.norm(terms.K3) > 0.0


def test_sector_by_sector_physical_contributions_sum_to_full() -> None:
    _, gens, f, U_fund, _ = _setup(seed=20260721)
    coords = _coords2()
    config = _config()
    full = assemble_physical_terms(U_fund, coords, gens, f, integration_policy=_policy(2), config=config)
    K1_sum = np.zeros_like(full.K1)
    K2_sum = np.zeros_like(full.K2)
    K3_sum = np.zeros_like(full.K3)
    for sector in ("KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ"):
        terms = assemble_physical_terms(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=_policy(2),
            config=config,
            sector_filter=(sector,),
        )
        K1_sum += terms.K1
        K2_sum += terms.K2
        K3_sum += terms.K3

    np.testing.assert_allclose(K1_sum, full.K1, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(K2_sum, full.K2, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(K3_sum, full.K3, atol=1e-12, rtol=1e-12)


def test_finite_difference_backend_alias_and_analytic_status() -> None:
    _, gens, f, U_fund, _ = _setup(seed=20260722)
    coords = _coords2()
    config = _config()
    finite = compute_physical_coefficient_derivatives(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=_policy(2),
        config=config,
        sector_filter=("KJJSJ",),
        backend="finite_difference",
    )
    diagnostic = compute_physical_coefficient_derivatives(
        U_fund,
        coords,
        gens,
        f,
        integration_policy=_policy(2),
        config=config,
        sector_filter=("KJJSJ",),
        backend="diagnostic",
    )

    assert finite["dK2"].shape == (16,)
    assert finite["dK3_first"]["LC_K3_ABC"].shape == (16, 16)
    assert finite["d2K3"].shape == (16,)
    np.testing.assert_allclose(finite["dK2"], diagnostic["dK2"], atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(finite["d2K3"], diagnostic["d2K3"], atol=1e-12, rtol=1e-12)
    assert finite["metadata"]["effective_backend"] == "finite_difference"

    with pytest.raises(NotImplementedError):
        compute_physical_coefficient_derivatives(
            U_fund,
            coords,
            gens,
            f,
            integration_policy=_policy(2),
            config=config,
            backend="analytic",
        )
