from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.coordinate_kernels import (  # noqa: E402
    KernelSingularityError,
    cross2,
    dot,
    norm2,
    pairwise_dist2,
    safe_inv,
    validate_coords,
    vec,
)
from nlo_current.physical_kernels import (  # noqa: E402
    KJSJIntegrationPolicy,
    KJSJ_unbarred_local_value,
    KJSJ_unbarred_tilde_integral_value,
    KJJSSJ_unbarred_value,
    KJJSJ_unbarred_value,
    KJSSJ_unbarred_value,
    KJSJ_unbarred_value,
    Kqbarq_unbarred_value,
    build_KJJSSJ_unbarred,
    build_KJJSJ_unbarred,
    build_KJSSJ_unbarred,
    build_KJSJ_unbarred,
    build_Kqbarq_unbarred,
    tilde_K_JJSSJ_unbarred_value,
)


PARAMS = {"Nc": 3, "nf": 2, "alpha_s": 0.3}


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


def _coords5() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0],
            [1.0, 0.2],
            [0.3, 1.1],
            [1.4, 1.3],
            [2.1, 0.7],
        ],
        dtype=float,
    )


def _kjsj_policy(ncoords: int, *, mu: float = 1.3, weight_scale: float = 1.0):
    return KJSJIntegrationPolicy(
        quadrature_weights=weight_scale * np.ones(ncoords) / ncoords,
        mu=mu,
        exclude_coincident_labels=("x", "y", "z"),
        description="unit-test explicit finite z-prime sum",
    )


def test_coordinate_utility_behavior_and_singularity_policy() -> None:
    coords = _coords4()
    validate_coords(coords)
    np.testing.assert_allclose(vec(coords, 1, 0), np.array([1.0, 0.2]))
    assert norm2(np.array([3.0, 4.0])) == 25.0
    assert dot(np.array([1.0, 2.0]), np.array([3.0, 4.0])) == 11.0
    assert cross2(np.array([1.0, 0.0]), np.array([0.0, 2.0])) == 2.0
    assert pairwise_dist2(coords).shape == (4, 4)

    with pytest.raises(ValueError):
        validate_coords(np.ones((3, 3)))
    with pytest.raises(KernelSingularityError):
        safe_inv(0.0)
    assert np.isnan(safe_inv(0.0, singularity_policy="nan"))
    assert safe_inv(0.0, singularity_policy="eps", eps=0.5) == 2.0


def test_kjsj_requires_explicit_integration_policy_by_default() -> None:
    with pytest.raises(ValueError, match="K_JSJ requires an explicit"):
        KJSJ_unbarred_value(_coords4(), 0, 1, 2)


def test_dense_builder_shapes_for_implemented_kernels() -> None:
    coords = _coords4()
    params = dict(PARAMS, singularity_policy="nan")
    eps_params = dict(PARAMS, singularity_policy="eps", eps=1e-6, integration_policy=_kjsj_policy(4))
    assert build_KJSJ_unbarred(coords, **eps_params).shape == (4, 4, 4)
    assert build_KJSSJ_unbarred(coords, **params).shape == (4, 4, 4, 4)
    assert build_Kqbarq_unbarred(coords, **params).shape == (4, 4, 4, 4)
    assert build_KJJSJ_unbarred(coords, **params).shape == (4, 4, 4, 4)
    assert build_KJJSSJ_unbarred(coords, **params).shape == (4, 4, 4, 4, 4)


def test_pointwise_values_are_finite_away_from_singularities() -> None:
    coords = _coords5()
    values = [
        KJSJ_unbarred_value(
            coords,
            0,
            1,
            2,
            singularity_policy="eps",
            eps=1e-6,
            integration_policy=_kjsj_policy(5),
            **PARAMS,
        ),
        KJSSJ_unbarred_value(coords, 0, 1, 2, 3, **PARAMS),
        Kqbarq_unbarred_value(coords, 0, 1, 2, 3, **PARAMS),
        KJJSJ_unbarred_value(coords, 0, 1, 2, 3, **PARAMS),
        KJJSSJ_unbarred_value(coords, 0, 1, 2, 3, 4, **PARAMS),
        tilde_K_JJSSJ_unbarred_value(coords, 0, 1, 2, 3, **PARAMS),
    ]
    assert all(np.isfinite(v) for v in values)


def test_kernel_singularity_policy_raise_nan_and_eps() -> None:
    coords = _coords4()
    with pytest.raises(KernelSingularityError):
        Kqbarq_unbarred_value(coords, 0, 1, 2, 2, **PARAMS)
    assert np.isnan(
        Kqbarq_unbarred_value(coords, 0, 1, 2, 2, singularity_policy="nan", **PARAMS)
    )
    eps_value = Kqbarq_unbarred_value(
        coords,
        0,
        1,
        2,
        2,
        singularity_policy="eps",
        eps=1.0e-6,
        **PARAMS,
    )
    assert np.isfinite(eps_value)
    policy = _kjsj_policy(4)
    with pytest.raises(KernelSingularityError):
        KJSJ_unbarred_value(coords, 0, 1, 0, integration_policy=policy, **PARAMS)
    assert np.isfinite(
        KJSJ_unbarred_value(
            coords,
            0,
            1,
            0,
            singularity_policy="eps",
            eps=1e-6,
            integration_policy=policy,
            **PARAMS,
        )
    )


def test_confirmed_kernel_symmetries() -> None:
    coords = _coords5()
    x, y, z, zp = 0, 1, 2, 3
    kjsj_policy = _kjsj_policy(5)
    np.testing.assert_allclose(
        KJSJ_unbarred_value(
            coords,
            x,
            y,
            z,
            singularity_policy="eps",
            eps=1e-6,
            integration_policy=kjsj_policy,
            **PARAMS,
        ),
        KJSJ_unbarred_value(
            coords,
            y,
            x,
            z,
            singularity_policy="eps",
            eps=1e-6,
            integration_policy=kjsj_policy,
            **PARAMS,
        ),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        Kqbarq_unbarred_value(coords, x, y, z, zp, **PARAMS),
        Kqbarq_unbarred_value(coords, y, x, z, zp, **PARAMS),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        Kqbarq_unbarred_value(coords, x, y, z, zp, **PARAMS),
        Kqbarq_unbarred_value(coords, x, y, zp, z, **PARAMS),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        KJSSJ_unbarred_value(coords, x, y, z, zp, **PARAMS),
        KJSSJ_unbarred_value(coords, y, x, z, zp, **PARAMS),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        KJSSJ_unbarred_value(coords, x, y, z, zp, **PARAMS),
        KJSSJ_unbarred_value(coords, x, y, zp, z, **PARAMS),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        KJJSJ_unbarred_value(coords, 0, 1, 2, 3, **PARAMS),
        -KJJSJ_unbarred_value(coords, 0, 2, 1, 3, **PARAMS),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        KJJSSJ_unbarred_value(coords, 0, 1, 2, 3, 4, **PARAMS),
        -KJJSSJ_unbarred_value(coords, 0, 2, 1, 4, 3, **PARAMS),
        atol=1e-14,
    )


def test_tilde_k_helper_matches_definition() -> None:
    coords = _coords5()
    x, y, z, zp = 0, 1, 2, 3
    direct = tilde_K_JJSSJ_unbarred_value(coords, x, y, z, zp, **PARAMS)
    expected = 0.5j * (
        KJJSSJ_unbarred_value(coords, x, x, y, z, zp, **PARAMS)
        - KJJSSJ_unbarred_value(coords, y, x, y, z, zp, **PARAMS)
        - KJJSSJ_unbarred_value(coords, x, y, x, z, zp, **PARAMS)
        + KJJSSJ_unbarred_value(coords, y, y, x, z, zp, **PARAMS)
    )
    np.testing.assert_allclose(direct, expected, atol=1e-14)


def test_kjsj_diagonal_zero_condition() -> None:
    coords = _coords5()
    value = KJSJ_unbarred_value(
        coords,
        0,
        0,
        2,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_kjsj_policy(5),
        **PARAMS,
    )
    assert value == 0.0


def test_kjsj_tilde_integral_is_linear_in_quadrature_weights() -> None:
    coords = _coords5()
    base = KJSJ_unbarred_tilde_integral_value(
        coords,
        0,
        1,
        2,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_kjsj_policy(5, weight_scale=1.0),
        **PARAMS,
    )
    doubled = KJSJ_unbarred_tilde_integral_value(
        coords,
        0,
        1,
        2,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_kjsj_policy(5, weight_scale=2.0),
        **PARAMS,
    )
    np.testing.assert_allclose(doubled, 2.0 * base, atol=1e-14)


def test_kjsj_is_sensitive_to_mu_in_local_term() -> None:
    coords = _coords5()
    value_mu1 = KJSJ_unbarred_value(
        coords,
        0,
        1,
        2,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_kjsj_policy(5, mu=1.0),
        **PARAMS,
    )
    value_mu2 = KJSJ_unbarred_value(
        coords,
        0,
        1,
        2,
        singularity_policy="eps",
        eps=1e-6,
        integration_policy=_kjsj_policy(5, mu=2.0),
        **PARAMS,
    )
    assert abs(value_mu2 - value_mu1) > 1e-7


def test_kjsj_quadrature_refinement_is_convergent_for_diagnostic_ring() -> None:
    base = np.array([[0.0, 0.0], [1.0, 0.2], [0.3, 1.1]], dtype=float)

    def value(nring: int):
        angles = np.linspace(0.0, 2.0 * np.pi, nring, endpoint=False) + 0.13
        ring = np.column_stack([2.0 + 1.4 * np.cos(angles), 0.7 + 1.1 * np.sin(angles)])
        coords = np.vstack([base, ring])
        weights = np.zeros(coords.shape[0])
        weights[3:] = 2.0 * np.pi / nring
        policy = KJSJIntegrationPolicy(
            quadrature_weights=weights,
            mu=1.2,
            exclude_coincident_labels=("x", "y", "z"),
            description=f"ring quadrature n={nring}",
        )
        return KJSJ_unbarred_value(
            coords,
            0,
            1,
            2,
            singularity_policy="eps",
            eps=1e-7,
            integration_policy=policy,
            **PARAMS,
        )

    coarse_delta = abs(value(16) - value(8))
    refined_delta = abs(value(32) - value(16))
    assert refined_delta < coarse_delta


def test_kjsj_local_sign_matches_worknlo_prefactor_for_positive_bracket() -> None:
    coords = _coords5()
    value = KJSJ_unbarred_local_value(
        coords,
        0,
        1,
        2,
        Nc=3,
        nf=0,
        alpha_s=0.3,
        mu=3.0,
    )
    assert value < 0.0
