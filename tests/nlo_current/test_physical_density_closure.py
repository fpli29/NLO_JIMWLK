from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.nlo_current_skeleton import NLOCurrentTerms  # noqa: E402
from nlo_current.physical_current_divergence import evaluate_current_divergence  # noqa: E402
from nlo_current.physical_density_closure import (  # noqa: E402
    compare_physical_density_closure,
    raw_cubic_terms_transform,
)
from nlo_current.physical_density_operator import evaluate_direct_density_operator  # noqa: E402
from nlo_current.physical_kernels import KJSJIntegrationPolicy  # noqa: E402
from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.test_densities import (  # noqa: E402
    compute_test_density_hessian_score,
    compute_test_density_score,
    density_log_weight,
    evaluate_test_density,
)


def _setup(seed=20260725, nsite=1):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U = np.stack([random_su3(rng) for _ in range(nsite)])
    return rng, gens, f, U


def _policy(nsite):
    return KJSJIntegrationPolicy(
        quadrature_weights=np.ones(nsite) / nsite,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
        description="density closure unit-test policy",
    )


def _config(**overrides):
    params = {
        "nf": 2,
        "alpha_s": 0.3,
        "singularity_policy": "eps",
        "eps": 1.0e-6,
        "fd_eps_first": 2.0e-5,
        "fd_eps_second": 5.0e-4,
    }
    params.update(overrides)
    return PhysicalNLOCurrentConfig(**params)


def _sparse_terms(dim=8, *, K1=True, K2=True, K3=True):
    k1 = np.zeros(dim)
    k2 = np.zeros((dim, dim))
    k3 = np.zeros((dim, dim, dim))
    if K1:
        k1[0] = 0.2
        k1[3] = -0.1
    if K2:
        k2[0, 1] = 0.05
        k2[2, 3] = -0.04
    if K3:
        k3[0, 1, 2] = 0.03
        k3[2, 3, 4] = -0.02
    return NLOCurrentTerms(k1, k2, k3, {"sectors": {"synthetic": {"nonproduction": True}}})


def _zero_derivatives(dim):
    return {
        "dK2": np.zeros(dim),
        "dK3_first": {"LC_K3_ABC": np.zeros((dim, dim)), "LB_K3_ABC": np.zeros((dim, dim))},
        "d2K3": np.zeros(dim),
    }


def _density_builder(gens, name="single_link_trace"):
    return lambda V: evaluate_test_density(V, gens, name, score_eps=2.0e-5, hessian_eps=5.0e-4)


def _synthetic_closure(U, gens, terms, *, active=(0, 2), omit_hessian=False, transform=None):
    density = _density_builder(gens)
    base_density = density(U)
    derivs = _zero_derivatives(terms.dim)
    direct = evaluate_direct_density_operator(
        U,
        density,
        terms,
        gens=gens,
        fd_eps=1.0e-3,
        active_outer_indices=active,
    )
    current = evaluate_current_divergence(
        U,
        density,
        base_density.score,
        base_density.hessian_score,
        terms,
        derivs,
        gens=gens,
        fd_eps=1.0e-3,
        density_builder=density,
        active_outer_indices=active,
        omit_hessian_score=omit_hessian,
        terms_transform=transform,
    )
    return direct, current


def test_positive_test_densities_return_finite_positive_weights() -> None:
    _, gens, _, U3 = _setup(seed=20260726, nsite=3)
    for name in ("single_link_trace", "dipole_trace", "multilink_nonlinear", "constant"):
        result = evaluate_test_density(U3, gens, name)
        assert np.isfinite(result.log_weight)
        assert np.isfinite(result.weight)
        assert result.weight > 0.0
        assert np.all(np.isfinite(result.score))
        assert np.all(np.isfinite(result.hessian_score))


def test_score_and_ordered_hessian_are_finite_difference_consistent() -> None:
    _, gens, _, U = _setup(seed=20260727, nsite=1)
    score = compute_test_density_score(U, gens, "single_link_trace", eps=2.0e-5)
    hessian = compute_test_density_hessian_score(U, gens, "single_link_trace", eps=5.0e-4)
    result = evaluate_test_density(U, gens, "single_link_trace", score_eps=2.0e-5, hessian_eps=5.0e-4)
    np.testing.assert_allclose(result.score, score, atol=1.0e-12, rtol=1.0e-12)
    np.testing.assert_allclose(result.hessian_score, hessian, atol=1.0e-12, rtol=1.0e-12)
    assert np.linalg.norm(result.hessian_score - result.hessian_score.T) > 1.0e-4


def test_first_order_only_closure() -> None:
    _, gens, _, U = _setup(seed=20260728)
    direct, current = _synthetic_closure(U, gens, _sparse_terms(K2=False, K3=False), active=(0, 3))
    np.testing.assert_allclose(current.value, direct.value, atol=2.0e-9, rtol=2.0e-6)


def test_second_order_closure_with_K3_zero() -> None:
    _, gens, _, U = _setup(seed=20260729)
    direct, current = _synthetic_closure(U, gens, _sparse_terms(K3=False), active=(0, 2))
    np.testing.assert_allclose(current.value, direct.value, atol=2.0e-9, rtol=2.0e-6)


def test_cubic_synthetic_closure_uses_nonzero_hessian_score() -> None:
    _, gens, _, U = _setup(seed=20260730)
    terms = _sparse_terms()
    direct, current = _synthetic_closure(U, gens, terms, active=(0, 2))
    np.testing.assert_allclose(current.value, direct.value, atol=2.0e-9, rtol=2.0e-6)

    direct_no_h, current_no_h = _synthetic_closure(U, gens, terms, active=(0, 2), omit_hessian=True)
    assert abs(direct_no_h.value - current_no_h.value) > abs(direct.value - current.value) + 1.0e-7


def test_physical_second_order_projected_closure_on_smallest_setup() -> None:
    _, gens, f, U = _setup(seed=20260731, nsite=2)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)
    result = compare_physical_density_closure(
        U,
        "dipole_trace",
        coords=coords,
        gens=gens,
        f=f,
        physical_policy=_policy(2),
        config=_config(),
        sector_mask=("KJSSJ", "Kqbarq"),
        fd_eps=1.0e-3,
        active_outer_indices=(0,),
        include_by_sector=True,
    )
    assert result.abs_residual < 5.0e-8
    assert result.rel_residual < 5.0e-5
    assert set(result.sector_residuals) == {"KJSSJ", "Kqbarq"}


def test_omitting_coefficient_derivatives_fails_generically() -> None:
    _, gens, f, U = _setup(seed=20260732, nsite=2)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)
    baseline = compare_physical_density_closure(
        U,
        "dipole_trace",
        coords=coords,
        gens=gens,
        f=f,
        physical_policy=_policy(2),
        config=_config(),
        sector_mask=("KJSSJ", "Kqbarq"),
        fd_eps=1.0e-3,
        active_outer_indices=(0,),
    )
    omitted = compare_physical_density_closure(
        U,
        "dipole_trace",
        coords=coords,
        gens=gens,
        f=f,
        physical_policy=_policy(2),
        config=_config(),
        sector_mask=("KJSSJ", "Kqbarq"),
        fd_eps=1.0e-3,
        active_outer_indices=(0,),
        omit_coefficient_derivatives=True,
    )
    assert omitted.abs_residual > baseline.abs_residual + 1.0e-7


def test_raw_cubic_normalization_fails_and_preserves_complex_diagnostics() -> None:
    _, gens, _, U = _setup(seed=20260733)
    terms = _sparse_terms(K1=False, K2=False, K3=True)
    direct, current = _synthetic_closure(U, gens, terms, active=(0, 2))
    direct_raw, current_raw = _synthetic_closure(
        U,
        gens,
        terms,
        active=(0, 2),
        transform=raw_cubic_terms_transform,
    )
    assert abs(direct.value - current.value) < 2.0e-9
    assert abs(direct_raw.value - current_raw.value) > 1.0e-7
    assert abs(current_raw.value.imag) > 1.0e-8


def test_commutator_omission_changes_physical_closure_when_exercised() -> None:
    _, gens, _, U = _setup(seed=20260734)
    complete_terms = _sparse_terms(K1=False, K2=True, K3=True)
    no_comm_terms = _sparse_terms(K1=False, K2=False, K3=True)
    density = _density_builder(gens)
    base_density = density(U)
    direct = evaluate_direct_density_operator(
        U,
        density,
        complete_terms,
        gens=gens,
        fd_eps=1.0e-3,
        active_outer_indices=(0,),
    )
    current_complete = evaluate_current_divergence(
        U,
        density,
        base_density.score,
        base_density.hessian_score,
        complete_terms,
        _zero_derivatives(complete_terms.dim),
        gens=gens,
        fd_eps=1.0e-3,
        density_builder=density,
        active_outer_indices=(0,),
    )
    current_no_comm = evaluate_current_divergence(
        U,
        density,
        base_density.score,
        base_density.hessian_score,
        no_comm_terms,
        _zero_derivatives(no_comm_terms.dim),
        gens=gens,
        fd_eps=1.0e-3,
        density_builder=density,
        active_outer_indices=(0,),
    )
    assert abs(direct.value - current_complete.value) < 2.0e-9
    assert abs(direct.value - current_no_comm.value) > 1.0e-7


def test_constant_density_limit_and_finite_difference_window() -> None:
    _, gens, _, U = _setup(seed=20260735)
    density = _density_builder(gens, "constant")
    terms = _sparse_terms(K3=False)
    direct = evaluate_direct_density_operator(
        U,
        density,
        terms,
        gens=gens,
        fd_eps=1.0e-3,
        active_outer_indices=(0, 2),
    )
    current = evaluate_current_divergence(
        U,
        density,
        density(U).score,
        density(U).hessian_score,
        terms,
        _zero_derivatives(terms.dim),
        gens=gens,
        fd_eps=1.0e-3,
        density_builder=density,
        active_outer_indices=(0, 2),
    )
    np.testing.assert_allclose(direct.value, 0.0, atol=1.0e-12)
    np.testing.assert_allclose(current.value, 0.0, atol=1.0e-12)

    residuals = []
    for eps in (2.0e-3, 1.0e-3, 5.0e-4):
        density_nonconstant = _density_builder(gens)
        terms_nonconstant = _sparse_terms()
        d = evaluate_direct_density_operator(
            U,
            density_nonconstant,
            terms_nonconstant,
            gens=gens,
            fd_eps=eps,
            active_outer_indices=(0, 2),
        )
        c = evaluate_current_divergence(
            U,
            density_nonconstant,
            density_nonconstant(U).score,
            density_nonconstant(U).hessian_score,
            terms_nonconstant,
            _zero_derivatives(terms_nonconstant.dim),
            gens=gens,
            fd_eps=eps,
            density_builder=density_nonconstant,
            active_outer_indices=(0, 2),
        )
        residuals.append(abs(d.value - c.value))
    assert min(residuals) < 1.0e-8


def test_no_complex_warning_or_silent_imaginary_loss() -> None:
    _, gens, _, U = _setup(seed=20260736)
    terms = _sparse_terms(K1=False, K2=False, K3=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        direct, current = _synthetic_closure(
            U,
            gens,
            terms,
            active=(0, 2),
            transform=raw_cubic_terms_transform,
        )
    assert not any("ComplexWarning" in item.category.__name__ for item in caught)
    assert abs(current.value.imag) > 1.0e-8
    assert isinstance(direct.value, complex)
    assert isinstance(current.value, complex)


def test_density_log_weight_constant_is_exact() -> None:
    _, _, _, U = _setup(seed=20260737)
    assert density_log_weight(U, "constant") == 0.0
