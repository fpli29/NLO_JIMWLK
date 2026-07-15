from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.coefficient_derivatives import (  # noqa: E402
    compute_all_coefficient_derivatives_fd,
    velocity_from_coeff_derivative_backend,
)
from nlo_current.nlo_current_skeleton import NLOCurrentTerms  # noqa: E402
from nlo_current.nlo_velocity_evaluator import evaluate_velocity_from_terms  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    su3_generators_fundamental,
)


def _S_builder(gens):
    return lambda U_fund: np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])


def _setup(seed: int = 52001):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng)])
    S_builder = _S_builder(gens)
    dim = gens.shape[0]
    return rng, gens, U_fund, S_builder, dim


def _nonconstant_callbacks(rng, dim):
    Q = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    M2 = 0.2 * rng.normal(size=(dim, dim))
    M3 = 0.05 * rng.normal(size=(dim, dim, dim))

    def phi(U_fund):
        return float(np.real(np.trace(Q @ U_fund[0])))

    def K2_callback(U_fund, _S_adj):
        return phi(U_fund) * M2

    def K3_callback(U_fund, _S_adj):
        return phi(U_fund) * M3

    return K2_callback, K3_callback


def test_computed_derivative_arrays_remove_velocity_omission_warnings() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=52002)
    K2_callback, K3_callback = _nonconstant_callbacks(rng, dim)
    terms = NLOCurrentTerms(
        K1=np.zeros(dim),
        K2=K2_callback(U_fund, S_builder(U_fund)),
        K3=K3_callback(U_fund, S_builder(U_fund)),
        metadata={},
    )
    derivatives = compute_all_coefficient_derivatives_fd(
        K2_callback,
        K3_callback,
        U_fund,
        S_builder,
        gens,
    )

    result = velocity_from_coeff_derivative_backend(
        terms,
        score=rng.normal(size=dim),
        hessian_score=rng.normal(size=(dim, dim)),
        derivatives=derivatives,
    )

    assert result["diagnostics"]["warnings"] == []
    assert result["diagnostics"]["coefficient_derivatives_omitted"] is False


def test_omitting_nonconstant_coefficient_derivatives_changes_velocity() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=52003)
    K2_callback, K3_callback = _nonconstant_callbacks(rng, dim)
    terms = NLOCurrentTerms(
        K1=np.zeros(dim),
        K2=K2_callback(U_fund, S_builder(U_fund)),
        K3=K3_callback(U_fund, S_builder(U_fund)),
        metadata={},
    )
    score = rng.normal(size=dim)
    hessian_score = rng.normal(size=(dim, dim))
    derivatives = compute_all_coefficient_derivatives_fd(
        K2_callback,
        K3_callback,
        U_fund,
        S_builder,
        gens,
    )

    without_derivatives = evaluate_velocity_from_terms(terms, score, hessian_score)
    with_derivatives = velocity_from_coeff_derivative_backend(
        terms,
        score,
        hessian_score,
        derivatives,
    )

    delta = np.linalg.norm(with_derivatives["velocity"] - without_derivatives["velocity"])
    assert delta > 1e-8


def test_constant_coefficients_match_derivative_free_velocity_values() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=52004)
    K2_const = rng.normal(size=(dim, dim))
    K3_const = rng.normal(size=(dim, dim, dim))
    terms = NLOCurrentTerms(
        K1=rng.normal(size=dim),
        K2=K2_const,
        K3=K3_const,
        metadata={},
    )
    derivatives = compute_all_coefficient_derivatives_fd(
        lambda _U, _S: K2_const,
        lambda _U, _S: K3_const,
        U_fund,
        S_builder,
        gens,
    )
    score = rng.normal(size=dim)
    hessian_score = rng.normal(size=(dim, dim))

    without_derivatives = evaluate_velocity_from_terms(terms, score, hessian_score)
    with_derivatives = velocity_from_coeff_derivative_backend(
        terms,
        score,
        hessian_score,
        derivatives,
    )

    np.testing.assert_allclose(
        with_derivatives["velocity"],
        without_derivatives["velocity"],
        atol=1e-12,
        rtol=1e-12,
    )
    assert without_derivatives["diagnostics"]["coefficient_derivatives_omitted"] is True
    assert with_derivatives["diagnostics"]["coefficient_derivatives_omitted"] is False
