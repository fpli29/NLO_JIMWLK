from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.coefficient_derivatives import (  # noqa: E402
    compute_all_coefficient_derivatives_fd,
    compute_dK2_fd,
    product_rule_K2_rhs,
    product_rule_K3_rhs,
    validate_coefficient_shapes,
)
from nlo_current.finite_difference_scores import (  # noqa: E402
    fd_left_derivative_scalar,
    fd_left_second_derivative_scalar,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    su3_generators_fundamental,
)


def _S_builder(gens):
    return lambda U_fund: np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])


def _setup(seed: int = 51001, nsite: int = 1):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng) for _ in range(nsite)])
    dim = nsite * gens.shape[0]
    return rng, gens, U_fund, _S_builder(gens), dim


def _unflatten(index: int, n_color: int = 8) -> tuple[int, int]:
    return index // n_color, index % n_color


def _linear_trace_factor(Q):
    def phi(U_fund):
        return float(np.real(np.trace(Q @ U_fund[0])))

    return phi


def _score_and_hessian(logW, U_fund, gens, eps_score=2e-5, eps_hessian=5e-4):
    dim = len(U_fund) * gens.shape[0]
    U_list = [np.array(U, copy=True) for U in U_fund]
    score = np.zeros(dim)
    hessian = np.zeros((dim, dim))
    for b in range(dim):
        site_b, color_b = _unflatten(b, gens.shape[0])
        score[b] = fd_left_derivative_scalar(logW, U_list, site_b, color_b, gens, eps_score)
        for c in range(dim):
            site_c, color_c = _unflatten(c, gens.shape[0])
            hessian[b, c] = fd_left_second_derivative_scalar(
                logW,
                U_list,
                site_b,
                color_b,
                site_c,
                color_c,
                gens,
                eps_hessian,
            )
    return score, hessian


def test_constant_callbacks_give_zero_coefficient_derivatives() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=51002)
    K2_const = rng.normal(size=(dim, dim))
    K3_const = rng.normal(size=(dim, dim, dim))

    derivatives = compute_all_coefficient_derivatives_fd(
        lambda _U, _S: K2_const,
        lambda _U, _S: K3_const,
        U_fund,
        S_builder,
        gens,
    )

    np.testing.assert_allclose(derivatives["dK2"], 0.0, atol=1e-12)
    np.testing.assert_allclose(derivatives["dK3_first"]["LC_K3_ABC"], 0.0, atol=1e-12)
    np.testing.assert_allclose(derivatives["dK3_first"]["LB_K3_ABC"], 0.0, atol=1e-12)
    np.testing.assert_allclose(derivatives["d2K3"], 0.0, atol=1e-12)


def test_coefficient_derivative_shapes_for_two_sites() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=51003, nsite=2)
    K2_const = rng.normal(size=(dim, dim))
    K3_const = rng.normal(size=(dim, dim, dim))

    assert validate_coefficient_shapes(K2_const, K3_const) == 16
    derivatives = compute_all_coefficient_derivatives_fd(
        lambda _U, _S: K2_const,
        lambda _U, _S: K3_const,
        U_fund,
        S_builder,
        gens,
    )

    assert derivatives["dK2"].shape == (16,)
    assert derivatives["dK3_first"]["LC_K3_ABC"].shape == (16, 16)
    assert derivatives["dK3_first"]["LB_K3_ABC"].shape == (16, 16)
    assert derivatives["d2K3"].shape == (16,)


def test_linear_callback_gives_nonzero_stable_first_derivative() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=51004)
    Q = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    M = rng.normal(size=(dim, dim))
    phi = _linear_trace_factor(Q)

    def K2_callback(U, _S):
        return phi(U) * M

    dK2_eps1 = compute_dK2_fd(K2_callback, U_fund, S_builder, gens, eps=1e-5)
    dK2_eps2 = compute_dK2_fd(K2_callback, U_fund, S_builder, gens, eps=5e-6)

    assert np.linalg.norm(dK2_eps1) > 1e-8
    np.testing.assert_allclose(dK2_eps1, dK2_eps2, atol=1e-7, rtol=1e-5)


def test_product_rule_K2_contraction_matches_finite_difference_density_side() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=51005)
    Q_coeff = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    Q_density = 0.15 * (rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)))
    M = 0.2 * rng.normal(size=(dim, dim))
    phi = _linear_trace_factor(Q_coeff)

    def K2_callback(U, _S):
        return phi(U) * M

    def logW(Vs):
        return float(np.real(np.trace(Q_density @ Vs[0])))

    K2 = K2_callback(U_fund, S_builder(U_fund))
    dK2 = compute_dK2_fd(K2_callback, U_fund, S_builder, gens, eps=2e-5)
    score, _ = _score_and_hessian(logW, U_fund, gens)
    rhs = product_rule_K2_rhs(K2, dK2, score)

    U_list = [np.array(U, copy=True) for U in U_fund]
    W0 = np.exp(logW(U_list))
    lhs = np.zeros(dim)
    for a in range(dim):
        for b in range(dim):
            site_b, color_b = _unflatten(b, gens.shape[0])

            def density_term(Vs, a=a, b=b):
                V = np.stack(Vs)
                return float(K2_callback(V, S_builder(V))[a, b] * np.exp(logW(Vs)))

            lhs[a] += fd_left_derivative_scalar(
                density_term,
                U_list,
                site_b,
                color_b,
                gens,
                eps=2e-5,
            )
    lhs /= W0

    np.testing.assert_allclose(lhs, rhs, atol=5e-5, rtol=5e-5)


def test_product_rule_K3_contraction_matches_selected_finite_difference_indices() -> None:
    rng, gens, U_fund, S_builder, dim = _setup(seed=51006)
    Q_coeff = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    Q_density = 0.08 * (rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)))
    T = 0.04 * rng.normal(size=(dim, dim, dim))
    phi = _linear_trace_factor(Q_coeff)

    def K3_callback(U, _S):
        return phi(U) * T

    def logW(Vs):
        return float(np.real(np.trace(Q_density @ Vs[0])))

    eps_first = 2e-5
    eps_second = 7e-4
    derivatives = compute_all_coefficient_derivatives_fd(
        lambda _U, _S: np.zeros((dim, dim)),
        K3_callback,
        U_fund,
        S_builder,
        gens,
        eps_first=eps_first,
        eps_second=eps_second,
    )
    score, hessian = _score_and_hessian(
        logW,
        U_fund,
        gens,
        eps_score=eps_first,
        eps_hessian=eps_second,
    )
    K3 = K3_callback(U_fund, S_builder(U_fund))
    rhs = product_rule_K3_rhs(
        K3,
        derivatives["dK3_first"],
        derivatives["d2K3"],
        score,
        hessian,
    )

    U_list = [np.array(U, copy=True) for U in U_fund]
    W0 = np.exp(logW(U_list))
    for a in (0, 3):
        lhs_a = 0.0
        for b in range(dim):
            site_b, color_b = _unflatten(b, gens.shape[0])
            for c in range(dim):
                site_c, color_c = _unflatten(c, gens.shape[0])

                def density_term(Vs, a=a, b=b, c=c):
                    V = np.stack(Vs)
                    return float(K3_callback(V, S_builder(V))[a, b, c] * np.exp(logW(Vs)))

                lhs_a += fd_left_second_derivative_scalar(
                    density_term,
                    U_list,
                    site_b,
                    color_b,
                    site_c,
                    color_c,
                    gens,
                    eps_second,
                )
        lhs_a /= W0
        np.testing.assert_allclose(lhs_a, rhs[a], atol=2e-3, rtol=2e-3)
