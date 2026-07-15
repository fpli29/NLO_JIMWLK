from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.analytic_lie_derivatives import (  # noqa: E402
    TraceFactor,
    left_derivative_adjoint,
    left_derivative_fundamental,
    left_derivative_fundamental_dagger,
    left_derivative_trace_word,
    second_left_derivative_adjoint,
    second_left_derivative_trace_word,
    trace_word_value,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    finite_diff_left_derivative,
    left_perturb,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


def _setup(seed=20260738, nsite=2):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U = np.stack([random_su3(rng) for _ in range(nsite)])
    S = np.stack([adjoint_from_fundamental(V, gens) for V in U])
    return rng, gens, f, U, S


def _fd_array(callback, U, site, color, gens, eps=1.0e-6):
    plus = np.array(U, copy=True)
    minus = np.array(U, copy=True)
    plus[site] = left_perturb(plus[site], color, eps, gens)
    minus[site] = left_perturb(minus[site], color, -eps, gens)
    return (callback(plus) - callback(minus)) / (2.0 * eps)


def _fd_second_array(callback, U, first_site, first_color, second_site, second_color, gens, eps=2.0e-4):
    return _fd_array(
        lambda V: _fd_array(callback, V, second_site, second_color, gens, eps=eps),
        U,
        first_site,
        first_color,
        gens,
        eps=eps,
    )


def test_left_derivative_fundamental_matches_finite_difference() -> None:
    _, gens, _, U, _ = _setup()
    analytic = left_derivative_fundamental(U, gens, site=0, color=2, target_site=0)
    finite = _fd_array(lambda V: V[0], U, 0, 2, gens)
    np.testing.assert_allclose(analytic, finite, atol=1.0e-10, rtol=1.0e-8)
    np.testing.assert_allclose(left_derivative_fundamental(U, gens, 1, 2, 0), 0.0, atol=1.0e-14)


def test_left_derivative_fundamental_dagger_matches_finite_difference() -> None:
    _, gens, _, U, _ = _setup(seed=20260739)
    analytic = left_derivative_fundamental_dagger(U, gens, site=1, color=5, target_site=1)
    finite = _fd_array(lambda V: V[1].conj().T, U, 1, 5, gens)
    np.testing.assert_allclose(analytic, finite, atol=1.0e-10, rtol=1.0e-8)
    np.testing.assert_allclose(left_derivative_fundamental_dagger(U, gens, 0, 5, 1), 0.0, atol=1.0e-14)


def test_left_derivative_adjoint_matches_finite_difference() -> None:
    _, gens, f, U, S = _setup(seed=20260740)
    analytic = left_derivative_adjoint(S, f, site=0, color=3, target_site=0)
    finite = _fd_array(lambda V: adjoint_from_fundamental(V[0], gens), U, 0, 3, gens)
    np.testing.assert_allclose(analytic, finite, atol=2.0e-10, rtol=2.0e-8)
    np.testing.assert_allclose(left_derivative_adjoint(S, f, 1, 3, 0), 0.0, atol=1.0e-14)


def test_ordered_second_left_derivative_adjoint_matches_finite_difference() -> None:
    _, gens, f, U, S = _setup(seed=20260741)
    analytic = second_left_derivative_adjoint(S, f, 0, 2, 0, 4, 0)
    finite = _fd_second_array(lambda V: adjoint_from_fundamental(V[0], gens), U, 0, 2, 0, 4, gens)
    np.testing.assert_allclose(analytic, finite, atol=2.0e-7, rtol=2.0e-5)


def test_trace_word_first_and_second_derivatives_match_finite_difference() -> None:
    rng, gens, _, U, _ = _setup(seed=20260742)
    probe = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    factors = [
        TraceFactor("Udag", site=0),
        TraceFactor("matrix", matrix=probe),
        TraceFactor("U", site=1),
        TraceFactor("matrix", matrix=gens[2]),
    ]
    analytic_first = left_derivative_trace_word(U, gens, factors, 0, 1)
    finite_first = finite_diff_left_derivative(
        lambda Vs: trace_word_value(np.stack(Vs), factors),
        [U[0], U[1]],
        0,
        1,
        eps=1.0e-6,
        gens=gens,
    )
    np.testing.assert_allclose(analytic_first, finite_first, atol=1.0e-10, rtol=1.0e-8)

    analytic_second = second_left_derivative_trace_word(U, gens, factors, 1, 3, 0, 1)
    finite_second = _fd_second_array(lambda V: trace_word_value(V, factors), U, 1, 3, 0, 1, gens)
    np.testing.assert_allclose(analytic_second, finite_second, atol=2.0e-7, rtol=2.0e-5)


def test_same_site_commutator_and_distinct_site_commutativity() -> None:
    rng, gens, f, U, _ = _setup(seed=20260743)
    probe = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))

    def scalar(Vs):
        return float(np.real(np.trace(probe @ Vs[0])))

    Us = [U[0], U[1]]
    a, b = 2, 5
    lhs = finite_diff_left_derivative(
        lambda Vs: finite_diff_left_derivative(scalar, Vs, 0, b, eps=2.0e-4, gens=gens),
        Us,
        0,
        a,
        eps=2.0e-4,
        gens=gens,
    ) - finite_diff_left_derivative(
        lambda Vs: finite_diff_left_derivative(scalar, Vs, 0, a, eps=2.0e-4, gens=gens),
        Us,
        0,
        b,
        eps=2.0e-4,
        gens=gens,
    )
    rhs = sum(
        f[a, b, c] * finite_diff_left_derivative(scalar, Us, 0, c, eps=2.0e-4, gens=gens)
        for c in range(8)
    )
    np.testing.assert_allclose(lhs, rhs, atol=2.0e-8, rtol=2.0e-6)

    distinct = finite_diff_left_derivative(
        lambda Vs: finite_diff_left_derivative(scalar, Vs, 1, b, eps=2.0e-4, gens=gens),
        Us,
        0,
        a,
        eps=2.0e-4,
        gens=gens,
    ) - finite_diff_left_derivative(
        lambda Vs: finite_diff_left_derivative(scalar, Vs, 0, a, eps=2.0e-4, gens=gens),
        Us,
        1,
        b,
        eps=2.0e-4,
        gens=gens,
    )
    np.testing.assert_allclose(distinct, 0.0, atol=1.0e-9)


def test_primitives_preserve_complex_values() -> None:
    rng, gens, _, U, _ = _setup(seed=20260744)
    probe = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    factors = [TraceFactor("matrix", matrix=probe), TraceFactor("U", site=0)]
    value = left_derivative_trace_word(U, gens, factors, 0, 1)
    assert isinstance(value, complex)
    assert abs(value.imag) > 1.0e-12

