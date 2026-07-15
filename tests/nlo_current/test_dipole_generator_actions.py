from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.dipole_observable import (  # noqa: E402
    apply_generator_word_to_dipole,
    dipole,
    fd_generator_action_on_dipole,
    fd_generator_word_action_on_dipole,
    left_generator_action_on_dipole,
    right_generator_action_on_dipole,
)
from nlo_current.su3_adjoint import random_su3, su3_generators_fundamental  # noqa: E402


def _setup(seed: int = 61001, nsite: int = 3):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng) for _ in range(nsite)])
    return rng, gens, U_fund


def test_dipole_normalization_for_equal_endpoints() -> None:
    _, _, U_fund = _setup(seed=61002)

    for site in range(U_fund.shape[0]):
        np.testing.assert_allclose(dipole(U_fund, site, site), 1.0, atol=1e-12, rtol=1e-12)


def test_single_left_generator_action_matches_finite_difference() -> None:
    _, gens, U_fund = _setup(seed=61003)

    analytic = left_generator_action_on_dipole(U_fund, 0, 2, 0, 1, gens)
    finite_diff = fd_generator_action_on_dipole(U_fund, 0, 2, 0, 1, gens, side="L")

    np.testing.assert_allclose(analytic, finite_diff, atol=1e-9, rtol=1e-7)


def test_single_right_generator_action_matches_finite_difference() -> None:
    _, gens, U_fund = _setup(seed=61004)

    analytic = right_generator_action_on_dipole(U_fund, 1, 5, 0, 1, gens)
    finite_diff = fd_generator_action_on_dipole(U_fund, 1, 5, 0, 1, gens, side="R")

    np.testing.assert_allclose(analytic, finite_diff, atol=1e-9, rtol=1e-7)


def test_two_generator_words_match_finite_differences() -> None:
    _, gens, U_fund = _setup(seed=61005)
    cases = [
        (((0, 1), (1, 2)), ("L", "L")),
        (((0, 1), (1, 2)), ("L", "R")),
        (((0, 1), (1, 2)), ("R", "L")),
        (((0, 1), (1, 2)), ("R", "R")),
    ]

    for word, sides in cases:
        analytic = apply_generator_word_to_dipole(U_fund, word, 0, 1, gens, sides)
        finite_diff = fd_generator_word_action_on_dipole(
            U_fund,
            word,
            0,
            1,
            gens,
            sides,
            eps=2e-5,
        )
        np.testing.assert_allclose(analytic, finite_diff, atol=5e-7, rtol=5e-5)


def test_three_generator_words_match_representative_finite_differences() -> None:
    _, gens, U_fund = _setup(seed=61006)
    cases = [
        (((0, 1), (1, 2), (0, 3)), ("L", "L", "R")),
        (((0, 1), (1, 2), (0, 3)), ("L", "R", "R")),
        (((0, 1), (1, 2), (0, 3)), ("L", "L", "L")),
        (((0, 1), (1, 2), (0, 3)), ("R", "R", "R")),
    ]

    for word, sides in cases:
        analytic = apply_generator_word_to_dipole(U_fund, word, 0, 1, gens, sides)
        finite_diff = fd_generator_word_action_on_dipole(
            U_fund,
            word,
            0,
            1,
            gens,
            sides,
            eps=8e-4,
        )
        np.testing.assert_allclose(analytic, finite_diff, atol=2e-5, rtol=2e-3)
