from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.finite_difference_scores import fd_left_derivative_scalar  # noqa: E402
from nlo_current.lie_word_algebra import canonicalize_word  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


def _setup():
    rng = np.random.default_rng(12121)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_list = [random_su3(rng), random_su3(rng)]
    q0 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    q1 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))

    def scalar(Vs):
        return float(
            np.real(np.trace(q0 @ Vs[0]))
            + 0.3 * np.real(np.trace(Vs[0] @ Vs[1].conj().T))
            + 0.2 * np.real(np.trace(q1 @ Vs[1]))
        )

    return gens, f, U_list, scalar


def _apply_word(word, scalar_func, U_list, gens, eps):
    def rec(Vs, labels):
        if not labels:
            return scalar_func(Vs)
        site, color = labels[0]
        return fd_left_derivative_scalar(
            lambda Xs: rec(Xs, labels[1:]), Vs, site, color, gens, eps
        )

    return rec(U_list, tuple(word))


def _canonical_value(word, scalar_func, U_list, gens, f, eps):
    terms = canonicalize_word(tuple(word), f)
    return sum(coeff * _apply_word(term_word, scalar_func, U_list, gens, eps) for term_word, coeff in terms.items())


def _residuals_for_patterns(patterns):
    gens, f, U_list, scalar = _setup()
    eps_values = (1e-3, 1e-4, 1e-5)
    stable_eps_values = (1e-2, 5e-3, 2e-3, 1e-3)
    rows = []
    stable_max = {}
    for name, word in patterns.items():
        stable = []
        for eps in stable_eps_values:
            direct = _apply_word(word, scalar, U_list, gens, eps)
            canonical = _canonical_value(word, scalar, U_list, gens, f, eps)
            stable.append(abs(direct - canonical) / max(abs(direct), abs(canonical), 1.0))
        stable_max[name] = max(stable)
        for eps in eps_values:
            direct = _apply_word(word, scalar, U_list, gens, eps)
            canonical = _canonical_value(word, scalar, U_list, gens, f, eps)
            residual = abs(direct - canonical) / max(abs(direct), abs(canonical), 1.0)
            rows.append((name, eps, direct, canonical, residual))
    return rows, stable_max


def test_same_site_second_and_third_derivative_canonicalization_fd() -> None:
    patterns = {
        "same_site_second": ((0, 2), (0, 1)),
        "same_site_third": ((0, 2), (0, 1), (0, 0)),
    }
    rows, stable_max = _residuals_for_patterns(patterns)
    report = ROOT / "reports" / "nlo_current" / "lie_word_canonicalization_fd_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Lie Word Canonicalization Finite-Difference Report",
        "",
        "Residuals at eps=1e-5 are expected to be roundoff dominated for nested third derivatives.",
        "",
        "| pattern | eps | direct | canonical | relative residual |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, eps, direct, canonical, residual in rows:
        lines.append(f"| {name} | {eps:.0e} | {direct:.16e} | {canonical:.16e} | {residual:.16e} |")
    lines.append("")
    for name, residual in stable_max.items():
        lines.append(f"stable_max_{name}: {residual:.16e}")
    lines.append("")
    report.write_text("\n".join(lines), encoding="utf-8")

    assert stable_max["same_site_second"] < 2e-6
    assert stable_max["same_site_third"] < 5e-6


def test_mixed_coincident_patterns_canonicalization_fd() -> None:
    patterns = {
        "x_eq_y_ne_w": ((0, 2), (0, 1), (1, 0)),
        "x_eq_w_ne_y": ((0, 2), (1, 0), (0, 1)),
    }
    _, stable_max = _residuals_for_patterns(patterns)
    assert stable_max["x_eq_y_ne_w"] < 1e-6
    assert stable_max["x_eq_w_ne_y"] < 1e-6
