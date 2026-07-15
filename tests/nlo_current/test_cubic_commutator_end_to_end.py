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
    rng = np.random.default_rng(13131)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_list = [random_su3(rng), random_su3(rng)]
    q0 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    q1 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    r0 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))

    def scalar_F(Vs):
        log_w = 0.05 * np.real(np.trace(q0 @ Vs[0])) + 0.04 * np.real(np.trace(q1 @ Vs[1]))
        coeff = (
            0.17
            + 0.03 * np.real(np.trace(r0 @ Vs[0]))
            + 0.02 * np.real(np.trace(Vs[0] @ Vs[1].conj().T))
        )
        return float(coeff * np.exp(log_w))

    return gens, f, U_list, scalar_F


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
    return sum(
        coeff * _apply_word(term_word, scalar_func, U_list, gens, eps)
        for term_word, coeff in canonicalize_word(tuple(word), f).items()
    )


def test_cubic_commutator_end_to_end_patterns() -> None:
    gens, f, U_list, scalar_F = _setup()
    patterns = {
        "A_x_eq_y_ne_w": ((0, 2), (0, 1), (1, 0)),
        "B_x_eq_w_ne_y": ((0, 2), (1, 0), (0, 1)),
        "C_y_eq_w_ne_x": ((0, 0), (1, 2), (1, 1)),
        "D_x_eq_y_eq_w": ((0, 2), (0, 1), (0, 0)),
    }
    requested_eps = (1e-3, 1e-4, 1e-5)
    stable_eps = (1e-2, 5e-3, 2e-3, 1e-3)
    rows = []
    stable_max = {}

    for name, word in patterns.items():
        stable_residuals = []
        for eps in stable_eps:
            direct = _apply_word(word, scalar_F, U_list, gens, eps)
            canonical = _canonical_value(word, scalar_F, U_list, gens, f, eps)
            stable_residuals.append(abs(direct - canonical) / max(abs(direct), abs(canonical), 1.0))
        stable_max[name] = max(stable_residuals)

        for eps in requested_eps:
            direct = _apply_word(word, scalar_F, U_list, gens, eps)
            canonical = _canonical_value(word, scalar_F, U_list, gens, f, eps)
            residual = abs(direct - canonical) / max(abs(direct), abs(canonical), 1.0)
            rows.append((name, eps, direct, canonical, residual))

    report = ROOT / "reports" / "nlo_current" / "cubic_commutator_end_to_end_fd_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Cubic Commutator End-to-End Finite-Difference Report",
        "",
        "Direct ordered cubic derivatives are compared with symbolic canonicalized expressions.",
        "The eps=1e-5 entries are included for the workflow but are often roundoff dominated.",
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

    assert stable_max["A_x_eq_y_ne_w"] < 2e-6
    assert stable_max["B_x_eq_w_ne_y"] < 2e-6
    assert stable_max["C_y_eq_w_ne_x"] < 2e-6
    assert stable_max["D_x_eq_y_eq_w"] < 6e-6

