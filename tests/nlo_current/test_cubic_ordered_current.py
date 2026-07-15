from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.finite_difference_scores import (  # noqa: E402
    fd_hessian_score,
    fd_left_derivative_scalar,
    fd_left_second_derivative_scalar,
    fd_score,
    toy_log_density,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    finite_diff_left_derivative,
    finite_diff_right_derivative,
    random_su3,
    su3_generators_fundamental,
)


def _toy_setup(seed: int = 9101):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    U_list = [random_su3(rng) for _ in range(3)]
    params = {
        "q_mats": [rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)) for _ in range(3)],
        "lambda": 0.08,
        "eta": 0.13,
    }
    probes = [rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)) for _ in range(3)]
    return rng, gens, U_list, params, probes


def _make_scalar_coeff(probes, c0=0.03):
    def coeff(Vs):
        return float(
            c0
            + 0.020 * np.real(np.trace(probes[0] @ Vs[0]))
            - 0.017 * np.real(np.trace(probes[1] @ Vs[1]))
            + 0.011 * np.real(np.trace(probes[2] @ Vs[2]))
            + 0.015 * np.real(np.trace(Vs[0] @ Vs[1].conj().T))
        )

    return coeff


def test_toy_density_has_nonzero_offdiagonal_hessian_score() -> None:
    _, gens, U_list, params, _ = _toy_setup(seed=9001)
    logW = lambda Vs: toy_log_density(Vs, params)
    values = [
        fd_hessian_score(logW, U_list, 0, a, 1, b, gens, eps=1e-4)
        for a in range(3)
        for b in range(3)
    ]
    assert max(abs(v) for v in values) > 1e-4


def test_llr_cubic_expansion_identity() -> None:
    _, gens, U_list, params, probes = _toy_setup(seed=9101)
    logW = lambda Vs: toy_log_density(Vs, params)
    A = _make_scalar_coeff(probes)
    W = lambda Vs: float(np.exp(logW(Vs)))
    x, d = 0, 0
    y, e = 1, 1
    residuals = []

    for eps in (1e-2, 5e-3, 1e-3):
        lhs = (
            fd_left_second_derivative_scalar(lambda Vs: A(Vs) * W(Vs), U_list, y, e, x, d, gens, eps)
            / W(U_list)
        )
        rhs = (
            fd_left_second_derivative_scalar(A, U_list, y, e, x, d, gens, eps)
            + fd_left_derivative_scalar(A, U_list, x, d, gens, eps)
            * fd_score(logW, U_list, y, e, gens, eps)
            + fd_left_derivative_scalar(A, U_list, y, e, gens, eps)
            * fd_score(logW, U_list, x, d, gens, eps)
            + A(U_list)
            * (
                fd_hessian_score(logW, U_list, y, e, x, d, gens, eps)
                + fd_score(logW, U_list, y, e, gens, eps)
                * fd_score(logW, U_list, x, d, gens, eps)
            )
        )
        residuals.append(abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0))

    report = ROOT / "reports" / "nlo_current" / "kjjsj_llr_expansion_fd_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# K_JJSJ LLR Expansion Finite-Difference Report\n\n"
        "| eps | relative residual |\n"
        "|---:|---:|\n"
        + "\n".join(
            f"| {eps:.0e} | {residual:.16e} |"
            for eps, residual in zip((1e-2, 5e-3, 1e-3), residuals)
        )
        + "\n",
        encoding="utf-8",
    )

    assert residuals[-1] < residuals[0]
    assert residuals[-1] < 1e-8


def test_llr_current_sign_distinct_sites() -> None:
    _, gens, U_list, params, probes = _toy_setup(seed=9202)
    logW = lambda Vs: toy_log_density(Vs, params)
    A = _make_scalar_coeff(probes)
    W = lambda Vs: float(np.exp(logW(Vs)))
    x, d = 0, 0
    y, e = 1, 1
    w, a = 2, 2
    eps = 5e-3

    def inner_ll(Vs):
        return fd_left_second_derivative_scalar(lambda Xs: A(Xs) * W(Xs), Vs, y, e, x, d, gens, eps)

    direct = finite_diff_right_derivative(inner_ll, U_list, w, a, eps, gens)
    div = 0.0
    for h in range(8):
        def current_h(Vs, h=h):
            S_w = adjoint_from_fundamental(Vs[w], gens)
            return -S_w[h, a] * inner_ll(Vs)

        div -= finite_diff_left_derivative(current_h, U_list, w, h, eps, gens)

    residual = abs(direct - div) / max(abs(direct), abs(div), 1.0)
    assert residual < 1e-8


def test_lrr_current_sign_distinct_sites() -> None:
    _, gens, U_list, params, probes = _toy_setup(seed=9303)
    logW = lambda Vs: toy_log_density(Vs, params)
    B = _make_scalar_coeff(probes, c0=-0.02)
    W = lambda Vs: float(np.exp(logW(Vs)))
    w, a = 2, 2
    x, d = 0, 0
    y, e = 1, 1
    eps = 5e-3

    BW = lambda Vs: B(Vs) * W(Vs)
    Lw_BW = lambda Vs: finite_diff_left_derivative(BW, Vs, w, a, eps, gens)
    Rx = lambda func, Vs: finite_diff_right_derivative(func, Vs, x, d, eps, gens)
    Ry = lambda func, Vs: finite_diff_right_derivative(func, Vs, y, e, eps, gens)

    direct = Ry(lambda Vs: Rx(Lw_BW, Vs), U_list)
    right_right_BW = lambda Vs: Ry(lambda Zs: Rx(BW, Zs), Vs)
    current = lambda Vs: -right_right_BW(Vs)
    div = -finite_diff_left_derivative(current, U_list, w, a, eps, gens)

    residual = abs(direct - div) / max(abs(direct), abs(div), 1.0)
    assert residual < 1e-10


def test_virtual_lll_current_sign_smoke() -> None:
    _, gens, U_list, params, probes = _toy_setup(seed=9404)
    logW = lambda Vs: toy_log_density(Vs, params)
    V = _make_scalar_coeff(probes, c0=0.01)
    W = lambda Vs: float(np.exp(logW(Vs)))
    x, d = 0, 0
    y, e = 1, 1
    w, b = 2, 2
    eps = 5e-3

    VW = lambda Vs: V(Vs) * W(Vs)
    direct = fd_left_derivative_scalar(
        lambda Vs: fd_left_second_derivative_scalar(VW, Vs, y, e, x, d, gens, eps),
        U_list,
        w,
        b,
        gens,
        eps,
    )
    current = lambda Vs: -fd_left_second_derivative_scalar(VW, Vs, y, e, x, d, gens, eps)
    div = -fd_left_derivative_scalar(current, U_list, w, b, gens, eps)

    residual = abs(direct - div) / max(abs(direct), abs(div), 1.0)
    assert residual < 1e-12


def test_virtual_rrr_current_sign_smoke_distinct_sites() -> None:
    _, gens, U_list, params, probes = _toy_setup(seed=9454)
    logW = lambda Vs: toy_log_density(Vs, params)
    V = _make_scalar_coeff(probes, c0=0.01)
    W = lambda Vs: float(np.exp(logW(Vs)))
    x, d = 0, 0
    y, e = 1, 1
    w, b = 2, 2
    eps = 5e-3

    VW = lambda Vs: V(Vs) * W(Vs)
    Rx = lambda func, Vs: finite_diff_right_derivative(func, Vs, x, d, eps, gens)
    Ry = lambda func, Vs: finite_diff_right_derivative(func, Vs, y, e, eps, gens)
    Rw = lambda func, Vs: finite_diff_right_derivative(func, Vs, w, b, eps, gens)

    direct = -Rw(lambda Vs: Ry(lambda Zs: Rx(VW, Zs), Vs), U_list)

    def inner_rights(Vs):
        return Ry(lambda Zs: Rx(VW, Zs), Vs)

    div = 0.0
    for h in range(8):
        def current_h(Vs, h=h):
            S_w = adjoint_from_fundamental(Vs[w], gens)
            return S_w[h, b] * inner_rights(Vs)

        div -= finite_diff_left_derivative(current_h, U_list, w, h, eps, gens)

    residual = abs(direct - div) / max(abs(direct), abs(div), 1.0)
    assert residual < 1e-8
