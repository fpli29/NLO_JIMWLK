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
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.three_generator_terms import (  # noqa: E402
    kjjssj_A_LLR_from_kernel,
    kjjssj_B_LRR_from_kernel,
    kjjssj_V_virtual_from_kernel,
    synthetic_kjjssj_kernel,
)


def _setup(seed: int = 11111):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_list = [random_su3(rng) for _ in range(3)]
    kernel = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)
    params = {
        "q_mats": [rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)) for _ in range(3)],
        "lambda": 0.08,
        "eta": 0.13,
    }
    return gens, f, U_list, kernel, params


def _S_adj_from_list(Vs, gens):
    return np.stack([adjoint_from_fundamental(U, gens) for U in Vs])


def _relative(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1.0)


def test_kjjssj_llr_current_sign_distinct_sites() -> None:
    gens, f, U_list, kernel, params = _setup(seed=11111)
    logW = lambda Vs: toy_log_density(Vs, params)
    W = lambda Vs: float(np.exp(logW(Vs)))
    x, d = 0, 0
    y, e = 1, 1
    w, a = 2, 2
    eps_values = (1e-2, 5e-3, 2e-3)
    residuals = []

    def A_scalar(Vs):
        return float(kjjssj_A_LLR_from_kernel(_S_adj_from_list(Vs, gens), kernel, f)[x, y, w, d, e, a])

    for eps in eps_values:
        def inner_ll(Vs):
            return fd_left_second_derivative_scalar(
                lambda Xs: A_scalar(Xs) * W(Xs), Vs, y, e, x, d, gens, eps
            )

        direct = finite_diff_right_derivative(inner_ll, U_list, w, a, eps, gens)
        div = 0.0
        for h in range(8):
            def current_h(Vs, h=h):
                S_w = adjoint_from_fundamental(Vs[w], gens)
                return -S_w[h, a] * inner_ll(Vs)

            div -= finite_diff_left_derivative(current_h, U_list, w, h, eps, gens)
        residuals.append(_relative(direct, div))

    report = ROOT / "reports" / "nlo_current" / "kjjssj_cubic_current_fd_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# K_JJSSJ Cubic Current Finite-Difference Report\n\n"
        "## LLR sign residuals\n\n"
        "| eps | relative residual |\n"
        "|---:|---:|\n"
        + "\n".join(f"| {eps:.0e} | {res:.16e} |" for eps, res in zip(eps_values, residuals))
        + "\n",
        encoding="utf-8",
    )

    assert residuals[-1] < residuals[0]
    assert residuals[-1] < 1e-6


def test_kjjssj_lrr_current_sign_distinct_sites() -> None:
    gens, f, U_list, kernel, params = _setup(seed=11212)
    logW = lambda Vs: toy_log_density(Vs, params)
    W = lambda Vs: float(np.exp(logW(Vs)))
    w, a = 2, 2
    x, d = 0, 0
    y, e = 1, 1
    eps = 5e-3

    def B_scalar(Vs):
        return float(kjjssj_B_LRR_from_kernel(_S_adj_from_list(Vs, gens), kernel, f)[w, x, y, a, d, e])

    BW = lambda Vs: B_scalar(Vs) * W(Vs)
    Lw_BW = lambda Vs: finite_diff_left_derivative(BW, Vs, w, a, eps, gens)
    Rx = lambda func, Vs: finite_diff_right_derivative(func, Vs, x, d, eps, gens)
    Ry = lambda func, Vs: finite_diff_right_derivative(func, Vs, y, e, eps, gens)

    direct = Ry(lambda Vs: Rx(Lw_BW, Vs), U_list)
    right_right_BW = lambda Vs: Ry(lambda Zs: Rx(BW, Zs), Vs)
    div = -finite_diff_left_derivative(lambda Vs: -right_right_BW(Vs), U_list, w, a, eps, gens)

    assert _relative(direct, div) < 1e-8


def test_kjjssj_virtual_lll_rrr_signs_distinct_sites() -> None:
    gens, f, U_list, kernel, params = _setup(seed=11313)
    logW = lambda Vs: toy_log_density(Vs, params)
    W = lambda Vs: float(np.exp(logW(Vs)))
    V = kjjssj_V_virtual_from_kernel(kernel, f)
    x, c = 0, 0
    y, b = 1, 1
    w, a = 2, 2
    eps = 5e-3

    VW = lambda Vs: float(V[x, y, w, c, b, a]) * W(Vs)
    lll_direct = fd_left_derivative_scalar(
        lambda Vs: fd_left_second_derivative_scalar(VW, Vs, y, b, x, c, gens, eps),
        U_list,
        w,
        a,
        gens,
        eps,
    )
    lll_div = -fd_left_derivative_scalar(
        lambda Vs: -fd_left_second_derivative_scalar(VW, Vs, y, b, x, c, gens, eps),
        U_list,
        w,
        a,
        gens,
        eps,
    )

    Rx = lambda func, Vs: finite_diff_right_derivative(func, Vs, x, c, eps, gens)
    Ry = lambda func, Vs: finite_diff_right_derivative(func, Vs, y, b, eps, gens)
    Rw = lambda func, Vs: finite_diff_right_derivative(func, Vs, w, a, eps, gens)
    rrr_direct = -Rw(lambda Vs: Ry(lambda Zs: Rx(VW, Zs), Vs), U_list)
    rr_VW = lambda Vs: Ry(lambda Zs: Rx(VW, Zs), Vs)
    rrr_div = 0.0
    for h in range(8):
        def current_h(Vs, h=h):
            return adjoint_from_fundamental(Vs[w], gens)[h, a] * rr_VW(Vs)

        rrr_div -= finite_diff_left_derivative(current_h, U_list, w, h, eps, gens)

    assert _relative(lll_direct, lll_div) < 1e-12
    assert _relative(rrr_direct, rrr_div) < 1e-8


def test_kjjssj_hessian_score_contribution_nonzero() -> None:
    gens, f, U_list, kernel, params = _setup(seed=11414)
    logW = lambda Vs: toy_log_density(Vs, params)
    x, d = 0, 0
    y, e = 1, 1
    w, a = 2, 2
    eps = 1e-3

    def A_scalar(Vs):
        return float(kjjssj_A_LLR_from_kernel(_S_adj_from_list(Vs, gens), kernel, f)[x, y, w, d, e, a])

    hessian_term = A_scalar(U_list) * fd_hessian_score(logW, U_list, y, e, x, d, gens, eps)
    score_product = (
        A_scalar(U_list)
        * fd_score(logW, U_list, y, e, gens, eps)
        * fd_score(logW, U_list, x, d, gens, eps)
    )
    assert abs(hessian_term) > 1e-10
    assert np.isfinite(score_product)

