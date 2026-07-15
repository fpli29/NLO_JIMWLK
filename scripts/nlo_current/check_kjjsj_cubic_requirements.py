#!/usr/bin/env python3
"""Diagnose K_JJSJ cubic current score/Hessian requirements on a tiny lattice."""

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
    kjjsj_A_LLR_from_kernel,
    kjjsj_B_LRR_from_kernel,
    kjjsj_V_virtual_from_kernel,
    synthetic_kjjsj_kernel,
)


def _toy_params(rng: np.random.Generator, nsite: int) -> dict:
    return {
        "q_mats": [rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)) for _ in range(nsite)],
        "lambda": 0.08,
        "eta": 0.13,
    }


def _toy_coeff(probes, c0):
    def coeff(Vs):
        return float(
            c0
            + 0.020 * np.real(np.trace(probes[0] @ Vs[0]))
            - 0.017 * np.real(np.trace(probes[1] @ Vs[1]))
            + 0.011 * np.real(np.trace(probes[2] @ Vs[2]))
            + 0.015 * np.real(np.trace(Vs[0] @ Vs[1].conj().T))
        )

    return coeff


def _term_components(coef_func, logW_func, U_list, gens, eps):
    x, d = 0, 0
    y, e = 1, 1
    second = fd_left_second_derivative_scalar(coef_func, U_list, y, e, x, d, gens, eps)
    score_linear = (
        fd_left_derivative_scalar(coef_func, U_list, x, d, gens, eps)
        * fd_score(logW_func, U_list, y, e, gens, eps)
        + fd_left_derivative_scalar(coef_func, U_list, y, e, gens, eps)
        * fd_score(logW_func, U_list, x, d, gens, eps)
    )
    hessian = coef_func(U_list) * fd_hessian_score(logW_func, U_list, y, e, x, d, gens, eps)
    score_product = (
        coef_func(U_list)
        * fd_score(logW_func, U_list, y, e, gens, eps)
        * fd_score(logW_func, U_list, x, d, gens, eps)
    )
    return {
        "second_derivative": abs(float(second)),
        "score_linear": abs(float(score_linear)),
        "hessian_score": abs(float(hessian)),
        "score_product": abs(float(score_product)),
    }


def _sign_residuals(U_list, gens, logW_func, probes):
    A = _toy_coeff(probes, 0.03)
    B = _toy_coeff(probes, -0.02)
    V = _toy_coeff(probes, 0.01)
    W = lambda Vs: float(np.exp(logW_func(Vs)))
    x, d = 0, 0
    y, e = 1, 1
    w, c = 2, 2
    eps = 5e-3

    def rel(a, b):
        return abs(a - b) / max(abs(a), abs(b), 1.0)

    def llr_inner(Vs):
        return fd_left_second_derivative_scalar(lambda Xs: A(Xs) * W(Xs), Vs, y, e, x, d, gens, eps)

    llr_direct = finite_diff_right_derivative(llr_inner, U_list, w, c, eps, gens)
    llr_div = 0.0
    for h in range(8):
        def current_h(Vs, h=h):
            return -adjoint_from_fundamental(Vs[w], gens)[h, c] * llr_inner(Vs)

        llr_div -= finite_diff_left_derivative(current_h, U_list, w, h, eps, gens)

    BW = lambda Vs: B(Vs) * W(Vs)
    Lw_BW = lambda Vs: finite_diff_left_derivative(BW, Vs, w, c, eps, gens)
    Rx = lambda func, Vs: finite_diff_right_derivative(func, Vs, x, d, eps, gens)
    Ry = lambda func, Vs: finite_diff_right_derivative(func, Vs, y, e, eps, gens)
    lrr_direct = Ry(lambda Vs: Rx(Lw_BW, Vs), U_list)
    rr_BW = lambda Vs: Ry(lambda Zs: Rx(BW, Zs), Vs)
    lrr_div = -finite_diff_left_derivative(lambda Vs: -rr_BW(Vs), U_list, w, c, eps, gens)

    VW = lambda Vs: V(Vs) * W(Vs)
    lll_direct = fd_left_derivative_scalar(
        lambda Vs: fd_left_second_derivative_scalar(VW, Vs, y, e, x, d, gens, eps),
        U_list,
        w,
        c,
        gens,
        eps,
    )
    lll_div = -fd_left_derivative_scalar(
        lambda Vs: -fd_left_second_derivative_scalar(VW, Vs, y, e, x, d, gens, eps),
        U_list,
        w,
        c,
        gens,
        eps,
    )

    Rw = lambda func, Vs: finite_diff_right_derivative(func, Vs, w, c, eps, gens)
    rrr_direct = -Rw(lambda Vs: Ry(lambda Zs: Rx(VW, Zs), Vs), U_list)
    rr_VW = lambda Vs: Ry(lambda Zs: Rx(VW, Zs), Vs)
    rrr_div = 0.0
    for h in range(8):
        def current_h(Vs, h=h):
            return adjoint_from_fundamental(Vs[w], gens)[h, c] * rr_VW(Vs)

        rrr_div -= finite_diff_left_derivative(current_h, U_list, w, h, eps, gens)

    return {
        "llr": rel(llr_direct, llr_div),
        "lrr": rel(lrr_direct, lrr_div),
        "virtual_lll": rel(lll_direct, lll_div),
        "virtual_rrr": rel(rrr_direct, rrr_div),
    }


def measure_case(seed: int, xy_symmetry: str):
    rng = np.random.default_rng(seed)
    nsite = 3
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_list = [random_su3(rng) for _ in range(nsite)]
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_list])
    kernel = synthetic_kjjsj_kernel(nsite, rng, xy_symmetry=xy_symmetry)
    A = kjjsj_A_LLR_from_kernel(S_adj, kernel, f)
    B = kjjsj_B_LRR_from_kernel(S_adj, kernel, f)
    V = kjjsj_V_virtual_from_kernel(kernel, f)
    params = _toy_params(rng, nsite)
    logW = lambda Vs: toy_log_density(Vs, params)
    eps = 1e-3

    x, d = 0, 0
    y, e = 1, 1
    w, a = 2, 2

    def A_scalar(Vs):
        S = np.stack([adjoint_from_fundamental(U, gens) for U in Vs])
        return float(kjjsj_A_LLR_from_kernel(S, kernel, f)[x, y, w, d, e, a])

    def B_scalar(Vs):
        S = np.stack([adjoint_from_fundamental(U, gens) for U in Vs])
        return float(kjjsj_B_LRR_from_kernel(S, kernel, f)[w, x, y, a, d, e])

    V_value = float(V[x, y, w, d, e, a])
    V_scalar = lambda Vs: V_value

    probes = [rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)) for _ in range(nsite)]
    sign_residuals = _sign_residuals(U_list, gens, logW, probes)

    return {
        "seed": seed,
        "xy_symmetry": xy_symmetry,
        "nsite": nsite,
        "block_norms": {
            "A_LLR": float(np.linalg.norm(A)),
            "B_LRR": float(np.linalg.norm(B)),
            "V_virtual": float(np.linalg.norm(V)),
        },
        "term_norms": {
            "LLR": _term_components(A_scalar, logW, U_list, gens, eps),
            "LRR": _term_components(B_scalar, logW, U_list, gens, eps),
            "virtual": _term_components(V_scalar, logW, U_list, gens, eps),
        },
        "hessian_probe": abs(fd_hessian_score(logW, U_list, y, e, x, d, gens, eps)),
        "sign_residuals": sign_residuals,
    }


def write_report(cases) -> Path:
    report = ROOT / "reports" / "nlo_current" / "kjjsj_cubic_requirements_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# K_JJSJ Cubic Requirements Report",
        "",
        "Dense small-lattice diagnostic for the first three-generator NLO block.",
        "The tested current signs are distinct-site ordered-block checks, not a",
        "production treatment of all coincident-site commutators.",
        "",
    ]
    for case in cases:
        lines.extend(
            [
                f"## Kernel {case['xy_symmetry']}",
                "",
                f"seed: {case['seed']}",
                f"nsite: {case['nsite']}",
                "",
                "Built coefficient blocks: A_LLR, B_LRR, V_virtual.",
                "",
                "| block | norm |",
                "|---|---:|",
            ]
        )
        for name, value in case["block_norms"].items():
            lines.append(f"| {name} | {value:.16e} |")
        lines.extend(["", "| component | second derivative | score-linear | Hessian-score | score-product |", "|---|---:|---:|---:|---:|"])
        for name, parts in case["term_norms"].items():
            lines.append(
                f"| {name} | {parts['second_derivative']:.16e} | "
                f"{parts['score_linear']:.16e} | {parts['hessian_score']:.16e} | "
                f"{parts['score_product']:.16e} |"
            )
        lines.extend(["", f"hessian_probe_abs: {case['hessian_probe']:.16e}", ""])
        lines.extend(["| sign check | relative residual |", "|---|---:|"])
        for name, value in case["sign_residuals"].items():
            lines.append(f"| {name} | {value:.16e} |")
        lines.append("")

    max_sign = max(max(case["sign_residuals"].values()) for case in cases)
    max_hessian = max(
        max(parts["hessian_score"] for parts in case["term_norms"].values()) for case in cases
    )
    lines.extend(
        [
            "## Conclusion",
            "",
            f"max_sign_residual: {max_sign:.16e}",
            f"max_hessian_score_component: {max_hessian:.16e}",
            "",
        ]
    )
    if max_sign < 1e-8 and max_hessian > 1e-8:
        lines.append(
            "K_JJSJ current can be represented as score + Hessian-score with the "
            "tested distinct-site ordered signs. Coincident-site production "
            "commutator handling remains outside this workflow."
        )
    else:
        lines.append("K_JJSJ current still has an unresolved ordering/sign issue.")
    lines.append("")
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main() -> None:
    cases = [
        measure_case(seed=20260724, xy_symmetry="symmetric"),
        measure_case(seed=20260725, xy_symmetry="antisymmetric"),
    ]
    report = write_report(cases)
    for case in cases:
        print(
            f"seed={case['seed']} xy_symmetry={case['xy_symmetry']} "
            f"hessian_probe={case['hessian_probe']:.16e} "
            f"max_sign_residual={max(case['sign_residuals'].values()):.16e}"
        )
    print(f"wrote {report}")


if __name__ == "__main__":
    main()

