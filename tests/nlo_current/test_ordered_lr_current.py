from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    left_perturb,
    random_su3,
    right_perturb,
    su3_generators_fundamental,
)


def test_ordered_lr_current_matches_left_divergence_form() -> None:
    rng = np.random.default_rng(7890)
    gens = su3_generators_fundamental()
    Us = [random_su3(rng), random_su3(rng)]
    base = 0.05 * rng.normal(size=(8, 8))
    color_weight = rng.normal(size=(8, 8))
    probe0 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    probe1 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    q0 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    q1 = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))

    def signal(Vs):
        return float(np.real(np.trace(probe0 @ Vs[0])) + 0.4 * np.real(np.trace(probe1 @ Vs[1])))

    def coefficient(Vs, a, b):
        sig = signal(Vs)
        return float(base[a, b] + 0.015 * color_weight[a, b] * sig + 0.01 * np.cos(sig))

    def W_toy(Vs):
        exponent = 0.04 * np.real(np.trace(q0 @ Vs[0])) + 0.03 * np.real(np.trace(q1 @ Vs[1]))
        return float(np.exp(exponent))

    def left_derivative_re_trace(matrix, U, a):
        return float(np.real(1.0j * np.trace(matrix @ gens[a] @ U)))

    def analytic_left_x_F(Vs, a, b):
        sig = signal(Vs)
        d_sig = left_derivative_re_trace(probe0, Vs[0], a)
        d_coeff = (0.015 * color_weight[a, b] - 0.01 * np.sin(sig)) * d_sig
        d_log_w = 0.04 * left_derivative_re_trace(q0, Vs[0], a)
        coeff = coefficient(Vs, a, b)
        return float(W_toy(Vs) * (d_coeff + coeff * d_log_w))

    def forward_left_derivative(func, Vs, site, a, eps):
        plus = [np.array(U, copy=True) for U in Vs]
        plus[site] = left_perturb(plus[site], a, eps, gens)
        return (func(plus) - func(Vs)) / eps

    def forward_right_derivative(func, Vs, site, a, eps):
        plus = [np.array(U, copy=True) for U in Vs]
        plus[site] = right_perturb(plus[site], a, eps, gens)
        return (func(plus) - func(Vs)) / eps

    def direct_lr_density(Vs, eps):
        total = 0.0
        for a in range(8):
            for b in range(8):
                total -= forward_right_derivative(
                    lambda Ws, a=a, b=b: analytic_left_x_F(Ws, a, b),
                    Vs,
                    1,
                    b,
                    eps,
                )
        return float(total)

    def divergence_current_density(Vs, eps):
        total = 0.0
        for h in range(8):
            def current_h(Ws, h=h):
                S_y = adjoint_from_fundamental(Ws[1], gens)
                subtotal = 0.0
                for b in range(8):
                    inner = 0.0
                    for a in range(8):
                        inner += analytic_left_x_F(Ws, a, b)
                    subtotal += S_y[h, b] * inner
                return float(subtotal)

            total -= forward_left_derivative(current_h, Vs, 1, h, eps)
        return float(total)

    residuals = []
    for eps in (1e-4, 1e-5, 1e-6):
        direct = direct_lr_density(Us, eps)
        divergence = divergence_current_density(Us, eps)
        residual = abs(direct - divergence) / max(abs(direct), abs(divergence), 1.0)
        residuals.append((eps, direct, divergence, residual))

    report_dir = ROOT / "reports" / "nlo_current"
    report_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Ordered LR Current Finite-Difference Check",
        "",
        "Checks -J_R^b J_L^a(A^{ab}W) against "
        "-L_y^h{S_y^{hb} L_x^a[A^{ab}W]} on a two-site SU(3) toy problem. "
        "The inner L_x derivative is evaluated analytically for the toy A and W; "
        "the outer left/right derivatives are forward finite differences.",
        "",
        "seed: 7890",
        "",
        "| eps | direct | divergence | relative residual |",
        "|---:|---:|---:|---:|",
    ]
    for eps, direct, divergence, residual in residuals:
        lines.append(f"| {eps:.0e} | {direct:.16e} | {divergence:.16e} | {residual:.16e} |")
    lines.append("")
    (report_dir / "ordered_lr_current_fd_report.md").write_text("\n".join(lines), encoding="utf-8")

    assert residuals[1][3] < residuals[0][3]
    assert residuals[2][3] < residuals[1][3]
    assert residuals[2][3] < 2e-8
