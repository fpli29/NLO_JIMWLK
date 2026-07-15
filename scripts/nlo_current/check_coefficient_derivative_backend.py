from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.coefficient_derivatives import (  # noqa: E402
    compute_all_coefficient_derivatives_fd,
    product_rule_K2_rhs,
    product_rule_K3_rhs,
    velocity_from_coeff_derivative_backend,
)
from nlo_current.finite_difference_scores import (  # noqa: E402
    fd_left_derivative_scalar,
    fd_left_second_derivative_scalar,
)
from nlo_current.nlo_current_skeleton import NLOCurrentTerms  # noqa: E402
from nlo_current.nlo_velocity_evaluator import evaluate_velocity_from_terms  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    su3_generators_fundamental,
)


SEED = 20260705
NSITE = 1
EPS_FIRST = 2e-5
EPS_SECOND = 7e-4


def _unflatten(index: int, n_color: int) -> tuple[int, int]:
    return index // n_color, index % n_color


def _S_builder(gens):
    return lambda U_fund: np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])


def _score_and_hessian(logW, U_fund, gens):
    dim = len(U_fund) * gens.shape[0]
    U_list = [np.array(U, copy=True) for U in U_fund]
    score = np.zeros(dim)
    hessian = np.zeros((dim, dim))
    for b in range(dim):
        site_b, color_b = _unflatten(b, gens.shape[0])
        score[b] = fd_left_derivative_scalar(logW, U_list, site_b, color_b, gens, EPS_FIRST)
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
                EPS_SECOND,
            )
    return score, hessian


def _K2_product_rule_residual(K2_callback, K2, dK2, score, logW, U_fund, S_builder, gens):
    dim = K2.shape[0]
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
                EPS_FIRST,
            )
    lhs /= W0
    rhs = product_rule_K2_rhs(K2, dK2, score)
    return float(np.linalg.norm(lhs - rhs))


def _selected_K3_product_rule_residuals(
    K3_callback,
    K3,
    dK3_first,
    d2K3,
    score,
    hessian,
    logW,
    U_fund,
    S_builder,
    gens,
):
    dim = K3.shape[0]
    rhs = product_rule_K3_rhs(K3, dK3_first, d2K3, score, hessian)
    U_list = [np.array(U, copy=True) for U in U_fund]
    W0 = np.exp(logW(U_list))
    residuals = {}
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
                    EPS_SECOND,
                )
        lhs_a /= W0
        residuals[a] = float(abs(lhs_a - rhs[a]))
    return residuals


def main() -> None:
    rng = np.random.default_rng(SEED)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng) for _ in range(NSITE)])
    S_builder = _S_builder(gens)
    dim = NSITE * gens.shape[0]

    Q_coeff = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    Q_density = 0.08 * (rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3)))
    M2 = 0.2 * rng.normal(size=(dim, dim))
    M3 = 0.04 * rng.normal(size=(dim, dim, dim))

    def phi(U):
        return float(np.real(np.trace(Q_coeff @ U[0])))

    def K2_callback(U, _S):
        return phi(U) * M2

    def K3_callback(U, _S):
        return phi(U) * M3

    def logW(Vs):
        return float(np.real(np.trace(Q_density @ Vs[0])))

    derivatives = compute_all_coefficient_derivatives_fd(
        K2_callback,
        K3_callback,
        U_fund,
        S_builder,
        gens,
        eps_first=EPS_FIRST,
        eps_second=EPS_SECOND,
    )
    K2 = K2_callback(U_fund, S_builder(U_fund))
    K3 = K3_callback(U_fund, S_builder(U_fund))
    terms = NLOCurrentTerms(K1=np.zeros(dim), K2=K2, K3=K3, metadata={})
    score, hessian = _score_and_hessian(logW, U_fund, gens)
    with_derivatives = velocity_from_coeff_derivative_backend(terms, score, hessian, derivatives)
    without_derivatives = evaluate_velocity_from_terms(terms, score, hessian)

    K2_residual = _K2_product_rule_residual(
        K2_callback,
        K2,
        derivatives["dK2"],
        score,
        logW,
        U_fund,
        S_builder,
        gens,
    )
    K3_residuals = _selected_K3_product_rule_residuals(
        K3_callback,
        K3,
        derivatives["dK3_first"],
        derivatives["d2K3"],
        score,
        hessian,
        logW,
        U_fund,
        S_builder,
        gens,
    )

    velocity_delta = with_derivatives["velocity"] - without_derivatives["velocity"]
    report = ROOT / "reports" / "nlo_current" / "coefficient_derivative_backend_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# Coefficient-Derivative Backend Report\n\n"
        "This report was generated by `scripts/nlo_current/check_coefficient_derivative_backend.py`.\n\n"
        "## Run Configuration\n\n"
        f"- random_seed: {SEED}\n"
        f"- nsite: {NSITE}\n"
        f"- combined_dimension_D: {dim}\n"
        f"- eps_first: {EPS_FIRST:.16e}\n"
        f"- eps_second: {EPS_SECOND:.16e}\n"
        "- backend_scope: dense finite-difference diagnostic only\n\n"
        "## Derivative Norms\n\n"
        f"- dK2_norm: {np.linalg.norm(derivatives['dK2']):.16e}\n"
        f"- LC_K3_norm: {np.linalg.norm(derivatives['dK3_first']['LC_K3_ABC']):.16e}\n"
        f"- LB_K3_norm: {np.linalg.norm(derivatives['dK3_first']['LB_K3_ABC']):.16e}\n"
        f"- d2K3_norm: {np.linalg.norm(derivatives['d2K3']):.16e}\n\n"
        "## Velocity Comparison\n\n"
        f"- velocity_with_derivatives_norm: {np.linalg.norm(with_derivatives['velocity']):.16e}\n"
        f"- velocity_without_derivatives_norm: {np.linalg.norm(without_derivatives['velocity']):.16e}\n"
        f"- velocity_difference_norm: {np.linalg.norm(velocity_delta):.16e}\n"
        f"- omitted_derivative_warnings_without: {len(without_derivatives['diagnostics']['warnings'])}\n"
        f"- omitted_derivative_warnings_with: {len(with_derivatives['diagnostics']['warnings'])}\n\n"
        "## Product-Rule Residuals\n\n"
        f"- K2_all_indices_residual_norm: {K2_residual:.16e}\n"
        "- K3_selected_indices_absolute_residuals:\n"
        + "".join(f"  - A={a}: {residual:.16e}\n" for a, residual in K3_residuals.items())
        + "\n## Warning\n\n"
        "This backend is finite-difference diagnostic infrastructure only. It is not suitable "
        "for production NLO evolution or realistic lattice sizes.\n",
        encoding="utf-8",
    )

    print(f"wrote {report}")
    print(f"dK2_norm={np.linalg.norm(derivatives['dK2']):.16e}")
    print(f"LC_K3_norm={np.linalg.norm(derivatives['dK3_first']['LC_K3_ABC']):.16e}")
    print(f"LB_K3_norm={np.linalg.norm(derivatives['dK3_first']['LB_K3_ABC']):.16e}")
    print(f"d2K3_norm={np.linalg.norm(derivatives['d2K3']):.16e}")
    print(f"velocity_difference_norm={np.linalg.norm(velocity_delta):.16e}")


if __name__ == "__main__":
    main()
