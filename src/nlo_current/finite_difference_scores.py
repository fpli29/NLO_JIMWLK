"""Finite-difference score and Hessian-score helpers for tiny SU(3) tests."""

from __future__ import annotations

import numpy as np

from .su3_adjoint import left_perturb


def toy_log_density(U_list: list[np.ndarray], params: dict) -> float:
    """Return a coupled toy log density with nonzero Hessian-score components."""

    q_mats = params["q_mats"]
    lam = float(params.get("lambda", 0.1))
    eta = float(params.get("eta", 0.05))

    total = 0.0
    for U, Q in zip(U_list, q_mats):
        total += lam * np.real(np.trace(Q @ U))

    for i in range(len(U_list)):
        for j in range(i + 1, len(U_list)):
            total += eta * np.real(np.trace(U_list[i] @ U_list[j].conj().T))
    return float(total)


def fd_left_derivative_scalar(func, U_list, site, color, gens, eps):
    """Central finite-difference left derivative of a scalar function."""

    plus = [np.array(U, copy=True) for U in U_list]
    minus = [np.array(U, copy=True) for U in U_list]
    plus[site] = left_perturb(plus[site], color, eps, gens)
    minus[site] = left_perturb(minus[site], color, -eps, gens)
    return float((func(plus) - func(minus)) / (2.0 * eps))


def fd_left_second_derivative_scalar(func, U_list, site_a, color_a, site_b, color_b, gens, eps):
    """Central finite-difference ordered second left derivative L_A L_B func."""

    return fd_left_derivative_scalar(
        lambda Vs: fd_left_derivative_scalar(func, Vs, site_b, color_b, gens, eps),
        U_list,
        site_a,
        color_a,
        gens,
        eps,
    )


def fd_score(logW_func, U_list, site, color, gens, eps):
    """Finite-difference score s_A = L_A log W."""

    return fd_left_derivative_scalar(logW_func, U_list, site, color, gens, eps)


def fd_hessian_score(logW_func, U_list, site_a, color_a, site_b, color_b, gens, eps):
    """Finite-difference Hessian-score H_AB = L_A s_B = L_A L_B log W."""

    return fd_left_second_derivative_scalar(
        logW_func,
        U_list,
        site_a,
        color_a,
        site_b,
        color_b,
        gens,
        eps,
    )

