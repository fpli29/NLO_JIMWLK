"""Positive tiny-lattice test densities for density-side closure checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .finite_difference_scores import fd_left_derivative_scalar, fd_left_second_derivative_scalar


@dataclass(frozen=True)
class TestDensityResult:
    """Value, score, and ordered Hessian-score for a positive test density."""

    log_weight: float
    weight: float
    score: np.ndarray
    hessian_score: np.ndarray
    metadata: dict


def _unflatten(index: int, n_color: int) -> tuple[int, int]:
    return int(index) // int(n_color), int(index) % int(n_color)


def default_density_params(name: str) -> dict:
    """Return deterministic parameters for a named positive test density."""

    if name == "single_link_trace":
        return {"lambda1": 0.07, "x": 0}
    if name == "dipole_trace":
        return {"lambda1": 0.05, "lambda2": -0.04, "x": 0, "y": 1}
    if name == "multilink_nonlinear":
        return {"lambda1": 0.04, "lambda2": -0.03, "lambda3": 0.015, "x": 0, "y": 1, "z": 2}
    if name == "constant":
        return {}
    raise ValueError(f"unknown test density: {name}")


def density_log_weight(U_fund, name: str, params: dict | None = None) -> float:
    """Return log W for one of the positive diagnostic densities."""

    U = np.asarray(U_fund)
    p = default_density_params(name)
    if params:
        p.update(params)

    if name == "constant":
        return 0.0

    x = int(p.get("x", 0))
    if x < 0 or x >= len(U):
        raise ValueError("density site x is outside U_fund")
    total = float(p.get("lambda1", 0.0)) * float(np.real(np.trace(U[x])))

    if name in {"dipole_trace", "multilink_nonlinear"}:
        y = int(p.get("y", 1))
        if y < 0 or y >= len(U):
            raise ValueError("density site y is outside U_fund")
        total += float(p.get("lambda2", 0.0)) * float(np.real(np.trace(U[x].conj().T @ U[y])))

    if name == "multilink_nonlinear":
        z = int(p.get("z", 2))
        y = int(p.get("y", 1))
        if z < 0 or z >= len(U):
            raise ValueError("density site z is outside U_fund")
        yz_trace = float(np.real(np.trace(U[y].conj().T @ U[z])))
        total += float(p.get("lambda3", 0.0)) * yz_trace**2

    return float(total)


def compute_test_density_score(
    U_fund,
    gens,
    name: str,
    params: dict | None = None,
    *,
    eps: float = 2.0e-5,
) -> np.ndarray:
    """Compute score s_A = L_A log W by central finite differences."""

    U_list = [np.array(U, copy=True) for U in np.asarray(U_fund)]
    n_color = int(gens.shape[0])
    dim = len(U_list) * n_color

    def logW(Vs):
        return density_log_weight(np.stack(Vs), name, params)

    score = np.zeros(dim, dtype=float)
    for index in range(dim):
        site, color = _unflatten(index, n_color)
        score[index] = fd_left_derivative_scalar(logW, U_list, site, color, gens, eps)
    return score


def compute_test_density_hessian_score(
    U_fund,
    gens,
    name: str,
    params: dict | None = None,
    *,
    eps: float = 5.0e-4,
) -> np.ndarray:
    """Compute ordered Hessian-score H_AB = L_A s_B by central differences."""

    U_list = [np.array(U, copy=True) for U in np.asarray(U_fund)]
    n_color = int(gens.shape[0])
    dim = len(U_list) * n_color

    def logW(Vs):
        return density_log_weight(np.stack(Vs), name, params)

    hessian = np.zeros((dim, dim), dtype=float)
    for a in range(dim):
        site_a, color_a = _unflatten(a, n_color)
        for b in range(dim):
            site_b, color_b = _unflatten(b, n_color)
            hessian[a, b] = fd_left_second_derivative_scalar(
                logW,
                U_list,
                site_a,
                color_a,
                site_b,
                color_b,
                gens,
                eps,
            )
    return hessian


def evaluate_test_density(
    U_fund,
    gens,
    name: str,
    params: dict | None = None,
    *,
    score_eps: float = 2.0e-5,
    hessian_eps: float = 5.0e-4,
) -> TestDensityResult:
    """Evaluate a positive diagnostic density and its derivative data."""

    log_weight = density_log_weight(U_fund, name, params)
    score = compute_test_density_score(U_fund, gens, name, params, eps=score_eps)
    hessian = compute_test_density_hessian_score(U_fund, gens, name, params, eps=hessian_eps)
    return TestDensityResult(
        log_weight=log_weight,
        weight=float(np.exp(log_weight)),
        score=score,
        hessian_score=hessian,
        metadata={
            "density": name,
            "params": default_density_params(name) | dict(params or {}),
            "derivative_backend": "finite_difference",
            "score_eps": float(score_eps),
            "hessian_eps": float(hessian_eps),
            "positive_by_construction": True,
        },
    )

