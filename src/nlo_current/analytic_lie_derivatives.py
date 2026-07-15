"""Analytic local Lie-derivative primitives for tiny SU(3) diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TraceFactor:
    """One factor in a fundamental trace word."""

    kind: str
    site: int | None = None
    matrix: np.ndarray | None = None


def _zero_like_matrix(dtype=complex):
    return np.zeros((3, 3), dtype=dtype)


def left_derivative_fundamental(U, gens, site: int, color: int, target_site: int) -> np.ndarray:
    r"""Return \(L_site^color U_target_site\).

    The convention is
    \(L_x^a F(U)=dF(e^{i eps t^a}U_x)/d eps|_0\), hence
    \(L_x^a U_y = i delta_{xy} t^a U_y\).
    """

    U_arr = np.asarray(U)
    if int(site) != int(target_site):
        return _zero_like_matrix(np.result_type(U_arr, complex))
    return 1.0j * gens[color] @ U_arr[target_site]


def left_derivative_fundamental_dagger(U, gens, site: int, color: int, target_site: int) -> np.ndarray:
    r"""Return \(L_site^color U_target_site^\dagger\).

    With the left-perturbation convention,
    \(L_x^a U_y^\dagger=-i delta_{xy} U_y^\dagger t^a\).
    """

    U_arr = np.asarray(U)
    if int(site) != int(target_site):
        return _zero_like_matrix(np.result_type(U_arr, complex))
    return -1.0j * U_arr[target_site].conj().T @ gens[color]


def left_derivative_adjoint(
    S_adj,
    f,
    site: int,
    color: int,
    target_site: int,
) -> np.ndarray:
    r"""Return \(L_site^color S_A(target_site)\).

    For \(S_A^{ab}=2 Re Tr(t^a U t^b U^\dagger)\), the calibrated rule is
    \(L^h S_A^{ab}=f^{hac}S_A^{cb}\). The derivative acts on the first
    adjoint index in this convention.
    """

    S = np.asarray(S_adj)
    if int(site) != int(target_site):
        return np.zeros_like(S[target_site], dtype=np.result_type(S, f, float))
    return np.einsum("ac,cb->ab", f[color], S[target_site], optimize=True)


def second_left_derivative_adjoint(
    S_adj,
    f,
    first_site: int,
    first_color: int,
    second_site: int,
    second_color: int,
    target_site: int,
) -> np.ndarray:
    r"""Return ordered \(L_first L_second S_A(target_site)\).

    If all sites coincide,
    \(L^g L^h S^{ab}=f^{hac}f^{gcd}S^{db}\), preserving the written order.
    """

    S = np.asarray(S_adj)
    if int(first_site) != int(target_site) or int(second_site) != int(target_site):
        return np.zeros_like(S[target_site], dtype=np.result_type(S, f, float))
    return np.einsum(
        "ac,cd,db->ab",
        f[second_color],
        f[first_color],
        S[target_site],
        optimize=True,
    )


def _factor_value(factor: TraceFactor, U) -> np.ndarray:
    if factor.kind == "U":
        return np.asarray(U)[factor.site]
    if factor.kind == "Udag":
        return np.asarray(U)[factor.site].conj().T
    if factor.kind == "matrix":
        return np.asarray(factor.matrix)
    raise ValueError(f"unknown trace factor kind: {factor.kind}")


def _factor_derivative(factor: TraceFactor, U, gens, site: int, color: int) -> np.ndarray:
    if factor.kind == "U":
        return left_derivative_fundamental(U, gens, site, color, factor.site)
    if factor.kind == "Udag":
        return left_derivative_fundamental_dagger(U, gens, site, color, factor.site)
    if factor.kind == "matrix":
        return np.zeros_like(np.asarray(factor.matrix), dtype=np.result_type(factor.matrix, complex))
    raise ValueError(f"unknown trace factor kind: {factor.kind}")


def _factor_second_derivative(
    factor: TraceFactor,
    U,
    gens,
    first_site: int,
    first_color: int,
    second_site: int,
    second_color: int,
) -> np.ndarray:
    if first_site != second_site or factor.site != first_site:
        return np.zeros((3, 3), dtype=complex)
    U_arr = np.asarray(U)
    if factor.kind == "U":
        return -gens[first_color] @ gens[second_color] @ U_arr[factor.site]
    if factor.kind == "Udag":
        return -U_arr[factor.site].conj().T @ gens[first_color] @ gens[second_color]
    if factor.kind == "matrix":
        return np.zeros_like(np.asarray(factor.matrix), dtype=np.result_type(factor.matrix, complex))
    raise ValueError(f"unknown trace factor kind: {factor.kind}")


def trace_word_value(U, factors: list[TraceFactor]) -> complex:
    """Return ``Tr(prod factors)`` for a fundamental trace word."""

    product = np.eye(3, dtype=complex)
    for factor in factors:
        product = product @ _factor_value(factor, U)
    return complex(np.trace(product))


def left_derivative_trace_word(U, gens, factors: list[TraceFactor], site: int, color: int) -> complex:
    """Return the product-rule derivative of an explicit trace word."""

    total = 0.0 + 0.0j
    values = [_factor_value(factor, U) for factor in factors]
    for i, factor in enumerate(factors):
        deriv = _factor_derivative(factor, U, gens, site, color)
        if not np.any(deriv):
            continue
        product = np.eye(3, dtype=complex)
        for j, value in enumerate(values):
            product = product @ (deriv if i == j else value)
        total += np.trace(product)
    return complex(total)


def second_left_derivative_trace_word(
    U,
    gens,
    factors: list[TraceFactor],
    first_site: int,
    first_color: int,
    second_site: int,
    second_color: int,
) -> complex:
    r"""Return ordered product-rule \(L_first L_second Tr(word)\)."""

    total = 0.0 + 0.0j
    values = [_factor_value(factor, U) for factor in factors]
    first_derivs = [
        _factor_derivative(factor, U, gens, first_site, first_color) for factor in factors
    ]
    second_derivs = [
        _factor_derivative(factor, U, gens, second_site, second_color) for factor in factors
    ]
    double_derivs = [
        _factor_second_derivative(
            factor,
            U,
            gens,
            first_site,
            first_color,
            second_site,
            second_color,
        )
        for factor in factors
    ]

    for i, double in enumerate(double_derivs):
        if np.any(double):
            product = np.eye(3, dtype=complex)
            for k, value in enumerate(values):
                product = product @ (double if i == k else value)
            total += np.trace(product)

    for i, first in enumerate(first_derivs):
        if not np.any(first):
            continue
        for j, second in enumerate(second_derivs):
            if i == j or not np.any(second):
                continue
            product = np.eye(3, dtype=complex)
            for k, value in enumerate(values):
                if k == i:
                    product = product @ first
                elif k == j:
                    product = product @ second
                else:
                    product = product @ value
            total += np.trace(product)

    return complex(total)
