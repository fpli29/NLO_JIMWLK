"""Fundamental dipole observable and generator actions for tiny diagnostics."""

from __future__ import annotations

import numpy as np

from .su3_adjoint import left_perturb, right_perturb


def dipole(U_fund, u, v):
    """Return s(u,v) = (1/Nc) Tr(U_u^dagger U_v)."""

    U_fund = np.asarray(U_fund)
    nc = U_fund.shape[-1]
    return np.trace(U_fund[u].conj().T @ U_fund[v]) / nc


def left_generator_action_on_dipole(U_fund, site, color, u, v, gens):
    """Analytic J_L action on s(u,v) for L F(U)=dF(exp(i eps t)U)."""

    U_fund = np.asarray(U_fund)
    nc = U_fund.shape[-1]
    t = gens[color]
    out = 0.0j
    if site == u:
        out += -1.0j * np.trace(U_fund[u].conj().T @ t @ U_fund[v]) / nc
    if site == v:
        out += 1.0j * np.trace(U_fund[u].conj().T @ t @ U_fund[v]) / nc
    return out


def right_generator_action_on_dipole(U_fund, site, color, u, v, gens):
    """Analytic J_R action on s(u,v) for R F(U)=dF(U exp(i eps t))."""

    U_fund = np.asarray(U_fund)
    nc = U_fund.shape[-1]
    t = gens[color]
    out = 0.0j
    if site == u:
        out += -1.0j * np.trace(t @ U_fund[u].conj().T @ U_fund[v]) / nc
    if site == v:
        out += 1.0j * np.trace(U_fund[u].conj().T @ U_fund[v] @ t) / nc
    return out


def _initial_terms(U_fund, u, v):
    nc = U_fund.shape[-1]
    return [(1.0 / nc, [("Ud", u), ("U", v)])]


def _differentiate_factor(factor, side, site, color, gens):
    kind, factor_site = factor
    if kind == "M":
        return []
    if factor_site != site:
        return []
    t = np.asarray(gens[color])
    if kind == "Ud" and side == "L":
        return [(-1.0j, [("Ud", factor_site), ("M", t)])]
    if kind == "Ud" and side == "R":
        return [(-1.0j, [("M", t), ("Ud", factor_site)])]
    if kind == "U" and side == "L":
        return [(1.0j, [("M", t), ("U", factor_site)])]
    if kind == "U" and side == "R":
        return [(1.0j, [("U", factor_site), ("M", t)])]
    return []


def _apply_one_generator_to_terms(terms, side, site, color, gens):
    out = []
    for coeff, factors in terms:
        for idx, factor in enumerate(factors):
            for d_coeff, replacement in _differentiate_factor(factor, side, site, color, gens):
                new_factors = factors[:idx] + replacement + factors[idx + 1 :]
                out.append((coeff * d_coeff, new_factors))
    return out


def _evaluate_terms(terms, U_fund):
    total = 0.0j
    for coeff, factors in terms:
        product = np.eye(U_fund.shape[-1], dtype=complex)
        for kind, value in factors:
            if kind == "Ud":
                product = product @ U_fund[value].conj().T
            elif kind == "U":
                product = product @ U_fund[value]
            elif kind == "M":
                product = product @ value
            else:
                raise ValueError(f"unknown factor kind: {kind}")
        total += coeff * np.trace(product)
    return total


def apply_generator_word_to_dipole(U_fund, word, u, v, gens, side_labels):
    """Apply an ordered generator word to the dipole.

    The tuple order is the operator order as written: for word (A,B,C), the
    action is A B C s, so C acts first.
    """

    if len(word) != len(side_labels):
        raise ValueError("word and side_labels must have the same length")
    terms = _initial_terms(np.asarray(U_fund), u, v)
    for (site, color), side in reversed(tuple(zip(word, side_labels))):
        if side not in {"L", "R"}:
            raise ValueError("side labels must be 'L' or 'R'")
        terms = _apply_one_generator_to_terms(terms, side, site, color, gens)
        if not terms:
            return 0.0j
    return _evaluate_terms(terms, np.asarray(U_fund))


def fd_generator_action_on_dipole(U_fund, site, color, u, v, gens, side="L", eps=1e-6):
    """Finite-difference check for one generator action."""

    U_fund = np.asarray(U_fund)
    plus = np.array(U_fund, copy=True)
    minus = np.array(U_fund, copy=True)
    if side == "L":
        plus[site] = left_perturb(plus[site], color, eps, gens)
        minus[site] = left_perturb(minus[site], color, -eps, gens)
    elif side == "R":
        plus[site] = right_perturb(plus[site], color, eps, gens)
        minus[site] = right_perturb(minus[site], color, -eps, gens)
    else:
        raise ValueError("side must be 'L' or 'R'")
    return (dipole(plus, u, v) - dipole(minus, u, v)) / (2.0 * eps)


def _perturb_copy(U_fund, site, color, gens, side, eps):
    out = np.array(U_fund, copy=True)
    if side == "L":
        out[site] = left_perturb(out[site], color, eps, gens)
    elif side == "R":
        out[site] = right_perturb(out[site], color, eps, gens)
    else:
        raise ValueError("side must be 'L' or 'R'")
    return out


def fd_generator_word_action_on_dipole(U_fund, word, u, v, gens, side_labels, eps=1e-5):
    """Nested finite-difference action for small validation only."""

    if len(word) != len(side_labels):
        raise ValueError("word and side_labels must have the same length")

    def action(U_current, index):
        if index == len(word):
            return dipole(U_current, u, v)
        site, color = word[index]
        side = side_labels[index]
        plus = _perturb_copy(U_current, site, color, gens, side, eps)
        minus = _perturb_copy(U_current, site, color, gens, side, -eps)
        return (action(plus, index + 1) - action(minus, index + 1)) / (2.0 * eps)

    return action(np.asarray(U_fund), 0)
