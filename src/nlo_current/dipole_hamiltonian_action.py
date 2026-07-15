"""Direct observable-side NLO sector actions on the fundamental dipole."""

from __future__ import annotations

import numpy as np

from .dipole_observable import apply_generator_word_to_dipole
from .three_generator_terms import (
    kjjsj_A_LLR_from_kernel,
    kjjsj_B_LRR_from_kernel,
    kjjsj_V_virtual_from_kernel,
    kjjssj_A_LLR_from_kernel,
    kjjssj_B_LRR_from_kernel,
    kjjssj_V_virtual_from_kernel,
)
from .two_generator_terms import kjssj_A_from_kernel, kqbarq_A_from_kernel


NC = 3


def _n_color(gens):
    return int(gens.shape[0])


def action_KJSJ_direct(U_fund, S_adj, KJSJ, u, v, gens):
    """Return -H_KJSJ s(u,v) by direct observable-side generator action."""

    nsite = np.asarray(U_fund).shape[0]
    n_color = _n_color(gens)
    total = 0.0j
    for x in range(nsite):
        for y in range(nsite):
            for z in range(nsite):
                kernel = KJSJ[x, y, z]
                if kernel == 0.0:
                    continue
                block = 0.0j
                for a in range(n_color):
                    block += apply_generator_word_to_dipole(
                        U_fund,
                        ((x, a), (y, a)),
                        u,
                        v,
                        gens,
                        ("L", "L"),
                    )
                    block += apply_generator_word_to_dipole(
                        U_fund,
                        ((x, a), (y, a)),
                        u,
                        v,
                        gens,
                        ("R", "R"),
                    )
                    for b in range(n_color):
                        block += -2.0 * S_adj[z, a, b] * apply_generator_word_to_dipole(
                            U_fund,
                            ((x, a), (y, b)),
                            u,
                            v,
                            gens,
                            ("L", "R"),
                        )
                total += -kernel * block
    return total


def _action_LR_from_A(U_fund, A, u, v, gens):
    nsite = np.asarray(U_fund).shape[0]
    n_color = _n_color(gens)
    total = 0.0j
    for x in range(nsite):
        for y in range(nsite):
            for a in range(n_color):
                for b in range(n_color):
                    coeff = A[x, y, a, b]
                    if coeff == 0.0:
                        continue
                    total += -coeff * apply_generator_word_to_dipole(
                        U_fund,
                        ((x, a), (y, b)),
                        u,
                        v,
                        gens,
                        ("L", "R"),
                    )
    return total


def action_KJSSJ_direct(U_fund, S_adj, KJSSJ, u, v, f, gens):
    """Return -H_KJSSJ s(u,v) using the ordered KLM two-generator block."""

    A = kjssj_A_from_kernel(S_adj, KJSSJ, f)
    return _action_LR_from_A(U_fund, A, u, v, gens)


def action_Kqbarq_direct(U_fund, S_adj, Kqbarq, u, v, gens):
    """Return -H_Kqbarq s(u,v) by direct observable-side generator action."""

    A = kqbarq_A_from_kernel(U_fund, S_adj, Kqbarq, gens)
    return _action_LR_from_A(U_fund, A, u, v, gens)


def _action_cubic_from_blocks(U_fund, A_LLR, B_LRR, V, u, v, gens, virtual_scale=1.0):
    nsite = np.asarray(U_fund).shape[0]
    n_color = _n_color(gens)
    total = 0.0j

    for x in range(nsite):
        for y in range(nsite):
            for w in range(nsite):
                for d in range(n_color):
                    for e in range(n_color):
                        for a in range(n_color):
                            coeff = A_LLR[x, y, w, d, e, a]
                            if coeff != 0.0:
                                total += -coeff * apply_generator_word_to_dipole(
                                    U_fund,
                                    ((x, d), (y, e), (w, a)),
                                    u,
                                    v,
                                    gens,
                                    ("L", "L", "R"),
                                )

                            coeff = B_LRR[w, x, y, a, d, e]
                            if coeff != 0.0:
                                total += -coeff * apply_generator_word_to_dipole(
                                    U_fund,
                                    ((w, a), (x, d), (y, e)),
                                    u,
                                    v,
                                    gens,
                                    ("L", "R", "R"),
                                )

                            v_coeff = virtual_scale * V[x, y, w, d, e, a]
                            if v_coeff != 0.0:
                                total += -v_coeff * apply_generator_word_to_dipole(
                                    U_fund,
                                    ((x, d), (y, e), (w, a)),
                                    u,
                                    v,
                                    gens,
                                    ("L", "L", "L"),
                                )
                                total += v_coeff * apply_generator_word_to_dipole(
                                    U_fund,
                                    ((x, d), (y, e), (w, a)),
                                    u,
                                    v,
                                    gens,
                                    ("R", "R", "R"),
                                )
    return total


def action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens, virtual_scale=1.0):
    """Return -H_KJJSJ s(u,v) for LLR, LRR, and virtual pieces."""

    A = kjjsj_A_LLR_from_kernel(S_adj, KJJSJ, f)
    B = kjjsj_B_LRR_from_kernel(S_adj, KJJSJ, f)
    V = kjjsj_V_virtual_from_kernel(KJJSJ, f)
    return _action_cubic_from_blocks(U_fund, A, B, V, u, v, gens, virtual_scale)


def action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens, virtual_scale=1.0):
    """Return -H_KJJSSJ s(u,v) for LLR, LRR, and virtual pieces."""

    A = kjjssj_A_LLR_from_kernel(S_adj, KJJSSJ, f)
    B = kjjssj_B_LRR_from_kernel(S_adj, KJJSSJ, f)
    V = kjjssj_V_virtual_from_kernel(KJJSSJ, f)
    return _action_cubic_from_blocks(U_fund, A, B, V, u, v, gens, virtual_scale)


def action_all_sectors_direct(U_fund, S_adj, kernels, u, v, f, gens):
    """Return direct actions by sector and their total."""

    actions = {}
    if "KJSJ" in kernels:
        actions["KJSJ"] = action_KJSJ_direct(U_fund, S_adj, kernels["KJSJ"], u, v, gens)
    if "KJSSJ" in kernels:
        actions["KJSSJ"] = action_KJSSJ_direct(U_fund, S_adj, kernels["KJSSJ"], u, v, f, gens)
    if "Kqbarq" in kernels:
        actions["Kqbarq"] = action_Kqbarq_direct(U_fund, S_adj, kernels["Kqbarq"], u, v, gens)
    if "KJJSJ" in kernels:
        actions["KJJSJ"] = action_KJJSJ_direct(U_fund, S_adj, kernels["KJJSJ"], u, v, f, gens)
    if "KJJSSJ" in kernels:
        actions["KJJSSJ"] = action_KJJSSJ_direct(
            U_fund,
            S_adj,
            kernels["KJJSSJ"],
            u,
            v,
            f,
            gens,
        )
    actions["total"] = sum(actions.values(), 0.0j)
    return actions
