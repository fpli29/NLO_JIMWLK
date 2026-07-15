"""KLM Appendix A dipole targets available in local validated notes."""

from __future__ import annotations

import numpy as np

from .dipole_observable import dipole
from .su3_adjoint import su3_generators_fundamental


NC = 3

_AVAILABLE = {
    "KJSJ": True,
    "KJSSJ": True,
    "Kqbarq": True,
    "KJJSJ": True,
    "KJJSSJ": True,
}


def target_KJSJ_appendix(U_fund, KJSJ, u, v):
    """Appendix A target for H_KJSJ s(u,v)."""

    total = 0.0j
    for z in range(np.asarray(U_fund).shape[0]):
        total += KJSJ[u, v, z] * (dipole(U_fund, u, z) * dipole(U_fund, z, v) - dipole(U_fund, u, v))
    return 2.0 * NC * total


def _trace_word(U_fund, factors):
    product = np.eye(np.asarray(U_fund).shape[-1], dtype=complex)
    for dagger, site in factors:
        U_site = U_fund[site]
        product = product @ (U_site.conj().T if dagger else U_site)
    return np.trace(product)


def _kqbarq_endpoint_combo(Kqbarq, u, v, z, zp):
    return Kqbarq[v, v, z, zp] - Kqbarq[u, v, z, zp] - Kqbarq[v, u, z, zp] + Kqbarq[u, u, z, zp]


def _generator_trace_product(U_fund, left, right, gens):
    return np.einsum(
        "ij,ajk,kl,bli->ab",
        U_fund[left].conj().T,
        gens,
        U_fund[right],
        gens,
        optimize=True,
    )


def target_Kqbarq_trace_current_appendix(U_fund, Kqbarq, u, v, gens=None):
    """Appendix A subsection target for the trace-current part only.

    This implements the exact generator trace-product expression in
    WORKNLO.tex lines 1174--1177. The later compact reduced line is not used as
    the implementation source because finite-grid endpoint/contact policies can
    make that reduction ambiguous in diagnostics.
    """

    U_fund = np.asarray(U_fund)
    gens = su3_generators_fundamental() if gens is None else np.asarray(gens)
    total = 0.0j
    nsite = U_fund.shape[0]
    dipole_trace = _generator_trace_product(U_fund, u, v, gens)
    for z in range(nsite):
        for zp in range(nsite):
            kernel_combo = _kqbarq_endpoint_combo(Kqbarq, u, v, z, zp)
            qbarq_trace = _generator_trace_product(U_fund, z, zp, gens)
            total += (kernel_combo / NC) * 2.0 * np.einsum("ab,ab->", dipole_trace, qbarq_trace)
    return total


def target_Kqbarq_subtraction_appendix(U_fund, Kqbarq, u, v):
    """Appendix A contribution from the -J_L S_A J_R subtraction."""

    U_fund = np.asarray(U_fund)
    total = 0.0j
    nsite = U_fund.shape[0]
    s_uv = dipole(U_fund, u, v)
    for z in range(nsite):
        subtraction_bracket = NC * dipole(U_fund, u, z) * dipole(U_fund, z, v) - s_uv / NC
        for zp in range(nsite):
            total += -0.5 * _kqbarq_endpoint_combo(Kqbarq, u, v, z, zp) * subtraction_bracket
    return total


def _missing(sector):
    raise NotImplementedError(
        f"exact KLM Appendix A target for {sector} is unavailable locally; "
        "sector is internal-consistency-only, Appendix A target missing"
    )


def target_KJSSJ_appendix(U_fund, KJSSJ, u, v):
    """Appendix A full KJSSJ target for H_KJSSJ s(u,v).

    Source: the isolated KJSSJ part of WORKNLO.tex lines 1353--1355,
    excluding the tilde-K contribution.
    """

    U_fund = np.asarray(U_fund)
    total = 0.0j
    nsite = U_fund.shape[0]
    s_uv = dipole(U_fund, u, v)
    for z in range(nsite):
        for zp in range(nsite):
            bracket = (
                (NC**3) * dipole(U_fund, u, zp) * dipole(U_fund, zp, z) * dipole(U_fund, z, v)
                - _trace_word(U_fund, ((False, v), (True, z), (False, zp), (True, u), (False, z), (True, zp)))
                - (NC**3) * dipole(U_fund, u, z) * dipole(U_fund, z, v)
                + NC * s_uv
            )
            total += -(KJSSJ[u, v, z, zp] / NC) * bracket
    return total


def target_Kqbarq_appendix(U_fund, Kqbarq, u, v, gens=None):
    """Appendix A full Kqbarq target for H_Kqbarq s(u,v).

    This is the exact trace-current plus subtraction decomposition of the
    Kqbarq Hamiltonian bracket in WORKNLO.tex lines 268--269. The trace-current
    factor uses the source trace-product expression in lines 1174--1177. The
    trace-current subsection alone is not the full sector target.
    """

    return target_Kqbarq_trace_current_appendix(
        U_fund,
        Kqbarq,
        u,
        v,
        gens=gens,
    ) + target_Kqbarq_subtraction_appendix(U_fund, Kqbarq, u, v)


def klm_normalized_cubic_direct_action(raw_direct):
    """Convert current Hermitian-generator cubic direct action to KLM target convention."""

    return -1.0j * raw_direct


def target_KJJSJ_appendix_real(U_fund, KJJSJ, u, v):
    """Appendix A real KJJSJ target for H_KJJSJ s(u,v).

    Source: WORKNLO.tex lines 1208--1224, compact formula lines 1216--1223.
    Assumes the KLM antisymmetry K(w;x,y;z) = -K(w;y,x;z) used on line 1253.
    """

    U_fund = np.asarray(U_fund)
    total = 0.0j
    s_uv = dipole(U_fund, u, v)
    for z in range(U_fund.shape[0]):
        kernel_combo = KJJSJ[v, u, v, z] - KJJSJ[v, v, u, z] - KJJSJ[u, u, v, z] + KJJSJ[u, v, u, z]
        total += 0.5j * kernel_combo * (s_uv - (NC**2) * dipole(U_fund, u, z) * dipole(U_fund, z, v))
    return total


def target_KJJSJ_appendix_virtual(U_fund, KJJSJ, u, v):
    """Appendix A virtual KJJSJ target for H_KJJSJ s(u,v).

    Source: WORKNLO.tex lines 1244--1253, formula body lines 1246--1251.
    """

    U_fund = np.asarray(U_fund)
    total = 0.0j
    s_uv = dipole(U_fund, u, v)
    for z in range(U_fund.shape[0]):
        total += (KJJSJ[u, v, u, z] + KJJSJ[v, u, v, z]) * s_uv
    return 1.0j * (NC**2 - 1.0) * total / 3.0


def target_KJJSJ_appendix(U_fund, KJJSJ, u, v):
    """Appendix A real plus virtual KJJSJ target for H_KJJSJ s(u,v)."""

    return target_KJJSJ_appendix_real(U_fund, KJJSJ, u, v) + target_KJJSJ_appendix_virtual(
        U_fund,
        KJJSJ,
        u,
        v,
    )


def tilde_K_JJSSJ(KJJSSJ, x, y, z, zp):
    """Return the KLM tilde-K combination for KJJSSJ.

    Source: WORKNLO.tex lines 307--311.
    """

    return 0.5j * (
        KJJSSJ[x, x, y, z, zp]
        - KJJSSJ[y, x, y, z, zp]
        - KJJSSJ[x, y, x, z, zp]
        + KJJSSJ[y, y, x, z, zp]
    )


def _kjjssj_eight_kernel_combo(KJJSSJ, u, v, z, zp):
    return (
        KJJSSJ[u, u, u, z, zp]
        - KJJSSJ[u, v, u, z, zp]
        + KJJSSJ[u, v, v, z, zp]
        - KJJSSJ[u, u, v, z, zp]
        + KJJSSJ[v, u, v, z, zp]
        - KJJSSJ[v, u, u, z, zp]
        + KJJSSJ[v, v, u, z, zp]
        - KJJSSJ[v, v, v, z, zp]
    )


def target_KJJSSJ_appendix_real(U_fund, KJJSSJ, u, v):
    """Appendix A real KJJSSJ target for H_KJJSSJ s(u,v).

    Source: WORKNLO.tex lines 1304--1328, final formula lines 1322--1326.
    """

    U_fund = np.asarray(U_fund)
    total = 0.0j
    nsite = U_fund.shape[0]
    for z in range(nsite):
        for zp in range(nsite):
            triple = (NC**3) * dipole(U_fund, z, v) * dipole(U_fund, zp, z) * dipole(U_fund, u, zp)
            trace = _trace_word(U_fund, ((False, v), (True, z), (False, zp), (True, u), (False, z), (True, zp)))
            tilde = tilde_K_JJSSJ(KJJSSJ, u, v, z, zp)
            total += 0.5j * _kjjssj_eight_kernel_combo(KJJSSJ, u, v, z, zp) * triple / NC
            total += tilde * (triple - trace) / NC
    return total


def target_KJJSSJ_appendix_virtual(U_fund, KJJSSJ, u, v):
    """Appendix A virtual KJJSSJ target for H_KJJSSJ s(u,v).

    Source: WORKNLO.tex lines 1334--1341.
    """

    U_fund = np.asarray(U_fund)
    total = 0.0j
    s_uv = dipole(U_fund, u, v)
    for z in range(U_fund.shape[0]):
        for zp in range(U_fund.shape[0]):
            total += tilde_K_JJSSJ(KJJSSJ, u, v, z, zp) * s_uv
    return -((NC**2 - 1.0) / 3.0) * total


def target_KJJSSJ_appendix(U_fund, KJJSSJ, u, v):
    """Appendix A real plus virtual KJJSSJ target for H_KJJSSJ s(u,v)."""

    return target_KJJSSJ_appendix_real(U_fund, KJJSSJ, u, v) + target_KJJSSJ_appendix_virtual(
        U_fund,
        KJJSSJ,
        u,
        v,
    )


def appendix_target_available(sector_name):
    """Return whether an exact Appendix A target is locally implemented."""

    return bool(_AVAILABLE.get(sector_name, False))
