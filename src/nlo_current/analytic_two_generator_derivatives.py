"""Analytic local dK2 contractions for validated two-generator sectors."""

from __future__ import annotations

import numpy as np

from .analytic_lie_derivatives import left_derivative_adjoint
from .nlo_current_skeleton import flatten_index, unflatten_index
from .two_generator_terms import kjssj_A_from_kernel, kqbarq_A_from_kernel


def _flatten_K2_tensor(C: np.ndarray) -> np.ndarray:
    nsite, n_color = C.shape[0], C.shape[1]
    out = np.zeros((nsite * n_color, nsite * n_color), dtype=np.result_type(C, complex))
    for x in range(nsite):
        for a in range(n_color):
            row = flatten_index(x, a, n_color)
            for y in range(nsite):
                for b in range(n_color):
                    out[row, flatten_index(y, b, n_color)] = C[x, a, y, b]
    return np.real_if_close(out)


def _contract_dK2_from_local_derivative(S_adj, f, derivative_tensor_builder) -> np.ndarray:
    nsite, n_color = S_adj.shape[:2]
    dim = nsite * n_color
    dtype = np.result_type(S_adj, f, complex)
    dK2 = np.zeros(dim, dtype=dtype)
    for b_index in range(dim):
        site, color = unflatten_index(b_index, n_color)
        dC = derivative_tensor_builder(site, color)
        dK2 += _flatten_K2_tensor(dC)[:, b_index]
    return np.real_if_close(dK2)


def _dS_array(S_adj, f, site: int, color: int) -> np.ndarray:
    dS = np.zeros_like(S_adj, dtype=np.result_type(S_adj, f, float))
    dS[site] = left_derivative_adjoint(S_adj, f, site, color, site)
    return dS


def analytic_dK2_KJSJ(S_adj: np.ndarray, KJSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    r"""Return analytic \(L_B K2^{AB}\) for the K_JSJ sector."""

    S = np.asarray(S_adj)
    diff = S[:, None, :, :] - S[None, :, :, :]

    def dC(site, color):
        dS = _dS_array(S, f, site, color)
        ddiff = dS[:, None, :, :] - dS[None, :, :, :]
        dchi = -np.einsum("xyz,xzbd,yzcd->xbyc", KJSJ, ddiff, diff, optimize=True)
        dchi += -np.einsum("xyz,xzbd,yzcd->xbyc", KJSJ, diff, ddiff, optimize=True)
        return 2.0 * dchi

    return _contract_dK2_from_local_derivative(S, f, dC)


def _analytic_dK2_KJSSJ_like(S_adj: np.ndarray, KJSSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    S = np.asarray(S_adj)
    diff = S[None, :, :, :] - S[:, None, :, :]
    A = kjssj_A_from_kernel(S, KJSSJ, f)

    def dC(site, color):
        dS = _dS_array(S, f, site, color)
        ddiff = dS[None, :, :, :] - dS[:, None, :, :]
        dA = np.einsum("xyuv,adc,bef,ude,uvcf->xyab", KJSSJ, f, f, dS, diff, optimize=True)
        dA += np.einsum("xyuv,adc,bef,ude,uvcf->xyab", KJSSJ, f, f, S, ddiff, optimize=True)
        return (
            np.einsum("xyab,yhb->xayh", dA, S, optimize=True)
            + np.einsum("xyab,yhb->xayh", A, dS, optimize=True)
        )

    return _contract_dK2_from_local_derivative(S, f, dC)


def analytic_dK2_KJSSJ(S_adj: np.ndarray, KJSSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    r"""Return analytic \(L_B K2^{AB}\) for the K_JSSJ sector."""

    return _analytic_dK2_KJSSJ_like(S_adj, KJSSJ, f)


def _qbarq_trace_blocks(U_fund, gens) -> np.ndarray:
    U = np.asarray(U_fund)
    nsite = U.shape[0]
    n_color = gens.shape[0]
    blocks = np.zeros((nsite, nsite, n_color, n_color), dtype=complex)
    for z in range(nsite):
        Udag = U[z].conj().T
        for zp in range(nsite):
            for a in range(n_color):
                for b in range(n_color):
                    blocks[z, zp, a, b] = 2.0 * np.trace(Udag @ gens[a] @ U[zp] @ gens[b])
    return blocks


def _qbarq_dtrace_blocks(U_fund, gens, site: int, color: int) -> np.ndarray:
    U = np.asarray(U_fund)
    nsite = U.shape[0]
    n_color = gens.shape[0]
    out = np.zeros((nsite, nsite, n_color, n_color), dtype=complex)
    for z in range(nsite):
        Udag = U[z].conj().T
        dUdag = -1.0j * Udag @ gens[color] if z == site else np.zeros((3, 3), dtype=complex)
        for zp in range(nsite):
            dUzp = 1.0j * gens[color] @ U[zp] if zp == site else np.zeros((3, 3), dtype=complex)
            for a in range(n_color):
                for b in range(n_color):
                    out[z, zp, a, b] = 2.0 * np.trace(
                        dUdag @ gens[a] @ U[zp] @ gens[b]
                        + Udag @ gens[a] @ dUzp @ gens[b]
                    )
    return out


def _analytic_dK2_Kqbarq_component(
    U_fund: np.ndarray,
    S_adj: np.ndarray,
    Kqbarq: np.ndarray,
    gens: np.ndarray,
    f: np.ndarray,
    *,
    component: str,
) -> np.ndarray:
    S = np.asarray(S_adj)
    trace_blocks = _qbarq_trace_blocks(U_fund, gens)
    if component == "trace":
        block = trace_blocks
    elif component == "subtraction":
        block = -S[:, None, :, :]
    elif component == "full":
        block = trace_blocks - S[:, None, :, :]
    else:
        raise ValueError("component must be 'trace', 'subtraction', or 'full'")

    A = np.einsum("xyuv,uvab->xyab", Kqbarq, block, optimize=True)

    def dC(site, color):
        dS = _dS_array(S, f, site, color)
        if component == "trace":
            dblock = _qbarq_dtrace_blocks(U_fund, gens, site, color)
        elif component == "subtraction":
            dblock = -dS[:, None, :, :]
        else:
            dblock = _qbarq_dtrace_blocks(U_fund, gens, site, color) - dS[:, None, :, :]
        dA = np.einsum("xyuv,uvab->xyab", Kqbarq, dblock, optimize=True)
        return (
            np.einsum("xyab,yhb->xayh", dA, S, optimize=True)
            + np.einsum("xyab,yhb->xayh", A, dS, optimize=True)
        )

    return _contract_dK2_from_local_derivative(S, f, dC)


def analytic_dK2_Kqbarq_trace(
    U_fund: np.ndarray,
    S_adj: np.ndarray,
    Kqbarq: np.ndarray,
    gens: np.ndarray,
    f: np.ndarray,
) -> np.ndarray:
    r"""Return analytic trace-product contribution to \(L_BK2^{AB}\)."""

    return _analytic_dK2_Kqbarq_component(U_fund, S_adj, Kqbarq, gens, f, component="trace")


def analytic_dK2_Kqbarq_subtraction(
    U_fund: np.ndarray,
    S_adj: np.ndarray,
    Kqbarq: np.ndarray,
    gens: np.ndarray,
    f: np.ndarray,
) -> np.ndarray:
    r"""Return analytic subtraction contribution to \(L_BK2^{AB}\)."""

    return _analytic_dK2_Kqbarq_component(
        U_fund,
        S_adj,
        Kqbarq,
        gens,
        f,
        component="subtraction",
    )


def analytic_dK2_Kqbarq(
    U_fund: np.ndarray,
    S_adj: np.ndarray,
    Kqbarq: np.ndarray,
    gens: np.ndarray,
    f: np.ndarray,
) -> np.ndarray:
    r"""Return full analytic \(L_BK2^{AB}\) for K_qbarq."""

    return _analytic_dK2_Kqbarq_component(U_fund, S_adj, Kqbarq, gens, f, component="full")


def kqbarq_trace_A_from_kernel(U_fund, Kqbarq, gens):
    """Build the K_qbarq trace-product A block for oracle tests."""

    return np.einsum("xyuv,uvab->xyab", Kqbarq, _qbarq_trace_blocks(U_fund, gens), optimize=True)


def kqbarq_subtraction_A_from_kernel(S_adj, Kqbarq):
    """Build the K_qbarq subtraction A block for oracle tests."""

    return -np.einsum("xyuv,uab->xyab", Kqbarq, S_adj, optimize=True)
