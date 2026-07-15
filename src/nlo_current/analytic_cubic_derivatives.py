"""Analytic cubic coefficient derivatives for validated tiny diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analytic_lie_derivatives import left_derivative_adjoint, second_left_derivative_adjoint
from .cubic_commutator_terms import (
    canonicalize_cubic_block_terms,
    cubic_block_terms_from_LLR,
    cubic_block_terms_from_LRR,
    cubic_block_terms_from_virtual_LLL,
    cubic_block_terms_from_virtual_RRR,
    llr_right_to_left,
    lrr_right_to_left,
    virtual_rrr_right_to_left,
)
from .nlo_current_skeleton import flatten_index, unflatten_index
from .three_generator_terms import (
    kjjsj_A_LLR_from_kernel,
    kjjsj_B_LRR_from_kernel,
    kjjsj_V_virtual_from_kernel,
)


KJJSJ_BLOCKS = ("LLR", "LRR", "virtual_LLL", "virtual_RRR")


@dataclass(frozen=True)
class KJJSJSectorData:
    """Inputs needed for KJJSJ analytic derivative contractions."""

    S_adj: np.ndarray
    KJJSJ: np.ndarray
    f: np.ndarray


def _pending(name: str):
    raise NotImplementedError(
        f"{name} is pending: use backend='hybrid_local_fd' or 'finite_difference' explicitly"
    )


def _sector_data(sector_data) -> KJJSJSectorData:
    if isinstance(sector_data, KJJSJSectorData):
        return sector_data
    if not isinstance(sector_data, dict):
        raise ValueError("KJJSJ analytic derivatives require sector_data with S_adj, KJJSJ, and f")
    return KJJSJSectorData(
        S_adj=np.asarray(sector_data["S_adj"]),
        KJJSJ=np.asarray(sector_data["KJJSJ"]),
        f=np.asarray(sector_data["f"]),
    )


def _dim(data: KJJSJSectorData) -> int:
    return int(data.S_adj.shape[0]) * int(data.S_adj.shape[1])


def _zero_S_like(data: KJJSJSectorData):
    return np.zeros_like(data.S_adj, dtype=np.result_type(data.S_adj, data.f, float))


def _dS(data: KJJSJSectorData, flat_index: int) -> np.ndarray:
    n_color = int(data.S_adj.shape[1])
    site, color = unflatten_index(flat_index, n_color)
    out = _zero_S_like(data)
    out[site] = left_derivative_adjoint(data.S_adj, data.f, site, color, site)
    return out


def _ddS(data: KJJSJSectorData, first_index: int, second_index: int) -> np.ndarray:
    n_color = int(data.S_adj.shape[1])
    first_site, first_color = unflatten_index(first_index, n_color)
    second_site, second_color = unflatten_index(second_index, n_color)
    out = _zero_S_like(data)
    if first_site == second_site:
        out[first_site] = second_left_derivative_adjoint(
            data.S_adj,
            data.f,
            first_site,
            first_color,
            second_site,
            second_color,
            first_site,
        )
    return out


def _A(data: KJJSJSectorData, S=None):
    return kjjsj_A_LLR_from_kernel(data.S_adj if S is None else S, data.KJJSJ, data.f)


def _B(data: KJJSJSectorData, S=None):
    return kjjsj_B_LRR_from_kernel(data.S_adj if S is None else S, data.KJJSJ, data.f)


def _V(data: KJJSJSectorData):
    return kjjsj_V_virtual_from_kernel(data.KJJSJ, data.f)


def _llr(A, S):
    return np.einsum("xywdea,wha->xywdeh", A, S, optimize=True)


def _lrr(B, Sx, Sy):
    return np.einsum("wxyade,xpd,yqe->wxyapq", B, Sx, Sy, optimize=True)


def _vrrr(V, Sx, Sy, Sw):
    return -np.einsum("xywcba,xpc,yqb,wha->xywpqh", V, Sx, Sy, Sw, optimize=True)


def _base_block_arrays(data: KJJSJSectorData) -> dict[str, np.ndarray]:
    S = data.S_adj
    A = _A(data)
    B = _B(data)
    V = _V(data)
    return {
        "LLR": llr_right_to_left(A, S),
        "LRR": lrr_right_to_left(B, S),
        "virtual_LLL": V,
        "virtual_RRR": virtual_rrr_right_to_left(V, S, sign=-1.0),
    }


def _first_block_arrays(data: KJJSJSectorData, flat_index: int) -> dict[str, np.ndarray]:
    S = data.S_adj
    dS = _dS(data, flat_index)
    A = _A(data)
    B = _B(data)
    V = _V(data)
    dA = _A(data, dS)
    dB = _B(data, dS)
    zero_V = np.zeros_like(V, dtype=np.result_type(V, dS, complex))
    return {
        "LLR": _llr(dA, S) + _llr(A, dS),
        "LRR": _lrr(dB, S, S) + _lrr(B, dS, S) + _lrr(B, S, dS),
        "virtual_LLL": zero_V,
        "virtual_RRR": _vrrr(V, dS, S, S) + _vrrr(V, S, dS, S) + _vrrr(V, S, S, dS),
    }


def _second_block_arrays(
    data: KJJSJSectorData,
    first_index: int,
    second_index: int,
) -> dict[str, np.ndarray]:
    S = data.S_adj
    dS_first = _dS(data, first_index)
    dS_second = _dS(data, second_index)
    ddS = _ddS(data, first_index, second_index)
    A = _A(data)
    B = _B(data)
    V = _V(data)
    dA_first = _A(data, dS_first)
    dA_second = _A(data, dS_second)
    ddA = _A(data, ddS)
    dB_first = _B(data, dS_first)
    dB_second = _B(data, dS_second)
    ddB = _B(data, ddS)
    zero_V = np.zeros_like(V, dtype=np.result_type(V, ddS, complex))
    llr = (
        _llr(ddA, S)
        + _llr(dA_second, dS_first)
        + _llr(dA_first, dS_second)
        + _llr(A, ddS)
    )
    lrr = (
        _lrr(ddB, S, S)
        + _lrr(B, ddS, S)
        + _lrr(B, S, ddS)
        + _lrr(dB_second, dS_first, S)
        + _lrr(dB_second, S, dS_first)
        + _lrr(dB_first, dS_second, S)
        + _lrr(B, dS_second, dS_first)
        + _lrr(dB_first, S, dS_second)
        + _lrr(B, dS_first, dS_second)
    )
    vrrr = (
        _vrrr(V, ddS, S, S)
        + _vrrr(V, S, ddS, S)
        + _vrrr(V, S, S, ddS)
        + _vrrr(V, dS_second, dS_first, S)
        + _vrrr(V, dS_second, S, dS_first)
        + _vrrr(V, dS_first, dS_second, S)
        + _vrrr(V, S, dS_second, dS_first)
        + _vrrr(V, dS_first, S, dS_second)
        + _vrrr(V, S, dS_first, dS_second)
    )
    return {"LLR": llr, "LRR": lrr, "virtual_LLL": zero_V, "virtual_RRR": vrrr}


def _raw_terms_from_block(block: str, array: np.ndarray) -> dict:
    if block == "LLR":
        return cubic_block_terms_from_LLR(array)
    if block == "LRR":
        return cubic_block_terms_from_LRR(array)
    if block == "virtual_LLL":
        return cubic_block_terms_from_virtual_LLL(array)
    if block == "virtual_RRR":
        return cubic_block_terms_from_virtual_RRR(array)
    raise ValueError(f"unknown KJJSJ block: {block}")


def _arrays_from_terms(cubic_terms, quadratic_terms, linear_terms, dim, n_color=8):
    dtype = np.dtype(float)
    for terms in (cubic_terms, quadratic_terms, linear_terms):
        for coeff in terms.values():
            dtype = np.result_type(dtype, coeff)
    K1 = np.zeros(dim, dtype=dtype)
    K2 = np.zeros((dim, dim), dtype=dtype)
    K3 = np.zeros((dim, dim, dim), dtype=dtype)
    for word, coeff in linear_terms.items():
        K1[flatten_index(*word[0], n_color)] += coeff
    for word, coeff in quadratic_terms.items():
        K2[flatten_index(*word[0], n_color), flatten_index(*word[1], n_color)] += coeff
    for word, coeff in cubic_terms.items():
        K3[
            flatten_index(*word[0], n_color),
            flatten_index(*word[1], n_color),
            flatten_index(*word[2], n_color),
        ] += coeff
    return K1, K2, K3


def _canonical_arrays_for_blocks(data: KJJSJSectorData, block_arrays: dict[str, np.ndarray], blocks=None):
    selected = KJJSJ_BLOCKS if blocks is None else tuple(blocks)
    raw_terms = {}
    for block in selected:
        for word, coeff in _raw_terms_from_block(block, block_arrays[block]).items():
            raw_terms[word] = raw_terms.get(word, 0.0) + coeff
    cubic, quadratic, linear = canonicalize_cubic_block_terms(raw_terms, data.f)
    return _arrays_from_terms(cubic, quadratic, linear, _dim(data), data.S_adj.shape[1])


def kjjsj_terms_from_blocks(sector_data, blocks=None, derivative=None, first_index=None, second_index=None):
    """Return K1/K2/K3 arrays from selected KJJSJ blocks.

    ``derivative`` can be ``None``, ``"first"``, or ``"second"``.
    This helper is analytic and uses no finite differences.
    """

    data = _sector_data(sector_data)
    if derivative is None:
        block_arrays = _base_block_arrays(data)
    elif derivative == "first":
        block_arrays = _first_block_arrays(data, int(first_index))
    elif derivative == "second":
        block_arrays = _second_block_arrays(data, int(first_index), int(second_index))
    else:
        raise ValueError("derivative must be None, 'first', or 'second'")
    return _canonical_arrays_for_blocks(data, block_arrays, blocks=blocks)


def _contraction_result(total, by_block):
    return {"total": np.real_if_close(total), "by_block": {k: np.real_if_close(v) for k, v in by_block.items()}}


def _first_derivative_contractions(data: KJJSJSectorData, *, with_blocks: bool = False):
    dim = _dim(data)
    dtype = np.result_type(data.S_adj, data.KJJSJ, data.f, complex)
    LC = np.zeros((dim, dim), dtype=dtype)
    LB = np.zeros((dim, dim), dtype=dtype)
    LC_blocks = {block: np.zeros((dim, dim), dtype=dtype) for block in KJJSJ_BLOCKS} if with_blocks else {}
    LB_blocks = {block: np.zeros((dim, dim), dtype=dtype) for block in KJJSJ_BLOCKS} if with_blocks else {}
    dK2_comm = np.zeros(dim, dtype=dtype)
    dK2_blocks = {block: np.zeros(dim, dtype=dtype) for block in KJJSJ_BLOCKS} if with_blocks else {}

    for index in range(dim):
        block_arrays = _first_block_arrays(data, index)
        K1_d, K2_d, K3_d = _canonical_arrays_for_blocks(data, block_arrays)
        LC += K3_d[:, :, index]
        LB += K3_d[:, index, :]
        dK2_comm += K2_d[:, index]
        if with_blocks:
            for block in KJJSJ_BLOCKS:
                _, K2_b, K3_b = _canonical_arrays_for_blocks(data, block_arrays, blocks=(block,))
                LC_blocks[block] += K3_b[:, :, index]
                LB_blocks[block] += K3_b[:, index, :]
                dK2_blocks[block] += K2_b[:, index]

    return {
        "LC": _contraction_result(LC, LC_blocks),
        "LB": _contraction_result(LB, LB_blocks),
        "dK2_comm": _contraction_result(dK2_comm, dK2_blocks),
    }


def analytic_LC_K3_KJJSJ(U=None, physical_terms=None, *, sector_data=None, return_diagnostics=False):
    """Return ``(LC_K3)^{AB} = L_C K3^{ABC}`` for KJJSJ only."""

    data = _sector_data(sector_data)
    result = _first_derivative_contractions(data, with_blocks=return_diagnostics)["LC"]
    return result if return_diagnostics else result["total"]


def analytic_LB_K3_KJJSJ(U=None, physical_terms=None, *, sector_data=None, return_diagnostics=False):
    """Return ``(LB_K3)^{AC} = L_B K3^{ABC}`` for KJJSJ only."""

    data = _sector_data(sector_data)
    result = _first_derivative_contractions(data, with_blocks=return_diagnostics)["LB"]
    return result if return_diagnostics else result["total"]


def analytic_dK2_comm_KJJSJ(U=None, physical_terms=None, *, sector_data=None, return_diagnostics=False):
    """Return KJJSJ quadratic commutator contribution to ``L_B K2^{AB}``."""

    data = _sector_data(sector_data)
    result = _first_derivative_contractions(data, with_blocks=return_diagnostics)["dK2_comm"]
    return result if return_diagnostics else result["total"]


def analytic_d2K3_KJJSJ(U=None, physical_terms=None, *, sector_data=None, return_diagnostics=False):
    """Return ordered ``d2K3^A = L_B L_C K3^{ABC}`` for KJJSJ only."""

    data = _sector_data(sector_data)
    dim = _dim(data)
    dtype = np.result_type(data.S_adj, data.KJJSJ, data.f, complex)
    d2 = np.zeros(dim, dtype=dtype)
    d2_blocks = {block: np.zeros(dim, dtype=dtype) for block in KJJSJ_BLOCKS}

    for b_index in range(dim):
        for c_index in range(dim):
            block_arrays = _second_block_arrays(data, b_index, c_index)
            _, _, K3_dd = _canonical_arrays_for_blocks(data, block_arrays)
            d2 += K3_dd[:, b_index, c_index]
            for block in KJJSJ_BLOCKS:
                _, _, K3_block = _canonical_arrays_for_blocks(data, block_arrays, blocks=(block,))
                d2_blocks[block] += K3_block[:, b_index, c_index]

    result = _contraction_result(d2, d2_blocks)
    return result if return_diagnostics else result["total"]


def analytic_dK1_comm_KJJSJ(U=None, physical_terms=None, *, sector_data=None):
    """Classify the current KJJSJ linear commutator correction."""

    data = _sector_data(sector_data)
    K1, _, _ = _canonical_arrays_for_blocks(data, _base_block_arrays(data))
    norm = float(np.linalg.norm(K1))
    status = "structurally_zero" if norm == 0.0 else "nonzero"
    return {"K1_comm": np.real_if_close(K1), "norm": norm, "status": status}


def analytic_first_derivatives_KJJSJ(*args, sector_data=None, return_diagnostics=False, **kwargs):
    """Return both KJJSJ first-derivative contractions."""

    data = _sector_data(sector_data)
    contractions = _first_derivative_contractions(data, with_blocks=return_diagnostics)
    out = {
        "LC_K3": contractions["LC"]["total"],
        "LB_K3": contractions["LB"]["total"],
        "dK2_comm": contractions["dK2_comm"]["total"],
    }
    if return_diagnostics:
        out["by_block"] = {
            "LC_K3": contractions["LC"]["by_block"],
            "LB_K3": contractions["LB"]["by_block"],
            "dK2_comm": contractions["dK2_comm"]["by_block"],
        }
    return out


def analytic_first_derivatives_KJJSSJ(*args, **kwargs):
    """Pending analytic LC/LB derivatives for K_JJSSJ."""

    _pending("analytic_first_derivatives_KJJSSJ")


def analytic_d2K3_KJJSSJ(*args, **kwargs):
    """Pending analytic d2K3 derivative for K_JJSSJ."""

    _pending("analytic_d2K3_KJJSSJ")


def hybrid_local_fd_d2K3_KJJSJ(*args, **kwargs):
    """Explicitly unavailable placeholder for a hybrid local-FD diagnostic."""

    _pending("hybrid_local_fd_d2K3_KJJSJ")
