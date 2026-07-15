"""Dense K_JJSJ cubic coefficient builders for small-lattice checks."""

from __future__ import annotations

import numpy as np


def kjjsj_A_LLR_from_kernel(S_adj: np.ndarray, KJJSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build A_LLR[x,y,w,d,e,a] = sum_z K[w,x,y,z] f[b,d,e] S[z,b,a]."""

    return np.einsum("wxyz,bde,zba->xywdea", KJJSJ, f, S_adj, optimize=True)


def kjjsj_B_LRR_from_kernel(S_adj: np.ndarray, KJJSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build B_LRR[w,x,y,a,d,e] = -sum_z K[w,x,y,z] f[b,d,e] S[z,a,b]."""

    return -np.einsum("wxyz,bde,zab->wxyade", KJJSJ, f, S_adj, optimize=True)


def kjjsj_V_virtual_from_kernel(KJJSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build V[x,y,w,d,e,b] = (1/3) sum_z K[w,x,y,z] f[b,d,e]."""

    return (1.0 / 3.0) * np.einsum("wxyz,bde->xywdeb", KJJSJ, f, optimize=True)


def flatten_cubic_index(site: int, color: int, n_color: int = 8) -> int:
    """Return combined index for a site/color pair."""

    return site * n_color + color


def synthetic_kjjsj_kernel(
    nsite: int,
    rng: np.random.Generator,
    xy_symmetry: str = "antisymmetric",
) -> np.ndarray:
    """Generate a dense synthetic KJJSJ[w,x,y,z] kernel.

    xy_symmetry can be "antisymmetric", "symmetric", or "unconstrained".
    Both symmetry choices are tested because this workflow does not assume the
    production K_JJSJ coordinate symmetry.
    """

    raw = rng.normal(size=(nsite, nsite, nsite, nsite))
    if xy_symmetry == "antisymmetric":
        kernel = 0.5 * (raw - np.swapaxes(raw, 1, 2))
    elif xy_symmetry == "symmetric":
        kernel = 0.5 * (raw + np.swapaxes(raw, 1, 2))
    elif xy_symmetry == "unconstrained":
        kernel = raw
    else:
        raise ValueError(f"unknown xy_symmetry: {xy_symmetry}")

    norm = np.linalg.norm(kernel)
    if norm == 0.0:
        return kernel
    return kernel / norm


def kjjssj_A_LLR_from_kernel(S_adj: np.ndarray, KJJSSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build A_LLR[x,y,w,d,e,a] for the K_JJSSJ LLR block."""

    return np.einsum("wxyuv,acb,udc,veb->xywdea", KJJSSJ, f, S_adj, S_adj, optimize=True)


def kjjssj_B_LRR_from_kernel(S_adj: np.ndarray, KJJSSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build B_LRR[w,x,y,a,d,e] for the K_JJSSJ LRR block."""

    return -np.einsum("wxyuv,acb,ucd,vbe->wxyade", KJJSSJ, f, S_adj, S_adj, optimize=True)


def kjjssj_V_virtual_from_kernel(KJJSSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build V[x,y,w,c,b,a] = (1/3) sum_{z,z'} K[w,x,y,z,z'] f[a,c,b]."""

    return (1.0 / 3.0) * np.einsum("wxyuv,acb->xywcba", KJJSSJ, f, optimize=True)


def synthetic_kjjssj_kernel(
    nsite: int,
    rng: np.random.Generator,
    klm_antisym: bool | str = True,
) -> np.ndarray:
    """Generate a dense synthetic KJJSSJ[w,x,y,z,zp] kernel.

    If klm_antisym is True, impose K(w;x,y;z,z') = -K(w;y,x;z',z).
    If klm_antisym is False or "unconstrained", return an unconstrained random
    kernel. The string "symmetric" imposes the corresponding plus relation and
    is provided only for diagnostic stress tests.
    """

    raw = rng.normal(size=(nsite, nsite, nsite, nsite, nsite))
    transformed = np.transpose(raw, (0, 2, 1, 4, 3))

    if klm_antisym is True or klm_antisym == "klm_antisym":
        kernel = 0.5 * (raw - transformed)
    elif klm_antisym is False or klm_antisym == "unconstrained":
        kernel = raw
    elif klm_antisym == "symmetric":
        kernel = 0.5 * (raw + transformed)
    else:
        raise ValueError(f"unknown klm_antisym option: {klm_antisym}")

    norm = np.linalg.norm(kernel)
    if norm == 0.0:
        return kernel
    return kernel / norm


def dense_cubic_tensor_from_blocks(
    A_LLR: np.ndarray | None = None,
    V_LLL: np.ndarray | None = None,
) -> np.ndarray:
    """Convert selected ordered left-derivative blocks to a dense tensor.

    This optional helper is intentionally minimal. It fills a tensor C[A,B,C]
    indexed by combined site/color, using the ordered derivative index order
    represented by each block. It does not symmetrize the cubic tensor.
    """

    source = A_LLR if A_LLR is not None else V_LLL
    if source is None:
        raise ValueError("provide at least one cubic block")

    nsite = source.shape[0]
    n_color = 8
    tensor = np.zeros((nsite * n_color, nsite * n_color, nsite * n_color))

    if A_LLR is not None:
        for x in range(nsite):
            for y in range(nsite):
                for w in range(nsite):
                    for d in range(n_color):
                        for e in range(n_color):
                            for a in range(n_color):
                                i = flatten_cubic_index(x, d, n_color)
                                j = flatten_cubic_index(y, e, n_color)
                                k = flatten_cubic_index(w, a, n_color)
                                tensor[i, j, k] += A_LLR[x, y, w, d, e, a]

    if V_LLL is not None:
        for x in range(nsite):
            for y in range(nsite):
                for w in range(nsite):
                    for d in range(n_color):
                        for e in range(n_color):
                            for b in range(n_color):
                                i = flatten_cubic_index(x, d, n_color)
                                j = flatten_cubic_index(y, e, n_color)
                                k = flatten_cubic_index(w, b, n_color)
                                tensor[i, j, k] += V_LLL[x, y, w, d, e, b]

    return tensor
