"""Dense coefficient builders for two-generator NLO current checks."""

from __future__ import annotations

import numpy as np


def kjsj_chi_from_kernel(
    S_adj: np.ndarray,
    KJSJ: np.ndarray,
    use_barred: bool = True,
) -> np.ndarray:
    """Build the corrected dense K_JSJ chi tensor.

    The returned shape is (n_site, 8, n_site, 8), representing
    chi^{(x,b)(y,c)} = -sum_z K_JSJ(x,y;z)
        (S_x^{bd} - S_z^{bd})(S_y^{cd} - S_z^{cd}).

    The use_barred flag records intent only; callers provide the actual
    barred or unbarred kernel tensor.
    """

    _ = use_barred
    diff = S_adj[:, None, :, :] - S_adj[None, :, :, :]
    return -np.einsum("xyz,xzbd,yzcd->xbyc", KJSJ, diff, diff, optimize=True)


def kjssj_A_from_kernel(S_adj: np.ndarray, KJSSJ: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Build A_JSSJ^{ab}(x,y) from a dense small-lattice kernel.

    A^{ab}(x,y) = sum_{z,z'} K_JSSJ(x,y;z,z')
        f^{adc} f^{bef} S_z^{de} (S_z'^{cf} - S_z^{cf}).
    """

    diff = S_adj[None, :, :, :] - S_adj[:, None, :, :]
    return np.einsum(
        "xyuv,adc,bef,ude,uvcf->xyab",
        KJSSJ,
        f,
        f,
        S_adj,
        diff,
        optimize=True,
    )


def kjssj_C_left_from_A(A: np.ndarray, S_adj: np.ndarray) -> np.ndarray:
    """Convert A^{ab}(x,y) J_L^a(x) J_R^b(y) to left-basis C.

    The returned shape is (n_site, 8, n_site, 8), representing
    C^{(x,a)(y,h)} = A^{ab}(x,y) S_y^{hb}.
    """

    return np.einsum("xyab,yhb->xayh", A, S_adj, optimize=True)


def qbarq_trace_block(Uz: np.ndarray, Uzp: np.ndarray, gens: np.ndarray) -> np.ndarray:
    """Return B^{ab}(z,z') = 2 Tr(Uz^dagger t^a Uzp t^b)."""

    return np.real_if_close(
        np.einsum("ij,ajk,kl,bli->ab", Uz.conj().T, gens, Uzp, gens, optimize=True)
        * 2.0
    )


def kqbarq_A_from_kernel(
    U_fund: np.ndarray,
    S_adj: np.ndarray,
    Kqbarq: np.ndarray,
    gens: np.ndarray,
) -> np.ndarray:
    """Build A_qbarq^{ab}(x,y) from a dense small-lattice kernel.

    A^{ab}(x,y) = sum_{z,z'} K_qbarq(x,y;z,z')
        [2 Tr(U_z^dagger t^a U_z' t^b) - S_adj[z,a,b]].

    The result is real for the z<->z' symmetric real kernels used in the
    small-lattice checks. For nonsymmetric synthetic kernels, the trace block
    can retain a physical complex diagnostic component.
    """

    U_fund = np.asarray(U_fund)
    trace_blocks = np.empty(
        (U_fund.shape[0], U_fund.shape[0], gens.shape[0], gens.shape[0]),
        dtype=complex,
    )
    for z in range(U_fund.shape[0]):
        for zp in range(U_fund.shape[0]):
            trace_blocks[z, zp] = qbarq_trace_block(U_fund[z], U_fund[zp], gens)

    block = trace_blocks - S_adj[:, None, :, :]
    return np.real_if_close(np.einsum("xyuv,uvab->xyab", Kqbarq, block, optimize=True))


def kqbarq_C_left_from_A(A: np.ndarray, S_adj: np.ndarray) -> np.ndarray:
    """Alias for converting an ordered A^{ab} J_L^a J_R^b block to left basis."""

    return kjssj_C_left_from_A(A, S_adj)


def sym_asym_parts(C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return symmetric and antisymmetric parts under combined index exchange."""

    C_T = np.transpose(C, (2, 3, 0, 1))
    return 0.5 * (C + C_T), 0.5 * (C - C_T)
