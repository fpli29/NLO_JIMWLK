"""Dense non-production NLO current skeleton for tiny-lattice diagnostics."""

from __future__ import annotations

from dataclasses import dataclass, field
from copy import deepcopy

import numpy as np

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
from .three_generator_terms import (
    kjjsj_A_LLR_from_kernel,
    kjjsj_B_LRR_from_kernel,
    kjjsj_V_virtual_from_kernel,
    kjjssj_A_LLR_from_kernel,
    kjjssj_B_LRR_from_kernel,
    kjjssj_V_virtual_from_kernel,
)
from .two_generator_terms import (
    kjsj_chi_from_kernel,
    kjssj_A_from_kernel,
    kjssj_C_left_from_A,
    kqbarq_A_from_kernel,
    kqbarq_C_left_from_A,
)


@dataclass
class NLOCurrentTerms:
    """Dense normal-form current terms for non-production diagnostics."""

    K1: np.ndarray
    K2: np.ndarray
    K3: np.ndarray
    metadata: dict = field(default_factory=dict)

    @property
    def dim(self) -> int:
        return int(self.K1.shape[0])

    def validate_shapes(self) -> None:
        if self.K1.ndim != 1:
            raise ValueError("K1 must have shape (D,)")
        dim = self.K1.shape[0]
        if self.K2.shape != (dim, dim):
            raise ValueError(f"K2 must have shape {(dim, dim)}, got {self.K2.shape}")
        if self.K3.shape != (dim, dim, dim):
            raise ValueError(f"K3 must have shape {(dim, dim, dim)}, got {self.K3.shape}")

    def norms(self) -> dict:
        self.validate_shapes()
        return {
            "K1": float(np.linalg.norm(self.K1)),
            "K2": float(np.linalg.norm(self.K2)),
            "K3": float(np.linalg.norm(self.K3)),
        }


def combined_dim(nsite: int, n_color: int = 8) -> int:
    return int(nsite) * int(n_color)


def flatten_index(site: int, color: int, n_color: int = 8) -> int:
    return int(site) * int(n_color) + int(color)


def unflatten_index(index: int, n_color: int = 8) -> tuple[int, int]:
    return int(index) // int(n_color), int(index) % int(n_color)


def empty_terms(nsite: int, n_color: int = 8, dtype=float) -> NLOCurrentTerms:
    dim = combined_dim(nsite, n_color)
    terms = NLOCurrentTerms(
        K1=np.zeros(dim, dtype=dtype),
        K2=np.zeros((dim, dim), dtype=dtype),
        K3=np.zeros((dim, dim, dim), dtype=dtype),
        metadata={"nsite": nsite, "n_color": n_color, "sectors": {}, "warnings": []},
    )
    terms.validate_shapes()
    return terms


def _merge_metadata(lhs: dict, rhs: dict, label: str | None = None) -> dict:
    metadata = deepcopy(lhs)
    metadata.setdefault("sectors", {})
    metadata.setdefault("warnings", [])
    rhs_copy = deepcopy(rhs)
    for key, value in rhs_copy.get("sectors", {}).items():
        metadata["sectors"][key if label is None else f"{label}:{key}"] = value
    metadata["warnings"].extend(rhs_copy.get("warnings", []))
    metadata.setdefault("sources", [])
    metadata["sources"].extend(rhs_copy.get("sources", []))
    if label is not None:
        metadata.setdefault("labels", []).append(label)
    for key, value in rhs_copy.items():
        if key not in {"sectors", "warnings", "sources"}:
            metadata.setdefault("extra", {})[key if label is None else f"{label}:{key}"] = value
    return metadata


def add_terms(lhs: NLOCurrentTerms, rhs: NLOCurrentTerms, label: str | None = None) -> NLOCurrentTerms:
    lhs.validate_shapes()
    rhs.validate_shapes()
    if lhs.dim != rhs.dim:
        raise ValueError("cannot add terms with different dimensions")
    out = NLOCurrentTerms(
        K1=lhs.K1 + rhs.K1,
        K2=lhs.K2 + rhs.K2,
        K3=lhs.K3 + rhs.K3,
        metadata=_merge_metadata(lhs.metadata, rhs.metadata, label),
    )
    out.validate_shapes()
    return out


def symmetrize_K2(K2: np.ndarray) -> np.ndarray:
    """Diagnostic-only K2 symmetrization; never applied by default."""

    return 0.5 * (K2 + K2.T)


def symmetrize_K3(K3: np.ndarray) -> np.ndarray:
    """Diagnostic-only full K3 symmetrization; never applied by default."""

    return (
        K3
        + np.transpose(K3, (0, 2, 1))
        + np.transpose(K3, (1, 0, 2))
        + np.transpose(K3, (1, 2, 0))
        + np.transpose(K3, (2, 0, 1))
        + np.transpose(K3, (2, 1, 0))
    ) / 6.0


def _empty_like(S_adj: np.ndarray, dtype=float) -> NLOCurrentTerms:
    return empty_terms(S_adj.shape[0], S_adj.shape[1], dtype=dtype)


def _fill_K2_from_tensor(K2: np.ndarray, tensor: np.ndarray) -> None:
    nsite, n_color = tensor.shape[0], tensor.shape[1]
    for x in range(nsite):
        for a in range(n_color):
            i = flatten_index(x, a, n_color)
            for y in range(nsite):
                for b in range(n_color):
                    K2[i, flatten_index(y, b, n_color)] += tensor[x, a, y, b]


def _terms_dtype(*term_dicts):
    dtype = np.dtype(float)
    for terms in term_dicts:
        for coeff in terms.values():
            dtype = np.result_type(dtype, coeff)
    return dtype


def _word_terms_to_arrays(cubic_terms, quadratic_terms, linear_terms, dim, n_color=8):
    dtype = _terms_dtype(cubic_terms, quadratic_terms, linear_terms)
    K1 = np.zeros(dim, dtype=dtype)
    K2 = np.zeros((dim, dim), dtype=dtype)
    K3 = np.zeros((dim, dim, dim), dtype=dtype)
    for word, coeff in linear_terms.items():
        K1[flatten_index(*word[0], n_color)] += coeff
    for word, coeff in quadratic_terms.items():
        K2[
            flatten_index(*word[0], n_color),
            flatten_index(*word[1], n_color),
        ] += coeff
    for word, coeff in cubic_terms.items():
        K3[
            flatten_index(*word[0], n_color),
            flatten_index(*word[1], n_color),
            flatten_index(*word[2], n_color),
        ] += coeff
    return K1, K2, K3


def _sector_terms(nsite, n_color, name, K1=None, K2=None, K3=None, extra=None, warnings=None):
    dim = combined_dim(nsite, n_color)
    dtype = np.result_type(
        float,
        *(np.asarray(arr).dtype for arr in (K1, K2, K3) if arr is not None),
    )
    terms = NLOCurrentTerms(
        K1=np.zeros(dim, dtype=dtype) if K1 is None else K1,
        K2=np.zeros((dim, dim), dtype=dtype) if K2 is None else K2,
        K3=np.zeros((dim, dim, dim), dtype=dtype) if K3 is None else K3,
        metadata={"nsite": nsite, "n_color": n_color, "sectors": {}, "warnings": warnings or []},
    )
    terms.validate_shapes()
    info = {"norms": terms.norms(), "nonproduction": True}
    if extra:
        info.update(extra)
    terms.metadata["sectors"][name] = info
    return terms


def assemble_kjsj_terms(U_fund, S_adj, KJSJ, metadata_only=False) -> NLOCurrentTerms:
    """Build the K_JSJ diagnostic contribution to K2."""

    _ = U_fund
    nsite, n_color = S_adj.shape[:2]
    terms = _empty_like(S_adj)
    if not metadata_only:
        chi = kjsj_chi_from_kernel(S_adj, KJSJ)
        _fill_K2_from_tensor(terms.K2, 2.0 * chi)
    terms.metadata["sectors"]["KJSJ"] = {
        "normal_form": "K2",
        "metadata_only": metadata_only,
        "norms": terms.norms(),
    }
    return terms


def assemble_kjssj_terms(U_fund, S_adj, KJSSJ, f, metadata_only=False) -> NLOCurrentTerms:
    """Build the K_JSSJ diagnostic contribution to K2."""

    _ = U_fund
    nsite, n_color = S_adj.shape[:2]
    terms = _empty_like(S_adj)
    if not metadata_only:
        A = kjssj_A_from_kernel(S_adj, KJSSJ, f)
        C = kjssj_C_left_from_A(A, S_adj)
        _fill_K2_from_tensor(terms.K2, C)
    terms.metadata["sectors"]["KJSSJ"] = {
        "normal_form": "K2",
        "ordered_generic_second_order": True,
        "metadata_only": metadata_only,
        "norms": terms.norms(),
    }
    return terms


def assemble_kqbarq_terms(U_fund, S_adj, Kqbarq, gens, metadata_only=False) -> NLOCurrentTerms:
    """Build the K_qbarq diagnostic contribution to K2."""

    if metadata_only:
        terms = _empty_like(S_adj)
    else:
        A = kqbarq_A_from_kernel(U_fund, S_adj, Kqbarq, gens)
        C = kqbarq_C_left_from_A(A, S_adj)
        C = np.real_if_close(C)
        terms = _empty_like(S_adj, dtype=np.result_type(C, float))
        _fill_K2_from_tensor(terms.K2, C)
    terms.metadata["sectors"]["Kqbarq"] = {
        "normal_form": "K2",
        "ordered_generic_second_order": True,
        "metadata_only": metadata_only,
        "norms": terms.norms(),
    }
    return terms


def _assemble_cubic_from_raw_terms(S_adj, f, raw_terms, include_commutators, sector_name):
    nsite, n_color = S_adj.shape[:2]
    dim = combined_dim(nsite, n_color)
    cubic, quadratic, linear = canonicalize_cubic_block_terms(raw_terms, f)
    if not include_commutators:
        quadratic = {}
        linear = {}
    K1, K2, K3 = _word_terms_to_arrays(cubic, quadratic, linear, dim, n_color)
    return _sector_terms(
        nsite,
        n_color,
        sector_name,
        K1=K1,
        K2=K2,
        K3=K3,
        extra={
            "normal_form": "K3+K2_comm+K1_comm",
            "include_commutators": include_commutators,
            "raw_terms": len(raw_terms),
            "canonical_cubic_terms": len(cubic),
            "quadratic_comm_terms": len(quadratic),
            "linear_comm_terms": len(linear),
        },
    )


def assemble_kjjsj_terms(
    U_fund,
    S_adj,
    KJJSJ,
    f,
    include_commutators=True,
    metadata_only=False,
) -> NLOCurrentTerms:
    """Build K_JJSJ diagnostic K3 plus commutator corrections."""

    _ = U_fund
    if metadata_only:
        terms = _empty_like(S_adj)
        terms.metadata["sectors"]["KJJSJ"] = {
            "normal_form": "K3+K2_comm+K1_comm",
            "include_commutators": include_commutators,
            "metadata_only": True,
            "norms": terms.norms(),
        }
        return terms
    A = kjjsj_A_LLR_from_kernel(S_adj, KJJSJ, f)
    B = kjjsj_B_LRR_from_kernel(S_adj, KJJSJ, f)
    V = kjjsj_V_virtual_from_kernel(KJJSJ, f)
    raw_terms = {}
    for source in (
        cubic_block_terms_from_LLR(llr_right_to_left(A, S_adj)),
        cubic_block_terms_from_LRR(lrr_right_to_left(B, S_adj)),
        cubic_block_terms_from_virtual_LLL(V),
        cubic_block_terms_from_virtual_RRR(virtual_rrr_right_to_left(V, S_adj, sign=-1.0)),
    ):
        for word, coeff in source.items():
            raw_terms[word] = raw_terms.get(word, 0.0) + coeff
    return _assemble_cubic_from_raw_terms(S_adj, f, raw_terms, include_commutators, "KJJSJ")


def assemble_kjjssj_terms(
    U_fund,
    S_adj,
    KJJSSJ,
    f,
    include_commutators=True,
    metadata_only=False,
) -> NLOCurrentTerms:
    """Build K_JJSSJ diagnostic K3 plus commutator corrections."""

    _ = U_fund
    if metadata_only:
        terms = _empty_like(S_adj)
        terms.metadata["sectors"]["KJJSSJ"] = {
            "normal_form": "K3+K2_comm+K1_comm",
            "include_commutators": include_commutators,
            "metadata_only": True,
            "norms": terms.norms(),
        }
        return terms
    A = kjjssj_A_LLR_from_kernel(S_adj, KJJSSJ, f)
    B = kjjssj_B_LRR_from_kernel(S_adj, KJJSSJ, f)
    V = kjjssj_V_virtual_from_kernel(KJJSSJ, f)
    raw_terms = {}
    for source in (
        cubic_block_terms_from_LLR(llr_right_to_left(A, S_adj)),
        cubic_block_terms_from_LRR(lrr_right_to_left(B, S_adj)),
        cubic_block_terms_from_virtual_LLL(V),
        cubic_block_terms_from_virtual_RRR(virtual_rrr_right_to_left(V, S_adj, sign=-1.0)),
    ):
        for word, coeff in source.items():
            raw_terms[word] = raw_terms.get(word, 0.0) + coeff
    return _assemble_cubic_from_raw_terms(S_adj, f, raw_terms, include_commutators, "KJJSSJ")


def assemble_nlo_current_terms(
    U_fund,
    S_adj,
    kernels: dict,
    gens,
    f,
    include_commutators: bool = True,
    metadata_only: bool = False,
) -> NLOCurrentTerms:
    """Assemble all currently validated NLO current pieces into dense K1/K2/K3."""

    total = _empty_like(S_adj)
    total.metadata["include_commutators"] = include_commutators
    total.metadata["commutators"] = {
        "included": include_commutators,
        "scope": "coincident-site lower-order corrections for validated cubic sectors",
    }
    assemblers = {
        "KJSJ": lambda: assemble_kjsj_terms(U_fund, S_adj, kernels["KJSJ"], metadata_only),
        "KJSSJ": lambda: assemble_kjssj_terms(U_fund, S_adj, kernels["KJSSJ"], f, metadata_only),
        "Kqbarq": lambda: assemble_kqbarq_terms(U_fund, S_adj, kernels["Kqbarq"], gens, metadata_only),
        "KJJSJ": lambda: assemble_kjjsj_terms(
            U_fund, S_adj, kernels["KJJSJ"], f, include_commutators, metadata_only
        ),
        "KJJSSJ": lambda: assemble_kjjssj_terms(
            U_fund, S_adj, kernels["KJJSSJ"], f, include_commutators, metadata_only
        ),
    }
    for key, assembler in assemblers.items():
        if key not in kernels:
            total.metadata["sectors"][key] = {"missing": True, "norms": {"K1": 0.0, "K2": 0.0, "K3": 0.0}}
            total.metadata["warnings"].append(f"missing kernel: {key}")
            continue
        total = add_terms(total, assembler())
    total.metadata["total_norms"] = total.norms()
    total.metadata["nonproduction_only"] = True
    return total
