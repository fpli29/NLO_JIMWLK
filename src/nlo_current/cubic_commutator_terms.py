"""Map small dense cubic blocks to canonical Lie-word commutator terms."""

from __future__ import annotations

from typing import Dict

import numpy as np

from .lie_word_algebra import Word, canonicalize_terms, split_by_order

DROP_TOL = 1e-14


def llr_right_to_left(A: np.ndarray, S_adj: np.ndarray) -> np.ndarray:
    """Convert A[x,y,w,d,e,a] J_L J_L J_R to C[x,y,w,d,e,h]."""

    return np.einsum("xywdea,wha->xywdeh", A, S_adj, optimize=True)


def lrr_right_to_left(B: np.ndarray, S_adj: np.ndarray) -> np.ndarray:
    """Convert B[w,x,y,a,d,e] J_L J_R J_R to C[w,x,y,a,p,q]."""

    return np.einsum("wxyade,xpd,yqe->wxyapq", B, S_adj, S_adj, optimize=True)


def virtual_rrr_right_to_left(V: np.ndarray, S_adj: np.ndarray, sign: float = -1.0) -> np.ndarray:
    """Convert V[x,y,w,c,b,a] R_x R_y R_w to left-basis C[x,y,w,p,q,h].

    The default sign is -1 because the virtual Hamiltonian block is
    V(LLL - RRR).
    """

    return sign * np.einsum("xywcba,xpc,yqb,wha->xywpqh", V, S_adj, S_adj, S_adj, optimize=True)


def _add(out: Dict[Word, complex], word: Word, coeff) -> None:
    if abs(coeff) < DROP_TOL:
        return
    out[word] = out.get(word, 0.0) + coeff
    if abs(out[word]) < DROP_TOL:
        del out[word]


def cubic_block_terms_from_LLR(C_left: np.ndarray) -> Dict[Word, complex]:
    """Convert C[x,y,w,d,e,h] into word terms ((x,d),(y,e),(w,h))."""

    terms: Dict[Word, complex] = {}
    nsite, _, _, n_color, _, _ = C_left.shape
    for x in range(nsite):
        for y in range(nsite):
            for w in range(nsite):
                for d in range(n_color):
                    for e in range(n_color):
                        for h in range(n_color):
                            _add(terms, ((x, d), (y, e), (w, h)), C_left[x, y, w, d, e, h])
    return terms


def cubic_block_terms_from_LRR(C_left: np.ndarray) -> Dict[Word, complex]:
    """Convert C[w,x,y,a,p,q] into word terms ((w,a),(x,p),(y,q))."""

    terms: Dict[Word, complex] = {}
    nsite, _, _, n_color, _, _ = C_left.shape
    for w in range(nsite):
        for x in range(nsite):
            for y in range(nsite):
                for a in range(n_color):
                    for p in range(n_color):
                        for q in range(n_color):
                            _add(terms, ((w, a), (x, p), (y, q)), C_left[w, x, y, a, p, q])
    return terms


def cubic_block_terms_from_virtual_LLL(V: np.ndarray) -> Dict[Word, complex]:
    """Convert V[x,y,w,c,b,a] into word terms ((x,c),(y,b),(w,a))."""

    terms: Dict[Word, complex] = {}
    nsite, _, _, n_color, _, _ = V.shape
    for x in range(nsite):
        for y in range(nsite):
            for w in range(nsite):
                for c in range(n_color):
                    for b in range(n_color):
                        for a in range(n_color):
                            _add(terms, ((x, c), (y, b), (w, a)), V[x, y, w, c, b, a])
    return terms


def cubic_block_terms_from_virtual_RRR(V_left: np.ndarray) -> Dict[Word, complex]:
    """Convert left-basis RRR C[x,y,w,p,q,h] into word terms."""

    terms: Dict[Word, complex] = {}
    nsite, _, _, n_color, _, _ = V_left.shape
    for x in range(nsite):
        for y in range(nsite):
            for w in range(nsite):
                for p in range(n_color):
                    for q in range(n_color):
                        for h in range(n_color):
                            _add(terms, ((x, p), (y, q), (w, h)), V_left[x, y, w, p, q, h])
    return terms


def canonicalize_cubic_block_terms(terms: Dict[Word, float], f: np.ndarray):
    """Canonicalize cubic block terms and return cubic, quadratic, and linear parts."""

    canonical = canonicalize_terms(terms, f)
    cubic_terms, quadratic_terms, linear_terms, scalar_terms = split_by_order(canonical)
    if scalar_terms:
        raise ValueError("unexpected scalar terms from cubic commutator canonicalization")
    return cubic_terms, quadratic_terms, linear_terms
