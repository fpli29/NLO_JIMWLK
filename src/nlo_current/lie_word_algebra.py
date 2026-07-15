"""Small symbolic Lie-derivative word algebra for cubic commutator checks.

Labels are ``(site, color)`` pairs and words are tuples of labels. The
canonical order is lexicographic. For the implemented left perturbation
convention,

    [L_x^a, L_x^b] = f^{abc} L_x^c,

so an adjacent same-site swap uses

    L_x^a L_x^b = L_x^b L_x^a + f^{abc} L_x^c.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Tuple

import numpy as np

Label = Tuple[int, int]
Word = Tuple[Label, ...]

DROP_TOL = 1e-14


def is_canonical_word(word: Word) -> bool:
    """Return True if the derivative word is in lexicographic order."""

    return all(word[i] <= word[i + 1] for i in range(len(word) - 1))


def _add_term(out: Dict[Word, complex], word: Word, coeff) -> None:
    if abs(coeff) < DROP_TOL:
        return
    out[word] = out.get(word, 0.0) + coeff
    if abs(out[word]) < DROP_TOL:
        del out[word]


def canonicalize_word(word: Word, f: np.ndarray, n_color: int = 8) -> Dict[Word, complex]:
    """Canonicalize a word of length 1, 2, or 3.

    The sign convention is

        L_x^a L_x^b = L_x^b L_x^a + f^{abc} L_x^c

    when swapping an adjacent same-site pair. Different-site swaps commute with
    no lower-order term. The algorithm recursively canonicalizes until all
    words are ordered and drops coefficients with absolute value below 1e-14.
    """

    word = tuple(word)

    @lru_cache(maxsize=None)
    def _canonicalize(cached_word: Word) -> Tuple[Tuple[Word, complex], ...]:
        if len(cached_word) <= 1 or is_canonical_word(cached_word):
            return ((cached_word, 1.0),)

        swap_index = None
        for i in range(len(cached_word) - 1):
            if cached_word[i] > cached_word[i + 1]:
                swap_index = i
                break
        if swap_index is None:
            return ((cached_word, 1.0),)

        i = swap_index
        left = cached_word[:i]
        right = cached_word[i + 2 :]
        first = cached_word[i]
        second = cached_word[i + 1]
        swapped = left + (second, first) + right

        combined: Dict[Word, complex] = {}
        for sub_word, coeff in _canonicalize(swapped):
            _add_term(combined, sub_word, coeff)

        site_a, color_a = first
        site_b, color_b = second
        if site_a == site_b:
            for color_c in range(n_color):
                coeff = float(f[color_a, color_b, color_c])
                if abs(coeff) >= DROP_TOL:
                    comm_word = left + ((site_a, color_c),) + right
                    for sub_word, sub_coeff in _canonicalize(comm_word):
                        _add_term(combined, sub_word, coeff * sub_coeff)

        return tuple(sorted(combined.items(), key=lambda item: item[0]))

    return dict(_canonicalize(word))


def canonicalize_terms(
    terms: Dict[Word, complex],
    f: np.ndarray,
    n_color: int = 8,
) -> Dict[Word, complex]:
    """Canonicalize and combine a dictionary of word coefficients."""

    out: Dict[Word, complex] = {}
    for word, coeff in terms.items():
        if abs(coeff) < DROP_TOL:
            continue
        for canonical_word, canonical_coeff in canonicalize_word(word, f, n_color).items():
            _add_term(out, canonical_word, coeff * canonical_coeff)
    return dict(sorted(out.items(), key=lambda item: item[0]))


def split_by_order(terms: Dict[Word, complex]):
    """Split terms into cubic, quadratic, linear, and scalar dictionaries."""

    cubic_terms: Dict[Word, complex] = {}
    quadratic_terms: Dict[Word, complex] = {}
    linear_terms: Dict[Word, complex] = {}
    scalar_terms: Dict[Word, complex] = {}
    for word, coeff in terms.items():
        if len(word) == 3:
            cubic_terms[word] = coeff
        elif len(word) == 2:
            quadratic_terms[word] = coeff
        elif len(word) == 1:
            linear_terms[word] = coeff
        elif len(word) == 0:
            scalar_terms[word] = coeff
        else:
            raise ValueError(f"unsupported word length {len(word)} for {word}")
    return cubic_terms, quadratic_terms, linear_terms, scalar_terms
