from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.lie_word_algebra import (  # noqa: E402
    canonicalize_word,
    is_canonical_word,
    split_by_order,
)
from nlo_current.su3_adjoint import structure_constants, su3_generators_fundamental  # noqa: E402


def _assert_terms_close(left, right, atol=1e-12):
    assert set(left) == set(right)
    for word in left:
        assert abs(left[word] - right[word]) < atol


def _canonicalize_last_inversion(word, f, n_color=8):
    if is_canonical_word(word):
        return {word: 1.0}
    inversions = [i for i in range(len(word) - 1) if word[i] > word[i + 1]]
    i = inversions[-1]
    first = word[i]
    second = word[i + 1]
    left = word[:i]
    right = word[i + 2 :]
    out = {}

    def add_terms(terms, factor=1.0):
        for term_word, coeff in terms.items():
            out[term_word] = out.get(term_word, 0.0) + factor * coeff

    add_terms(_canonicalize_last_inversion(left + (second, first) + right, f, n_color))
    if first[0] == second[0]:
        for c in range(n_color):
            coeff = float(f[first[1], second[1], c])
            if abs(coeff) > 1e-14:
                add_terms(_canonicalize_last_inversion(left + ((first[0], c),) + right, f, n_color), coeff)
    return {word: coeff for word, coeff in out.items() if abs(coeff) > 1e-14}


def test_different_site_derivatives_commute_without_lower_order_term() -> None:
    f = structure_constants(su3_generators_fundamental())
    terms = canonicalize_word(((1, 3), (0, 5)), f)
    assert terms == {((0, 5), (1, 3)): 1.0}


def test_same_site_two_word_commutator() -> None:
    f = structure_constants(su3_generators_fundamental())
    terms = canonicalize_word(((0, 2), (0, 1)), f)
    expected = {((0, 1), (0, 2)): 1.0}
    for c in range(8):
        if abs(f[2, 1, c]) > 1e-14:
            expected[((0, c),)] = float(f[2, 1, c])
    _assert_terms_close(terms, expected)


def test_three_same_site_words_are_path_independent() -> None:
    f = structure_constants(su3_generators_fundamental())
    word = ((0, 4), (0, 2), (0, 1))
    first_path = canonicalize_word(word, f)
    last_path = _canonicalize_last_inversion(word, f)
    _assert_terms_close(first_path, last_path, atol=5e-13)


def test_already_canonical_word_unchanged() -> None:
    f = structure_constants(su3_generators_fundamental())
    word = ((0, 1), (0, 2), (1, 0))
    assert is_canonical_word(word)
    assert canonicalize_word(word, f) == {word: 1.0}


def test_split_by_order() -> None:
    terms = {
        ((0, 0), (0, 1), (1, 2)): 1.0,
        ((0, 0), (1, 2)): 2.0,
        ((1, 2),): 3.0,
        tuple(): 4.0,
    }
    cubic, quadratic, linear, scalar = split_by_order(terms)
    assert cubic == {((0, 0), (0, 1), (1, 2)): 1.0}
    assert quadratic == {((0, 0), (1, 2)): 2.0}
    assert linear == {((1, 2),): 3.0}
    assert scalar == {tuple(): 4.0}

