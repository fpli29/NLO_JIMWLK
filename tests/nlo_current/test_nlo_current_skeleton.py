from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.nlo_current_skeleton import (  # noqa: E402
    NLOCurrentTerms,
    add_terms,
    combined_dim,
    empty_terms,
)
from nlo_current.nlo_velocity_evaluator import (  # noqa: E402
    cubic_density_contraction,
    evaluate_velocity_from_terms,
)


def test_empty_terms_shape_for_three_sites() -> None:
    terms = empty_terms(nsite=3)

    assert combined_dim(3) == 24
    assert terms.dim == 24
    assert terms.K1.shape == (24,)
    assert terms.K2.shape == (24, 24)
    assert terms.K3.shape == (24, 24, 24)
    terms.validate_shapes()


def test_add_terms_adds_arrays_and_merges_metadata() -> None:
    rng = np.random.default_rng(20260704)
    dim = 4
    lhs = NLOCurrentTerms(
        K1=rng.normal(size=dim),
        K2=rng.normal(size=(dim, dim)),
        K3=rng.normal(size=(dim, dim, dim)),
        metadata={"sectors": {"lhs": {"norms": {}}}, "warnings": ["left"], "sources": ["a"]},
    )
    rhs = NLOCurrentTerms(
        K1=rng.normal(size=dim),
        K2=rng.normal(size=(dim, dim)),
        K3=rng.normal(size=(dim, dim, dim)),
        metadata={"sectors": {"rhs": {"norms": {}}}, "warnings": ["right"], "sources": ["b"]},
    )

    out = add_terms(lhs, rhs)

    np.testing.assert_allclose(out.K1, lhs.K1 + rhs.K1)
    np.testing.assert_allclose(out.K2, lhs.K2 + rhs.K2)
    np.testing.assert_allclose(out.K3, lhs.K3 + rhs.K3)
    assert set(out.metadata["sectors"]) == {"lhs", "rhs"}
    assert out.metadata["warnings"] == ["left", "right"]
    assert out.metadata["sources"] == ["a", "b"]
    out.validate_shapes()


def test_velocity_evaluator_output_shape_with_explicit_derivatives() -> None:
    rng = np.random.default_rng(31001)
    dim = 5
    terms = NLOCurrentTerms(
        K1=rng.normal(size=dim),
        K2=rng.normal(size=(dim, dim)),
        K3=rng.normal(size=(dim, dim, dim)),
        metadata={},
    )
    score = rng.normal(size=dim)
    hessian_score = rng.normal(size=(dim, dim))
    dK2 = rng.normal(size=dim)
    d2K3 = rng.normal(size=dim)
    dK3_first = {
        "LC_K3_ABC": rng.normal(size=(dim, dim)),
        "LB_K3_ABC": rng.normal(size=(dim, dim)),
    }

    result = evaluate_velocity_from_terms(
        terms,
        score,
        hessian_score,
        dK2=dK2,
        dK3_first=dK3_first,
        d2K3=d2K3,
    )

    assert result["velocity"].shape == (dim,)
    assert result["diagnostics"]["warnings"] == []


def test_cubic_density_contraction_matches_explicit_loops() -> None:
    rng = np.random.default_rng(31002)
    dim = 3
    K3 = rng.normal(size=(dim, dim, dim))
    score = rng.normal(size=dim)
    hessian_score = rng.normal(size=(dim, dim))

    expected = np.zeros(dim)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                expected[a] += K3[a, b, c] * (
                    hessian_score[b, c] + score[b] * score[c]
                )

    np.testing.assert_allclose(
        cubic_density_contraction(K3, score, hessian_score),
        expected,
        atol=1e-12,
        rtol=1e-12,
    )


def test_velocity_evaluator_warns_when_derivative_arrays_are_omitted() -> None:
    rng = np.random.default_rng(31003)
    dim = 4
    terms = NLOCurrentTerms(
        K1=rng.normal(size=dim),
        K2=rng.normal(size=(dim, dim)),
        K3=rng.normal(size=(dim, dim, dim)),
        metadata={},
    )

    result = evaluate_velocity_from_terms(
        terms,
        score=rng.normal(size=dim),
        hessian_score=rng.normal(size=(dim, dim)),
    )

    assert result["velocity"].shape == (dim,)
    assert result["diagnostics"]["coefficient_derivatives_omitted"] is True
    warnings = result["diagnostics"]["warnings"]
    assert any("dK2 omitted" in warning for warning in warnings)
    assert any("dK3_first omitted" in warning for warning in warnings)
    assert any("d2K3 omitted" in warning for warning in warnings)
