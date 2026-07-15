"""Diagnostic velocity evaluators for dense non-production NLO current terms."""

from __future__ import annotations

import numpy as np


def cubic_density_contraction(K3, score, hessian_score):
    """Return K3^{ABC}(H_BC + s_B s_C)."""

    score = np.asarray(score)
    hessian_score = np.asarray(hessian_score)
    if K3.shape[1:] != hessian_score.shape:
        raise ValueError("K3 and hessian_score shapes are incompatible")
    if score.shape != (K3.shape[0],):
        raise ValueError("score shape is incompatible with K3")
    density_pair = hessian_score + np.outer(score, score)
    return np.einsum("abc,bc->a", K3, density_pair, optimize=True)


def evaluate_velocity_from_terms(
    terms,
    score,
    hessian_score,
    dK2=None,
    dK3_first=None,
    d2K3=None,
):
    """Evaluate diagnostic velocity from explicit dense derivative arrays.

    Missing coefficient derivative arrays are treated as zero and recorded in
    the returned diagnostics. This is not a complete production evaluator.
    """

    terms.validate_shapes()
    dim = terms.dim
    score = np.asarray(score)
    hessian_score = np.asarray(hessian_score)
    if score.shape != (dim,):
        raise ValueError(f"score must have shape {(dim,)}, got {score.shape}")
    if hessian_score.shape != (dim, dim):
        raise ValueError(f"hessian_score must have shape {(dim, dim)}, got {hessian_score.shape}")

    warnings = []
    if dK2 is None:
        dK2 = np.zeros(dim)
        warnings.append("dK2 omitted; treated as zero")
    if d2K3 is None:
        d2K3 = np.zeros(dim)
        warnings.append("d2K3 omitted; treated as zero")
    if dK3_first is None:
        LC_K3_ABC = np.zeros((dim, dim))
        LB_K3_ABC = np.zeros((dim, dim))
        warnings.append("dK3_first omitted; first coefficient derivatives treated as zero")
    else:
        LC_K3_ABC = np.asarray(dK3_first.get("LC_K3_ABC"))
        LB_K3_ABC = np.asarray(dK3_first.get("LB_K3_ABC"))
        if LC_K3_ABC.shape != (dim, dim) or LB_K3_ABC.shape != (dim, dim):
            raise ValueError("dK3_first arrays must have shape (D,D)")

    dK2 = np.asarray(dK2)
    d2K3 = np.asarray(d2K3)
    if dK2.shape != (dim,) or d2K3.shape != (dim,):
        raise ValueError("dK2 and d2K3 must have shape (D,)")

    velocity = (
        terms.K1
        - 0.5 * (dK2 + np.einsum("ab,b->a", terms.K2, score, optimize=True))
        + (1.0 / 6.0)
        * (
            d2K3
            + np.einsum("ab,b->a", LC_K3_ABC, score, optimize=True)
            + np.einsum("ac,c->a", LB_K3_ABC, score, optimize=True)
            + cubic_density_contraction(terms.K3, score, hessian_score)
        )
    )
    return {
        "velocity": velocity,
        "diagnostics": {
            "warnings": warnings,
            "coefficient_derivatives_omitted": bool(warnings),
            "velocity_norm": float(np.linalg.norm(velocity)),
        },
    }


def evaluate_velocity_score_only(terms, score):
    """Convenience diagnostic ignoring coefficient derivatives and Hessian-score."""

    terms.validate_shapes()
    score = np.asarray(score)
    if score.shape != (terms.dim,):
        raise ValueError("score shape is incompatible with terms")
    cubic_score = np.einsum("abc,b,c->a", terms.K3, score, score, optimize=True)
    return terms.K1 - 0.5 * terms.K2 @ score + (1.0 / 6.0) * cubic_score

