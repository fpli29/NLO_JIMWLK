"""Dense finite-difference coefficient derivatives for tiny diagnostics."""

from __future__ import annotations

import numpy as np

from .nlo_velocity_evaluator import cubic_density_contraction, evaluate_velocity_from_terms
from .su3_adjoint import left_perturb


def validate_coefficient_shapes(K2, K3):
    """Validate K2 shape (D,D) and K3 shape (D,D,D), with matching D."""

    K2 = np.asarray(K2)
    K3 = np.asarray(K3)
    if K2.ndim != 2 or K2.shape[0] != K2.shape[1]:
        raise ValueError(f"K2 must have shape (D,D), got {K2.shape}")
    dim = K2.shape[0]
    if K3.shape != (dim, dim, dim):
        raise ValueError(f"K3 must have shape {(dim, dim, dim)}, got {K3.shape}")
    return dim


def _validate_K2_shape(K2, expected_dim=None):
    K2 = np.asarray(K2)
    if K2.ndim != 2 or K2.shape[0] != K2.shape[1]:
        raise ValueError(f"K2 must have shape (D,D), got {K2.shape}")
    if expected_dim is not None and K2.shape != (expected_dim, expected_dim):
        raise ValueError(f"K2 must have shape {(expected_dim, expected_dim)}, got {K2.shape}")
    return K2.shape[0]


def _validate_K3_shape(K3, expected_dim=None):
    K3 = np.asarray(K3)
    if K3.ndim != 3 or K3.shape[0] != K3.shape[1] or K3.shape[0] != K3.shape[2]:
        raise ValueError(f"K3 must have shape (D,D,D), got {K3.shape}")
    if expected_dim is not None and K3.shape != (expected_dim, expected_dim, expected_dim):
        raise ValueError(
            f"K3 must have shape {(expected_dim, expected_dim, expected_dim)}, got {K3.shape}"
        )
    return K3.shape[0]


def _combined_dim_from_U(U_fund, gens):
    return int(len(U_fund)) * int(gens.shape[0])


def _unflatten(index, n_color):
    return int(index) // int(n_color), int(index) % int(n_color)


def left_perturbed_copy(U_fund, site, color, gens, eps):
    """Return a copy of U_fund with one Wilson line left-perturbed."""

    out = np.array(U_fund, copy=True)
    out[site] = left_perturb(out[site], color, eps, gens)
    return out


def fd_left_derivative_array(callback, U_fund, S_adj_builder, site, color, gens, eps):
    """Central finite-difference left derivative of an array-valued callback."""

    plus = left_perturbed_copy(U_fund, site, color, gens, eps)
    minus = left_perturbed_copy(U_fund, site, color, gens, -eps)
    plus_value = np.asarray(callback(plus, S_adj_builder(plus)))
    minus_value = np.asarray(callback(minus, S_adj_builder(minus)))
    return (plus_value - minus_value) / (2.0 * eps)


def _fd_second_left_derivative_array(
    callback,
    U_fund,
    S_adj_builder,
    outer_site,
    outer_color,
    inner_site,
    inner_color,
    gens,
    eps,
):
    """Return ordered L_outer L_inner callback by nested central differences."""

    return fd_left_derivative_array(
        lambda V, _S: fd_left_derivative_array(
            callback,
            V,
            S_adj_builder,
            inner_site,
            inner_color,
            gens,
            eps,
        ),
        U_fund,
        S_adj_builder,
        outer_site,
        outer_color,
        gens,
        eps,
    )


def compute_dK2_fd(K2_callback, U_fund, S_adj_builder, gens, eps=1e-5):
    """Compute dK2^A = sum_B L_B K2^{A B} by dense finite differences."""

    n_color = gens.shape[0]
    dim = _combined_dim_from_U(U_fund, gens)
    base = np.asarray(K2_callback(U_fund, S_adj_builder(U_fund)))
    _validate_K2_shape(base, dim)
    dK2 = np.zeros(dim, dtype=np.result_type(base, float))
    for b_index in range(dim):
        site, color = _unflatten(b_index, n_color)
        derivative = fd_left_derivative_array(
            K2_callback,
            U_fund,
            S_adj_builder,
            site,
            color,
            gens,
            eps,
        )
        _validate_K2_shape(derivative, dim)
        dK2 += derivative[:, b_index]
    return np.real_if_close(dK2).real


def compute_dK3_first_fd(K3_callback, U_fund, S_adj_builder, gens, eps=1e-5):
    """Compute first-derivative contractions of K3 by dense finite differences."""

    n_color = gens.shape[0]
    dim = _combined_dim_from_U(U_fund, gens)
    base = np.asarray(K3_callback(U_fund, S_adj_builder(U_fund)))
    _validate_K3_shape(base, dim)
    LC_K3_ABC = np.zeros((dim, dim), dtype=np.result_type(base, float))
    LB_K3_ABC = np.zeros((dim, dim), dtype=np.result_type(base, float))
    for index in range(dim):
        site, color = _unflatten(index, n_color)
        derivative = fd_left_derivative_array(
            K3_callback,
            U_fund,
            S_adj_builder,
            site,
            color,
            gens,
            eps,
        )
        _validate_K3_shape(derivative, dim)
        LC_K3_ABC += derivative[:, :, index]
        LB_K3_ABC += derivative[:, index, :]
    return {
        "LC_K3_ABC": np.real_if_close(LC_K3_ABC).real,
        "LB_K3_ABC": np.real_if_close(LB_K3_ABC).real,
    }


def compute_d2K3_fd(K3_callback, U_fund, S_adj_builder, gens, eps=1e-4):
    """Compute d2K3^A = sum_{B,C} L_B L_C K3^{A B C}.

    This nested finite difference is intentionally dense and expensive. It is
    only meant for tiny diagnostic lattices.
    """

    n_color = gens.shape[0]
    dim = _combined_dim_from_U(U_fund, gens)
    base = np.asarray(K3_callback(U_fund, S_adj_builder(U_fund)))
    _validate_K3_shape(base, dim)
    d2K3 = np.zeros(dim, dtype=np.result_type(base, float))
    for b_index in range(dim):
        b_site, b_color = _unflatten(b_index, n_color)
        for c_index in range(dim):
            c_site, c_color = _unflatten(c_index, n_color)
            derivative = _fd_second_left_derivative_array(
                K3_callback,
                U_fund,
                S_adj_builder,
                b_site,
                b_color,
                c_site,
                c_color,
                gens,
                eps,
            )
            _validate_K3_shape(derivative, dim)
            d2K3 += derivative[:, b_index, c_index]
    return np.real_if_close(d2K3).real


def compute_all_coefficient_derivatives_fd(
    K2_callback,
    K3_callback,
    U_fund,
    S_adj_builder,
    gens,
    eps_first=1e-5,
    eps_second=1e-4,
):
    """Return all dense coefficient-derivative arrays for the velocity evaluator."""

    return {
        "dK2": compute_dK2_fd(K2_callback, U_fund, S_adj_builder, gens, eps=eps_first),
        "dK3_first": compute_dK3_first_fd(
            K3_callback,
            U_fund,
            S_adj_builder,
            gens,
            eps=eps_first,
        ),
        "d2K3": compute_d2K3_fd(K3_callback, U_fund, S_adj_builder, gens, eps=eps_second),
    }


def product_rule_K2_rhs(K2, dK2, score):
    """Return dK2^A + K2^{AB}s_B."""

    dim = _validate_K2_shape(K2)
    dK2 = np.asarray(dK2)
    score = np.asarray(score)
    if dK2.shape != (dim,) or score.shape != (dim,):
        raise ValueError("dK2 and score must have shape (D,)")
    return dK2 + np.einsum("ab,b->a", K2, score, optimize=True)


def product_rule_K3_rhs(K3, dK3_first, d2K3, score, hessian_score):
    """Return the K3 product-rule contraction used in the diagnostic velocity."""

    dim = _validate_K3_shape(K3)
    d2K3 = np.asarray(d2K3)
    score = np.asarray(score)
    hessian_score = np.asarray(hessian_score)
    LC_K3_ABC = np.asarray(dK3_first.get("LC_K3_ABC"))
    LB_K3_ABC = np.asarray(dK3_first.get("LB_K3_ABC"))
    if d2K3.shape != (dim,) or score.shape != (dim,):
        raise ValueError("d2K3 and score must have shape (D,)")
    if hessian_score.shape != (dim, dim):
        raise ValueError("hessian_score must have shape (D,D)")
    if LC_K3_ABC.shape != (dim, dim) or LB_K3_ABC.shape != (dim, dim):
        raise ValueError("dK3_first arrays must have shape (D,D)")
    return (
        d2K3
        + np.einsum("ab,b->a", LC_K3_ABC, score, optimize=True)
        + np.einsum("ac,c->a", LB_K3_ABC, score, optimize=True)
        + cubic_density_contraction(K3, score, hessian_score)
    )


def velocity_from_coeff_derivative_backend(terms, score, hessian_score, derivatives):
    """Evaluate velocity using coefficient derivatives returned by this backend."""

    return evaluate_velocity_from_terms(
        terms,
        score,
        hessian_score,
        dK2=derivatives["dK2"],
        dK3_first=derivatives["dK3_first"],
        d2K3=derivatives["d2K3"],
    )
