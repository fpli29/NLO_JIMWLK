r"""Dense SU(3) helpers for small-lattice NLO current checks.

The Lie derivative convention used here is

    L^a F(U) = d/deps F(exp(i eps t^a) U)|_{eps=0}
    R^a F(U) = d/deps F(U exp(i eps t^a))|_{eps=0}

with fundamental generators normalized by Tr(t^a t^b) = delta^{ab}/2.
The adjoint Wilson line is

    S_A^{ab}(U) = 2 Re Tr(t^a U t^b U^\dagger).

With these conventions, J_L^a = S_A^{ab} J_R^b and
J_R^a = S_A^{ba} J_L^b.
"""

from __future__ import annotations

import numpy as np


def su3_generators_fundamental() -> np.ndarray:
    """Return SU(3) fundamental generators t^a = lambda^a/2."""

    zero = 0.0
    one = 1.0
    i = 1.0j
    sqrt3 = np.sqrt(3.0)

    lambdas = [
        np.array([[zero, one, zero], [one, zero, zero], [zero, zero, zero]], dtype=complex),
        np.array([[zero, -i, zero], [i, zero, zero], [zero, zero, zero]], dtype=complex),
        np.array([[one, zero, zero], [zero, -one, zero], [zero, zero, zero]], dtype=complex),
        np.array([[zero, zero, one], [zero, zero, zero], [one, zero, zero]], dtype=complex),
        np.array([[zero, zero, -i], [zero, zero, zero], [i, zero, zero]], dtype=complex),
        np.array([[zero, zero, zero], [zero, zero, one], [zero, one, zero]], dtype=complex),
        np.array([[zero, zero, zero], [zero, zero, -i], [zero, i, zero]], dtype=complex),
        (1.0 / sqrt3) * np.array(
            [[one, zero, zero], [zero, one, zero], [zero, zero, -2.0]],
            dtype=complex,
        ),
    ]
    return 0.5 * np.stack(lambdas, axis=0)


def structure_constants(gens: np.ndarray) -> np.ndarray:
    """Return f[a,b,c] from [t^a,t^b] = i f^{abc} t^c."""

    n_gen = gens.shape[0]
    f = np.empty((n_gen, n_gen, n_gen), dtype=float)
    for a in range(n_gen):
        for b in range(n_gen):
            comm = gens[a] @ gens[b] - gens[b] @ gens[a]
            for c in range(n_gen):
                f[a, b, c] = np.real(-2.0j * np.trace(comm @ gens[c]))
    return f


def adjoint_from_fundamental(U: np.ndarray, gens: np.ndarray) -> np.ndarray:
    r"""Return S_A^{ab} = 2 Re Tr(t^a U t^b U^\dagger)."""

    Udag = U.conj().T
    n_gen = gens.shape[0]
    S = np.empty((n_gen, n_gen), dtype=float)
    for a in range(n_gen):
        for b in range(n_gen):
            S[a, b] = np.real(2.0 * np.trace(gens[a] @ U @ gens[b] @ Udag))
    return S


def random_su3(rng: np.random.Generator) -> np.ndarray:
    """Generate a random SU(3) matrix using QR projection."""

    z = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    q, r = np.linalg.qr(z)
    diag = np.diag(r)
    phases = diag / np.where(np.abs(diag) == 0.0, 1.0, np.abs(diag))
    q = q @ np.diag(phases.conj())
    q[:, 0] /= np.linalg.det(q)
    return q


def _fundamental_exponential(a: int, eps: float, gens: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eigh(gens[a])
    return (vecs * np.exp(1.0j * eps * vals)) @ vecs.conj().T


def left_perturb(U: np.ndarray, a: int, eps: float, gens: np.ndarray) -> np.ndarray:
    """Return exp(i eps t^a) U."""

    return _fundamental_exponential(a, eps, gens) @ U


def right_perturb(U: np.ndarray, a: int, eps: float, gens: np.ndarray) -> np.ndarray:
    """Return U exp(i eps t^a)."""

    return U @ _fundamental_exponential(a, eps, gens)


def finite_diff_left_derivative(
    func,
    U_list: list[np.ndarray],
    site: int,
    a: int,
    eps: float = 1e-6,
    gens: np.ndarray | None = None,
):
    """Central finite-difference left Lie derivative at one lattice site."""

    if gens is None:
        gens = su3_generators_fundamental()
    plus = [np.array(U, copy=True) for U in U_list]
    minus = [np.array(U, copy=True) for U in U_list]
    plus[site] = left_perturb(plus[site], a, eps, gens)
    minus[site] = left_perturb(minus[site], a, -eps, gens)
    return (func(plus) - func(minus)) / (2.0 * eps)


def finite_diff_right_derivative(
    func,
    U_list: list[np.ndarray],
    site: int,
    a: int,
    eps: float = 1e-6,
    gens: np.ndarray | None = None,
):
    """Central finite-difference right Lie derivative at one lattice site."""

    if gens is None:
        gens = su3_generators_fundamental()
    plus = [np.array(U, copy=True) for U in U_list]
    minus = [np.array(U, copy=True) for U in U_list]
    plus[site] = right_perturb(plus[site], a, eps, gens)
    minus[site] = right_perturb(minus[site], a, -eps, gens)
    return (func(plus) - func(minus)) / (2.0 * eps)
