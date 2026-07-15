"""Coordinate helpers for non-production physical-kernel diagnostics."""

from __future__ import annotations

import numpy as np


class KernelSingularityError(ZeroDivisionError):
    """Raised when an exact singular denominator is encountered."""


def validate_coords(coords: np.ndarray) -> np.ndarray:
    """Validate and return coordinates as a float array with shape (Nsite, 2)."""

    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"coords must have shape (Nsite, 2), got {arr.shape}")
    if arr.shape[0] == 0:
        raise ValueError("coords must contain at least one site")
    if not np.all(np.isfinite(arr)):
        raise ValueError("coords must be finite")
    return arr


def vec(coords: np.ndarray, i: int, j: int) -> np.ndarray:
    """Return coords[i] - coords[j]."""

    arr = validate_coords(coords)
    return arr[int(i)] - arr[int(j)]


def norm2(v: np.ndarray) -> float:
    """Return squared Euclidean norm."""

    a = np.asarray(v, dtype=float)
    return float(np.dot(a, a))


def dot(a: np.ndarray, b: np.ndarray) -> float:
    """Return Euclidean dot product."""

    return float(np.dot(np.asarray(a, dtype=float), np.asarray(b, dtype=float)))


def cross2(a: np.ndarray, b: np.ndarray) -> float:
    """Return the 2D scalar cross product a_x b_y - a_y b_x."""

    av = np.asarray(a, dtype=float)
    bv = np.asarray(b, dtype=float)
    if av.shape != (2,) or bv.shape != (2,):
        raise ValueError("cross2 expects two vectors of shape (2,)")
    return float(av[0] * bv[1] - av[1] * bv[0])


def safe_inv(
    x,
    singularity_policy: str = "raise",
    eps: float | None = None,
    name: str = "denominator",
):
    """Return 1/x with explicit exact-zero singularity handling."""

    arr = np.asarray(x, dtype=float)
    singular = arr == 0.0
    if singularity_policy == "raise":
        if np.any(singular):
            raise KernelSingularityError(f"singular {name}: exact zero denominator")
        out = 1.0 / arr
    elif singularity_policy == "nan":
        out = np.empty_like(arr, dtype=float)
        nonsingular = ~singular
        out[nonsingular] = 1.0 / arr[nonsingular]
        out[singular] = np.nan
    elif singularity_policy == "eps":
        if eps is None or eps <= 0.0:
            raise ValueError("eps policy requires a positive eps argument")
        regularized = np.where(singular, float(eps), arr)
        out = 1.0 / regularized
    else:
        raise ValueError("singularity_policy must be 'raise', 'nan', or 'eps'")
    if np.isscalar(x):
        return float(out)
    return out


def pairwise_dist2(coords: np.ndarray) -> np.ndarray:
    """Return matrix r_ij^2."""

    arr = validate_coords(coords)
    diff = arr[:, None, :] - arr[None, :, :]
    return np.einsum("ijk,ijk->ij", diff, diff)
