"""Synthetic tiny-lattice kernels for non-production NLO current tests."""

from __future__ import annotations

import numpy as np

from .three_generator_terms import synthetic_kjjsj_kernel, synthetic_kjjssj_kernel


def _normalize(kernel):
    norm = np.linalg.norm(kernel)
    return kernel if norm == 0.0 else kernel / norm


def synthetic_kjsj_kernel(nsite, rng):
    raw = rng.normal(size=(nsite, nsite, nsite))
    return _normalize(0.5 * (raw + np.swapaxes(raw, 0, 1)))


def synthetic_kjssj_kernel(nsite, rng):
    raw = rng.normal(size=(nsite, nsite, nsite, nsite))
    kernel = 0.5 * (raw + np.swapaxes(raw, 0, 1))
    kernel = 0.5 * (kernel + np.swapaxes(kernel, 2, 3))
    return _normalize(kernel)


def synthetic_kqbarq_kernel(nsite, rng):
    raw = rng.normal(size=(nsite, nsite, nsite, nsite))
    kernel = 0.5 * (raw + np.swapaxes(raw, 0, 1))
    kernel = 0.5 * (kernel + np.swapaxes(kernel, 2, 3))
    return _normalize(kernel)


def synthetic_kernels_all(nsite, rng):
    """Return all synthetic kernels used by the non-production skeleton."""

    return {
        "KJSJ": synthetic_kjsj_kernel(nsite, rng),
        "KJSSJ": synthetic_kjssj_kernel(nsite, rng),
        "Kqbarq": synthetic_kqbarq_kernel(nsite, rng),
        "KJJSJ": synthetic_kjjsj_kernel(nsite, rng, xy_symmetry="antisymmetric"),
        "KJJSSJ": synthetic_kjjssj_kernel(nsite, rng, klm_antisym=True),
    }

