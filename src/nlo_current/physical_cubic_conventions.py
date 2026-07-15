"""Explicit cubic-kernel convention helpers for physical diagnostics."""

from __future__ import annotations

import numpy as np


def klm_normalized_cubic_kernel(raw_kernel, *, require_real: bool = True, atol: float = 1.0e-12):
    """Return the KLM-normalized cubic coefficient array.

    The physical KLM formulas for K_JJSJ and K_JJSSJ carry an explicit ``-i``.
    The previously calibrated observable convention is

        TeX target = (-i) * raw Hermitian-generator direct action.

    Since the coefficient assembly is linear in the cubic kernel, the same
    convention layer maps raw physical cubic kernels to real normal-form
    coefficients by multiplying by ``-i``. If ``require_real`` is true, a
    non-negligible imaginary remainder is reported as an error instead of being
    silently discarded.
    """

    normalized = -1.0j * np.asarray(raw_kernel)
    if not require_real:
        return normalized
    max_imag = float(np.max(np.abs(np.imag(normalized)))) if normalized.size else 0.0
    if max_imag > atol:
        raise ValueError(
            "KLM-normalized cubic kernel has non-negligible imaginary part "
            f"{max_imag:.6e} > {atol:.6e}"
        )
    return np.real(normalized)


def cubic_kernel_convention_diagnostics(raw_kernel, normalized_kernel) -> dict[str, object]:
    """Return dtype and real/imaginary magnitude diagnostics."""

    raw = np.asarray(raw_kernel)
    normalized = np.asarray(normalized_kernel)
    return {
        "raw_dtype": str(raw.dtype),
        "normalized_dtype": str(normalized.dtype),
        "raw_max_real": float(np.max(np.abs(np.real(raw)))) if raw.size else 0.0,
        "raw_max_imag": float(np.max(np.abs(np.imag(raw)))) if raw.size else 0.0,
        "normalized_max_real": float(np.max(np.abs(np.real(normalized)))) if normalized.size else 0.0,
        "normalized_max_imag": float(np.max(np.abs(np.imag(normalized)))) if normalized.size else 0.0,
        "normalization": "KLM-normalized = (-1j) * raw physical cubic kernel",
    }
