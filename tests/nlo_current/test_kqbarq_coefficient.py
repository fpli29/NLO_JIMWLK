from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    su3_generators_fundamental,
)
from nlo_current.two_generator_terms import (  # noqa: E402
    kqbarq_A_from_kernel,
    kqbarq_C_left_from_A,
    qbarq_trace_block,
    sym_asym_parts,
)


def _fundamental_and_adjoint_lines(seed: int, n_site: int):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng) for _ in range(n_site)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    return rng, gens, U_fund, S_adj


def test_qbarq_zprime_equals_z_subtraction_matches_adjoint_convention() -> None:
    _, gens, U_fund, S_adj = _fundamental_and_adjoint_lines(seed=8101, n_site=5)
    max_abs_error = 0.0
    for z, U in enumerate(U_fund):
        residual = qbarq_trace_block(U, U, gens) - S_adj[z]
        max_abs_error = max(max_abs_error, float(np.max(np.abs(residual))))
    assert max_abs_error < 1e-10


def test_qbarq_coefficient_is_real_for_symmetric_real_kernel() -> None:
    rng, gens, U_fund, S_adj = _fundamental_and_adjoint_lines(seed=8202, n_site=3)
    raw = rng.normal(size=(3, 3, 3, 3))
    kernel = 0.5 * (raw + np.swapaxes(raw, 0, 1))
    kernel = 0.5 * (kernel + np.swapaxes(kernel, 2, 3))
    A = kqbarq_A_from_kernel(U_fund, S_adj, kernel, gens)
    assert float(np.max(np.abs(np.imag(A)))) < 1e-12


def test_qbarq_left_basis_conversion_shape_and_symmetry_split() -> None:
    rng, gens, U_fund, S_adj = _fundamental_and_adjoint_lines(seed=8303, n_site=3)
    raw = rng.normal(size=(3, 3, 3, 3))
    kernel = 0.5 * (raw + np.swapaxes(raw, 0, 1))
    kernel = 0.5 * (kernel + np.swapaxes(kernel, 2, 3))
    A = kqbarq_A_from_kernel(U_fund, S_adj, kernel, gens)
    C = kqbarq_C_left_from_A(A, S_adj)
    C_sym, C_asym = sym_asym_parts(C)

    assert A.shape == (3, 3, 8, 8)
    assert C.shape == (3, 8, 3, 8)
    assert C_sym.shape == C.shape
    assert C_asym.shape == C.shape
    np.testing.assert_allclose(C_sym + C_asym, C, atol=1e-12, rtol=1e-12)

