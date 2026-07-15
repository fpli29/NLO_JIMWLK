from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


def _b5_lhs_rhs(Sz: np.ndarray, Szp: np.ndarray, f: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nc = 3.0
    lhs = np.einsum("abc,def,be,cf->ad", f, f, Sz, Szp, optimize=True) - nc * Sz
    rhs = np.einsum("apc,def,pe,cf->ad", f, f, Sz, Szp - Sz, optimize=True)
    return lhs, rhs


def test_klm_appendix_b_identity_random_adjoint_lines() -> None:
    rng = np.random.default_rng(5678)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    for _ in range(5):
        Sz = adjoint_from_fundamental(random_su3(rng), gens)
        Szp = adjoint_from_fundamental(random_su3(rng), gens)
        lhs, rhs = _b5_lhs_rhs(Sz, Szp, f)
        rel = np.linalg.norm(lhs - rhs) / max(np.linalg.norm(lhs), np.linalg.norm(rhs), 1.0)
        assert rel < 1e-10


def test_klm_appendix_b_subtraction_vanishes_when_zprime_equals_z() -> None:
    rng = np.random.default_rng(6789)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    Sz = adjoint_from_fundamental(random_su3(rng), gens)
    lhs, rhs = _b5_lhs_rhs(Sz, Sz, f)
    np.testing.assert_allclose(lhs, np.zeros_like(lhs), atol=2e-12, rtol=0.0)
    np.testing.assert_allclose(rhs, np.zeros_like(rhs), atol=2e-12, rtol=0.0)

