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
from nlo_current.three_generator_terms import (  # noqa: E402
    kjjssj_A_LLR_from_kernel,
    kjjssj_B_LRR_from_kernel,
    kjjssj_V_virtual_from_kernel,
    synthetic_kjjssj_kernel,
)


def _coeff_setup(seed: int = 10101, nsite: int = 3):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    S_adj = np.stack([adjoint_from_fundamental(random_su3(rng), gens) for _ in range(nsite)])
    return rng, f, S_adj


def test_kjjssj_coefficient_shapes_and_realness() -> None:
    rng, f, S_adj = _coeff_setup()
    kernel = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)
    A = kjjssj_A_LLR_from_kernel(S_adj, kernel, f)
    B = kjjssj_B_LRR_from_kernel(S_adj, kernel, f)
    V = kjjssj_V_virtual_from_kernel(kernel, f)

    assert A.shape == (3, 3, 3, 8, 8, 8)
    assert B.shape == (3, 3, 3, 8, 8, 8)
    assert V.shape == (3, 3, 3, 8, 8, 8)
    assert np.isrealobj(A)
    assert np.isrealobj(B)
    assert np.isrealobj(V)


def test_kjjssj_klm_like_simultaneous_antisymmetry_and_induced_A_behavior() -> None:
    rng, f, S_adj = _coeff_setup(seed=10202)
    kernel = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)
    transformed = np.transpose(kernel, (0, 2, 1, 4, 3))
    kernel_residual = float(np.max(np.abs(kernel + transformed)))

    A = kjjssj_A_LLR_from_kernel(S_adj, kernel, f)
    combined_plus = np.max(np.abs(A + np.transpose(A, (1, 0, 2, 4, 3, 5))))
    combined_minus = np.max(np.abs(A - np.transpose(A, (1, 0, 2, 4, 3, 5))))

    report = ROOT / "reports" / "nlo_current" / "kjjssj_kernel_symmetry_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# K_JJSSJ KLM-Like Kernel Symmetry Check\n\n"
        f"kernel_simultaneous_antisymmetry_residual: {kernel_residual:.16e}\n\n"
        f"A_combined_plus_norm_max: {combined_plus:.16e}\n\n"
        f"A_combined_minus_norm_max: {combined_minus:.16e}\n\n"
        "The test enforces K(w;x,y;z,z') = -K(w;y,x;z',z). "
        "The induced A_LLR behavior under (x,d)<->(y,e) is measured rather "
        "than assumed.\n",
        encoding="utf-8",
    )

    assert kernel_residual < 1e-12


def test_kjjssj_virtual_one_third_factor() -> None:
    rng, f, _ = _coeff_setup(seed=10303)
    kernel = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)
    V = kjjssj_V_virtual_from_kernel(kernel, f)
    V_no_factor = np.einsum("wxyuv,acb->xywcba", kernel, f, optimize=True)
    ratio = np.linalg.norm(V_no_factor) / np.linalg.norm(V)
    np.testing.assert_allclose(ratio, 3.0, atol=1e-12, rtol=1e-12)


def test_kjjssj_zero_kernel_gives_zero_blocks() -> None:
    _, f, S_adj = _coeff_setup(seed=10404)
    kernel = np.zeros((3, 3, 3, 3, 3))
    np.testing.assert_allclose(kjjssj_A_LLR_from_kernel(S_adj, kernel, f), 0.0)
    np.testing.assert_allclose(kjjssj_B_LRR_from_kernel(S_adj, kernel, f), 0.0)
    np.testing.assert_allclose(kjjssj_V_virtual_from_kernel(kernel, f), 0.0)

