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
    kjjsj_A_LLR_from_kernel,
    kjjsj_B_LRR_from_kernel,
    kjjsj_V_virtual_from_kernel,
    synthetic_kjjsj_kernel,
)


def _coeff_setup(seed: int = 9505, nsite: int = 3):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    S_adj = np.stack([adjoint_from_fundamental(random_su3(rng), gens) for _ in range(nsite)])
    return rng, f, S_adj


def test_kjjsj_coefficient_shapes_and_realness() -> None:
    rng, f, S_adj = _coeff_setup()
    kernel = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")
    A = kjjsj_A_LLR_from_kernel(S_adj, kernel, f)
    B = kjjsj_B_LRR_from_kernel(S_adj, kernel, f)
    V = kjjsj_V_virtual_from_kernel(kernel, f)

    assert A.shape == (3, 3, 3, 8, 8, 8)
    assert B.shape == (3, 3, 3, 8, 8, 8)
    assert V.shape == (3, 3, 3, 8, 8, 8)
    assert np.isrealobj(A)
    assert np.isrealobj(B)
    assert np.isrealobj(V)


def test_kjjsj_xy_symmetry_candidates_for_llr_block() -> None:
    rng, f, S_adj = _coeff_setup(seed=9606)
    kernel_sym = synthetic_kjjsj_kernel(3, rng, xy_symmetry="symmetric")
    kernel_antisym = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")

    A_sym = kjjsj_A_LLR_from_kernel(S_adj, kernel_sym, f)
    A_antisym = kjjsj_A_LLR_from_kernel(S_adj, kernel_antisym, f)

    color_antisym_residual = np.max(np.abs(A_sym + np.swapaxes(A_sym, 3, 4)))
    combined_sym_residual = np.max(
        np.abs(A_antisym - np.transpose(A_antisym, (1, 0, 2, 4, 3, 5)))
    )

    report = ROOT / "reports" / "nlo_current" / "kjjsj_kernel_xy_symmetry_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# K_JJSJ Synthetic x,y Symmetry Check\n\n"
        f"symmetric_kernel_color_antisym_residual: {color_antisym_residual:.16e}\n\n"
        f"antisymmetric_kernel_combined_sym_residual: {combined_sym_residual:.16e}\n\n"
        "A symmetric synthetic kernel gives antisymmetry in color indices d,e at fixed x,y. "
        "An antisymmetric synthetic kernel gives symmetry under the combined exchange "
        "(x,d)<->(y,e), which is the convention used for the main cubic diagnostics.\n",
        encoding="utf-8",
    )

    assert color_antisym_residual < 1e-12
    assert combined_sym_residual < 1e-12


def test_kjjsj_virtual_one_third_factor() -> None:
    rng, f, _ = _coeff_setup(seed=9707)
    kernel = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")
    V = kjjsj_V_virtual_from_kernel(kernel, f)
    V_no_factor = np.einsum("wxyz,bde->xywdeb", kernel, f, optimize=True)
    ratio = np.linalg.norm(V_no_factor) / np.linalg.norm(V)
    np.testing.assert_allclose(ratio, 3.0, atol=1e-12, rtol=1e-12)

