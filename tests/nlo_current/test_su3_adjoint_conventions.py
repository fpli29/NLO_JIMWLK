from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    finite_diff_left_derivative,
    finite_diff_right_derivative,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


def test_adjoint_orthogonality() -> None:
    rng = np.random.default_rng(1234)
    gens = su3_generators_fundamental()
    eye = np.eye(8)
    for _ in range(5):
        S = adjoint_from_fundamental(random_su3(rng), gens)
        np.testing.assert_allclose(S @ S.T, eye, atol=2e-14, rtol=0.0)
        np.testing.assert_allclose(S.T @ S, eye, atol=2e-14, rtol=0.0)


def test_adjoint_structure_constant_invariance() -> None:
    rng = np.random.default_rng(2345)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    S = adjoint_from_fundamental(random_su3(rng), gens)

    lhs = np.einsum("abc,bd,ce->ade", f, S, S, optimize=True)
    rhs = np.einsum("af,fde->ade", S, f, optimize=True)
    rel = np.linalg.norm(lhs - rhs) / max(np.linalg.norm(rhs), 1.0)
    assert rel < 1e-12


def test_right_generator_to_left_generator_relation() -> None:
    rng = np.random.default_rng(3456)
    gens = su3_generators_fundamental()
    U = random_su3(rng)
    S = adjoint_from_fundamental(U, gens)
    probe = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))

    def scalar(Us):
        return float(np.real(np.trace(probe @ Us[0])))

    Us = [U]
    left = np.array(
        [finite_diff_left_derivative(scalar, Us, 0, b, eps=1e-6, gens=gens) for b in range(8)]
    )
    right = np.array(
        [finite_diff_right_derivative(scalar, Us, 0, a, eps=1e-6, gens=gens) for a in range(8)]
    )
    converted = np.einsum("ba,b->a", S, left, optimize=True)
    np.testing.assert_allclose(right, converted, atol=2e-9, rtol=2e-9)


def test_left_divergence_of_adjoint_column_is_zero() -> None:
    rng = np.random.default_rng(4567)
    gens = su3_generators_fundamental()
    U = random_su3(rng)
    Us = [U]
    values = []
    for b in range(8):
        div_b = 0.0
        for h in range(8):
            div_b += finite_diff_left_derivative(
                lambda Vs, h=h, b=b: adjoint_from_fundamental(Vs[0], gens)[h, b],
                Us,
                0,
                h,
                eps=1e-6,
                gens=gens,
            )
        values.append(div_b)

    values = np.array(values)
    report_dir = ROOT / "reports" / "nlo_current"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "left_divergence_identity_report.md").write_text(
        "# Left-Divergence Identity Check\n\n"
        "Convention tested: S_A^{ab}=2 Re Tr(t^a U t^b U^dagger), "
        "L^h F(U)=d/deps F(exp(i eps t^h)U).\n\n"
        f"seed: 4567\n\nmax_abs_sum_h_Lh_Shb: {np.max(np.abs(values)):.16e}\n\n"
        f"values_by_b: {np.array2string(values, precision=16)}\n",
        encoding="utf-8",
    )
    assert np.max(np.abs(values)) < 3e-10

