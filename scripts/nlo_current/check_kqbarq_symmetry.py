#!/usr/bin/env python3
"""Measure dense small-lattice K_qbarq left-basis coefficient asymmetry."""

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
from nlo_current.two_generator_terms import kqbarq_A_from_kernel, kqbarq_C_left_from_A  # noqa: E402


def synthetic_kernel(rng: np.random.Generator, n_site: int, sym_z_zprime: bool) -> np.ndarray:
    raw = rng.normal(size=(n_site, n_site, n_site, n_site))
    kernel = 0.5 * (raw + np.swapaxes(raw, 0, 1))
    if sym_z_zprime:
        kernel = 0.5 * (kernel + np.swapaxes(kernel, 2, 3))
    norm = np.linalg.norm(kernel)
    return kernel / norm


def measure_case(seed: int, n_site: int, sym_z_zprime: bool) -> dict[str, float | str | bool]:
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng) for _ in range(n_site)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    kernel = synthetic_kernel(rng, n_site, sym_z_zprime=sym_z_zprime)

    A = kqbarq_A_from_kernel(U_fund, S_adj, kernel, gens)
    C = kqbarq_C_left_from_A(A, S_adj)
    C_T = np.transpose(C, (2, 3, 0, 1))
    asym_norm = float(np.linalg.norm(C - C_T))
    total_norm = float(np.linalg.norm(C))
    ratio = asym_norm / max(total_norm, 1e-300)
    return {
        "seed": seed,
        "n_site": n_site,
        "symmetry": "K(x,y;z,z')=K(y,x;z,z')"
        + (" and K(x,y;z,z')=K(x,y;z',z)" if sym_z_zprime else ""),
        "sym_z_zprime": sym_z_zprime,
        "norm_C": total_norm,
        "norm_C_minus_CT": asym_norm,
        "r_asym": ratio,
        "max_abs_imag_A": float(np.max(np.abs(np.imag(A)))),
        "max_abs_imag_C": float(np.max(np.abs(np.imag(C)))),
    }


def write_report(cases: list[dict[str, float | str | bool]]) -> Path:
    report = ROOT / "reports" / "nlo_current" / "kqbarq_symmetry_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# K_qbarq Coefficient Symmetry Report",
        "",
        "This dense small-lattice diagnostic tests whether the left-basis",
        "coefficient C^{(x,a)(y,h)} built from the ordered K_qbarq A block is",
        "symmetric under combined index exchange (x,a)<->(y,h). It is not a proof.",
        "",
        "| seed | n_site | kernel symmetry | ||C|| | ||C-C^T|| | r_asym | max Im A | max Im C |",
        "|---:|---:|---|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        lines.append(
            f"| {case['seed']} | {case['n_site']} | {case['symmetry']} | "
            f"{case['norm_C']:.16e} | {case['norm_C_minus_CT']:.16e} | "
            f"{case['r_asym']:.16e} | {case['max_abs_imag_A']:.16e} | "
            f"{case['max_abs_imag_C']:.16e} |"
        )

    max_ratio = max(float(case["r_asym"]) for case in cases)
    lines.extend(["", "## Conclusion", "", f"max_r_asym: {max_ratio:.16e}", ""])
    if max_ratio > 1e-2:
        lines.append(
            "C_qbarq has an order-one antisymmetric component in these tests. "
            "Keep the generic ordered-current / commutator-drift representation."
        )
    elif max_ratio < 1e-10:
        lines.append(
            "C_qbarq appears symmetric in this synthetic test, but this remains "
            "a diagnostic and needs analytic confirmation before simplification."
        )
    else:
        lines.append(
            "C_qbarq has a visible nonzero antisymmetric component in this diagnostic. "
            "Do not simplify it without analytic control."
        )
    lines.append("")

    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main() -> None:
    cases = [
        measure_case(seed=20260714, n_site=3, sym_z_zprime=False),
        measure_case(seed=20260715, n_site=3, sym_z_zprime=True),
    ]
    report = write_report(cases)
    for case in cases:
        print(
            f"seed={case['seed']} n_site={case['n_site']} "
            f"sym_z_zprime={case['sym_z_zprime']} r_asym={case['r_asym']:.16e} "
            f"max_imag_C={case['max_abs_imag_C']:.16e}"
        )
    print(f"wrote {report}")


if __name__ == "__main__":
    main()

