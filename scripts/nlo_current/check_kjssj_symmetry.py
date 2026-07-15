#!/usr/bin/env python3
"""Measure dense small-lattice K_JSSJ left-basis coefficient asymmetry."""

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
from nlo_current.two_generator_terms import kjssj_A_from_kernel, kjssj_C_left_from_A  # noqa: E402


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
    f = structure_constants(gens)
    S_adj = np.stack([adjoint_from_fundamental(random_su3(rng), gens) for _ in range(n_site)])
    kernel = synthetic_kernel(rng, n_site, sym_z_zprime=sym_z_zprime)
    A = kjssj_A_from_kernel(S_adj, kernel, f)
    C = kjssj_C_left_from_A(A, S_adj)
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
    }


def write_report(cases: list[dict[str, float | str | bool]]) -> Path:
    report = ROOT / "reports" / "nlo_current" / "kjssj_symmetry_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# K_JSSJ Coefficient Symmetry Report",
        "",
        "This is a dense small-lattice diagnostic, not a proof. It tests whether the",
        "left-basis coefficient C^{(x,a)(y,h)} built from the Appendix B ordered",
        "A_JSSJ block is symmetric under combined index exchange (x,a)<->(y,h).",
        "",
        "| seed | n_site | kernel symmetry | ||C|| | ||C-C^T|| | r_asym |",
        "|---:|---:|---|---:|---:|---:|",
    ]
    for case in cases:
        lines.append(
            f"| {case['seed']} | {case['n_site']} | {case['symmetry']} | "
            f"{case['norm_C']:.16e} | {case['norm_C_minus_CT']:.16e} | "
            f"{case['r_asym']:.16e} |"
        )

    max_ratio = max(float(case["r_asym"]) for case in cases)
    lines.extend(
        [
            "",
            "## Conclusion",
            "",
            f"max_r_asym: {max_ratio:.16e}",
            "",
        ]
    )
    if max_ratio < 1e-10:
        lines.append(
            "Under the tested synthetic kernel assumptions, C_JSSJ appears symmetric."
        )
    else:
        lines.append(
            "C_JSSJ is not symmetric in these tests. Do not simplify it to a pure "
            "symmetric score-current; keep the antisymmetric/commutator drift."
        )
    lines.append("")
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main() -> None:
    cases = [
        measure_case(seed=20260704, n_site=3, sym_z_zprime=False),
        measure_case(seed=20260705, n_site=3, sym_z_zprime=True),
    ]
    report = write_report(cases)
    for case in cases:
        print(
            f"seed={case['seed']} n_site={case['n_site']} "
            f"sym_z_zprime={case['sym_z_zprime']} r_asym={case['r_asym']:.16e}"
        )
    print(f"wrote {report}")


if __name__ == "__main__":
    main()

