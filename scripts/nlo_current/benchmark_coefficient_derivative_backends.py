from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.physical_coefficient_derivatives import compute_physical_coefficient_derivatives  # noqa: E402
from nlo_current.physical_kernels import KJSJIntegrationPolicy  # noqa: E402
from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig  # noqa: E402
from nlo_current.su3_adjoint import random_su3, structure_constants, su3_generators_fundamental  # noqa: E402


def _policy(nsite):
    return KJSJIntegrationPolicy(
        quadrature_weights=np.ones(nsite) / nsite,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
        description="analytic derivative benchmark policy",
    )


def _config():
    return PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=2.0e-5,
        fd_eps_second=5.0e-4,
    )


def run_benchmark():
    rng = np.random.default_rng(20260755)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)
    U = np.stack([random_su3(rng) for _ in range(2)])
    sectors = ("KJSJ", "KJSSJ", "Kqbarq")
    rows = []
    for sector in sectors:
        timings = {}
        outputs = {}
        for backend in ("analytic", "finite_difference"):
            start = time.perf_counter()
            outputs[backend] = compute_physical_coefficient_derivatives(
                U,
                coords,
                gens,
                f,
                integration_policy=_policy(2),
                config=_config(),
                backend=backend,
                sector_filter=(sector,),
            )
            timings[backend] = time.perf_counter() - start
        diff = outputs["analytic"].dK2 - outputs["finite_difference"].dK2
        rows.append(
            {
                "sector": sector,
                "analytic_time": timings["analytic"],
                "fd_time": timings["finite_difference"],
                "speedup": timings["finite_difference"] / max(timings["analytic"], 1.0e-30),
                "dK2_max_abs": float(np.max(np.abs(diff))),
                "dK2_rel": float(
                    np.linalg.norm(diff) / (np.linalg.norm(outputs["finite_difference"].dK2) + 1.0e-30)
                ),
            }
        )
    cubic_rows = []
    cubic_outputs = {}
    cubic_timings = {}
    for backend in ("analytic", "finite_difference", "hybrid_local_fd"):
        start = time.perf_counter()
        cubic_outputs[backend] = compute_physical_coefficient_derivatives(
            U,
            coords,
            gens,
            f,
            integration_policy=_policy(2),
            config=_config(),
            backend=backend,
            sector_filter=("KJJSJ",),
        )
        cubic_timings[backend] = time.perf_counter() - start
    fd = cubic_outputs["finite_difference"]
    for backend in ("analytic", "hybrid_local_fd"):
        out = cubic_outputs[backend]
        cubic_rows.append(
            {
                "backend": backend,
                "time": cubic_timings[backend],
                "fd_time": cubic_timings["finite_difference"],
                "speedup": cubic_timings["finite_difference"] / max(cubic_timings[backend], 1.0e-30),
                "dK2_max_abs": float(np.max(np.abs(out.dK2 - fd.dK2))),
                "LC_max_abs": float(np.max(np.abs(out.LC_K3 - fd.LC_K3))),
                "LB_max_abs": float(np.max(np.abs(out.LB_K3 - fd.LB_K3))),
                "d2_max_abs": float(np.max(np.abs(out.d2K3 - fd.d2K3))),
                "fallback_used": bool(out.metadata.get("fallback_used", False)),
            }
        )
    return rows, cubic_rows


def write_report(rows, cubic_rows):
    report = ROOT / "reports" / "nlo_current" / "analytic_coefficient_derivative_benchmark.md"
    lines = [
        "# Analytic Coefficient Derivative Benchmark",
        "",
        "This benchmark compares implemented analytic derivatives with the",
        "finite-difference oracle on the smallest physical diagnostic setup. It",
        "does not claim asymptotic production scaling. The physical two-site",
        "`KJJSJ` kernel is expected-zero; nonzero KJJSJ residuals are reported in",
        "`reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`.",
        "",
        "| sector | analytic s | FD s | observed speedup | dK2 max residual | dK2 relative residual |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['sector']}` | {row['analytic_time']:.6e} | {row['fd_time']:.6e} | "
            f"{row['speedup']:.2f} | {row['dK2_max_abs']:.3e} | {row['dK2_rel']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## KJJSJ Physical Two-Site Backend Benchmark",
            "",
            "| backend | time s | FD s | observed speedup | dK2 max | LC max | LB max | d2 max | fallback used |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in cubic_rows:
        lines.append(
            f"| `{row['backend']}` | {row['time']:.6e} | {row['fd_time']:.6e} | "
            f"{row['speedup']:.2f} | {row['dK2_max_abs']:.3e} | {row['LC_max_abs']:.3e} | "
            f"{row['LB_max_abs']:.3e} | {row['d2_max_abs']:.3e} | {row['fallback_used']} |"
        )
    lines.append("")
    lines.append("`KJJSJ` is analytic-complete for the validated diagnostic backend; `KJJSSJ` remains pending.")
    report.write_text("\n".join(lines), encoding="utf-8")


def main():
    rows, cubic_rows = run_benchmark()
    write_report(rows, cubic_rows)
    for row in rows:
        print(
            f"{row['sector']}: analytic={row['analytic_time']:.4e}s "
            f"fd={row['fd_time']:.4e}s speedup={row['speedup']:.2f} "
            f"rel={row['dK2_rel']:.3e}"
        )
    for row in cubic_rows:
        print(
            f"KJJSJ {row['backend']}: time={row['time']:.4e}s "
            f"fd={row['fd_time']:.4e}s speedup={row['speedup']:.2f} "
            f"d2_max={row['d2_max_abs']:.3e} fallback={row['fallback_used']}"
        )


if __name__ == "__main__":
    main()
