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
        description="analytic derivative oracle comparison policy",
    )


def _config(fd_eps_first):
    return PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=fd_eps_first,
        fd_eps_second=5.0e-4,
    )


def _residual(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    diff = a - b
    return {
        "max_abs": float(np.max(np.abs(diff))) if diff.size else 0.0,
        "rel": float(np.linalg.norm(diff) / (np.linalg.norm(b) + 1.0e-30)),
        "imag": float(max(np.max(np.abs(np.imag(a))), np.max(np.abs(np.imag(b))))) if diff.size else 0.0,
        "expected_zero": bool(np.linalg.norm(b) < 1.0e-24),
    }


def run_check():
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)
    sectors = ("KJSJ", "KJSSJ", "Kqbarq")
    fd_steps = (2.0e-5, 1.0e-5)
    seeds = (20260753, 20260754)
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        U = np.stack([random_su3(rng) for _ in range(2)])
        for fd_eps in fd_steps:
            config = _config(fd_eps)
            for sector in sectors:
                start = time.perf_counter()
                analytic = compute_physical_coefficient_derivatives(
                    U,
                    coords,
                    gens,
                    f,
                    integration_policy=_policy(2),
                    config=config,
                    backend="analytic",
                    sector_filter=(sector,),
                )
                analytic_time = time.perf_counter() - start
                start = time.perf_counter()
                fd = compute_physical_coefficient_derivatives(
                    U,
                    coords,
                    gens,
                    f,
                    integration_policy=_policy(2),
                    config=config,
                    backend="finite_difference",
                    sector_filter=(sector,),
                )
                fd_time = time.perf_counter() - start
                rows.append(
                    {
                        "seed": seed,
                        "fd_eps": fd_eps,
                        "sector": sector,
                        "dK2": _residual(analytic.dK2, fd.dK2),
                        "LC": _residual(analytic.LC_K3, fd.LC_K3),
                        "LB": _residual(analytic.LB_K3, fd.LB_K3),
                        "d2": _residual(analytic.d2K3, fd.d2K3),
                        "analytic_time": analytic_time,
                        "fd_time": fd_time,
                        "speedup": fd_time / max(analytic_time, 1.0e-30),
                    }
                )
    return rows


def write_report(rows):
    report = ROOT / "reports" / "nlo_current" / "analytic_coefficient_derivative_validation_report.md"
    lines = [
        "# Analytic Coefficient Derivative Validation Report",
        "",
        "## Scope",
        "",
        "This report compares implemented analytic/local physical coefficient derivatives",
        "against the preserved finite-difference oracle. It is non-production and does",
        "not claim regulator independence or physical positivity.",
        "",
        "## Implemented Analytic Sectors",
        "",
        "- `KJSJ`: analytic `dK2`.",
        "- `KJSSJ`: analytic `dK2`.",
        "- `Kqbarq`: analytic trace, subtraction, and full `dK2`.",
        "- `KJJSJ`: analytic `dK2_comm`, `LC_K3`, `LB_K3`, and ordered `d2K3`.",
        "",
        "## Pending Sectors",
        "",
        "- `KJJSSJ`: cubic `LC_K3`, `LB_K3`, and `d2K3` pending.",
        "",
        "`KJJSJ` is validated by the dedicated nonzero oracle report",
        "`reports/nlo_current/kjjsj_analytic_cubic_validation_report.md`.",
        "`KJJSSJ` is not marked analytic-complete; use `finite_difference` or",
        "`hybrid_local_fd` explicitly for that sector.",
        "",
        "## Residual Table",
        "",
        "| seed | fd eps | sector | dK2 max | dK2 rel | LC max | LB max | d2 max | analytic s | FD s | speedup |",
        "|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {seed} | {fd_eps:.0e} | `{sector}` | {dK2_max:.3e} | {dK2_rel:.3e} | "
            "{LC_max:.3e} | {LB_max:.3e} | {d2_max:.3e} | {analytic_time:.3e} | "
            "{fd_time:.3e} | {speedup:.2f} |".format(
                seed=row["seed"],
                fd_eps=row["fd_eps"],
                sector=row["sector"],
                dK2_max=row["dK2"]["max_abs"],
                dK2_rel=row["dK2"]["rel"],
                LC_max=row["LC"]["max_abs"],
                LB_max=row["LB"]["max_abs"],
                d2_max=row["d2"]["max_abs"],
                analytic_time=row["analytic_time"],
                fd_time=row["fd_time"],
                speedup=row["speedup"],
            )
        )
    lines.extend(
        [
            "",
            "## Dtype and Complex Checks",
            "",
            "The analytic two-generator paths preserve complex intermediate values. The",
            "tested physical outputs are real after `np.real_if_close`; no `ComplexWarning`",
            "was emitted by the validation tests.",
            "",
            "## Closure Status",
            "",
            "`tests/nlo_current/test_analytic_coefficient_derivatives.py` checks analytic",
            "two-generator velocity agreement with the FD reference and projected density",
            "closure for `KJSSJ` plus `Kqbarq`.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def main():
    rows = run_check()
    write_report(rows)
    for row in rows:
        print(
            f"{row['sector']} seed={row['seed']} fd={row['fd_eps']:.0e} "
            f"dK2_max={row['dK2']['max_abs']:.3e} dK2_rel={row['dK2']['rel']:.3e} "
            f"speedup={row['speedup']:.2f}"
        )


if __name__ == "__main__":
    main()
