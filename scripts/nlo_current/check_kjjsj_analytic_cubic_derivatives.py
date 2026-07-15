from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.analytic_cubic_derivatives import (  # noqa: E402
    KJJSJ_BLOCKS,
    KJJSJSectorData,
    analytic_d2K3_KJJSJ,
    analytic_dK1_comm_KJJSJ,
    analytic_first_derivatives_KJJSJ,
    kjjsj_terms_from_blocks,
)
from nlo_current.coefficient_derivatives import compute_all_coefficient_derivatives_fd  # noqa: E402
from nlo_current.nlo_current_skeleton import NLOCurrentTerms  # noqa: E402
from nlo_current.nlo_velocity_evaluator import evaluate_velocity_from_terms  # noqa: E402
from nlo_current.physical_current_divergence import evaluate_current_divergence  # noqa: E402
from nlo_current.physical_density_operator import evaluate_direct_density_operator  # noqa: E402
from nlo_current.physical_coefficient_derivatives import (  # noqa: E402
    compute_physical_coefficient_derivatives,
)
from nlo_current.physical_kernels import KJSJIntegrationPolicy  # noqa: E402
from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.test_densities import evaluate_test_density  # noqa: E402


def _S_builder(gens):
    return lambda U: np.stack([adjoint_from_fundamental(V, gens) for V in U])


def _data_from_U(U, K, gens, f):
    return KJJSJSectorData(S_adj=_S_builder(gens)(U), KJJSJ=K, f=f)


def _terms_from_data(data: KJJSJSectorData, *, blocks=None) -> NLOCurrentTerms:
    K1, K2, K3 = kjjsj_terms_from_blocks(data, blocks=blocks)
    return NLOCurrentTerms(K1=K1, K2=K2, K3=K3, metadata={"sectors": {"KJJSJ": {}}})


def _residual(analytic, fd):
    analytic = np.asarray(analytic)
    fd = np.asarray(fd)
    diff = analytic - fd
    return {
        "max_abs": float(np.max(np.abs(diff))) if diff.size else 0.0,
        "rel": float(np.linalg.norm(diff) / (np.linalg.norm(fd) + 1.0e-30)),
        "real_max": float(np.max(np.abs(np.real(diff)))) if diff.size else 0.0,
        "imag_max": float(np.max(np.abs(np.imag(diff)))) if diff.size else 0.0,
        "analytic_norm": float(np.linalg.norm(analytic)),
        "fd_norm": float(np.linalg.norm(fd)),
        "expected_zero": bool(np.linalg.norm(fd) < 1.0e-24),
    }


def _compare_derivatives(analytic_first, analytic_d2, fd):
    return {
        "dK2_comm": _residual(analytic_first["dK2_comm"], fd["dK2"]),
        "LC_K3": _residual(analytic_first["LC_K3"], fd["dK3_first"]["LC_K3_ABC"]),
        "LB_K3": _residual(analytic_first["LB_K3"], fd["dK3_first"]["LB_K3_ABC"]),
        "d2K3": _residual(analytic_d2["total"], fd["d2K3"]),
    }


def _synthetic_setup(seed: int):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U = np.stack([random_su3(rng) for _ in range(2)])
    K = 0.03 * rng.normal(size=(2, 2, 2, 2))
    return gens, f, U, K


def _fd_for_blocks(U, K, gens, f, *, blocks=None, eps_first=2.0e-5, eps_second=5.0e-4):
    return compute_all_coefficient_derivatives_fd(
        lambda V, S: _terms_from_data(KJJSJSectorData(S, K, f), blocks=blocks).K2,
        lambda V, S: _terms_from_data(KJJSJSectorData(S, K, f), blocks=blocks).K3,
        U,
        _S_builder(gens),
        gens,
        eps_first=eps_first,
        eps_second=eps_second,
    )


def run_synthetic_scan():
    seeds = (20260756, 20260758)
    step_pairs = ((2.0e-5, 1.0e-3), (2.0e-5, 5.0e-4))
    rows = []
    for seed in seeds:
        gens, f, U, K = _synthetic_setup(seed)
        data = _data_from_U(U, K, gens, f)
        start = time.perf_counter()
        first = analytic_first_derivatives_KJJSJ(sector_data=data, return_diagnostics=True)
        d2 = analytic_d2K3_KJJSJ(sector_data=data, return_diagnostics=True)
        analytic_time = time.perf_counter() - start
        for eps_first, eps_second in step_pairs:
            start = time.perf_counter()
            fd = _fd_for_blocks(U, K, gens, f, eps_first=eps_first, eps_second=eps_second)
            fd_time = time.perf_counter() - start
            rows.append(
                {
                    "seed": seed,
                    "eps_first": eps_first,
                    "eps_second": eps_second,
                    "analytic_time": analytic_time,
                    "fd_time": fd_time,
                    "speedup": fd_time / max(analytic_time, 1.0e-30),
                    "residuals": _compare_derivatives(first, d2, fd),
                    "k1_status": analytic_dK1_comm_KJJSJ(sector_data=data),
                }
            )
    return rows


def run_block_scan():
    seed = 20260756
    gens, f, U, K = _synthetic_setup(seed)
    data = _data_from_U(U, K, gens, f)
    first = analytic_first_derivatives_KJJSJ(sector_data=data, return_diagnostics=True)
    d2 = analytic_d2K3_KJJSJ(sector_data=data, return_diagnostics=True)
    rows = []
    for block in KJJSJ_BLOCKS:
        start = time.perf_counter()
        fd = _fd_for_blocks(U, K, gens, f, blocks=(block,), eps_first=2.0e-5, eps_second=5.0e-4)
        fd_time = time.perf_counter() - start
        rows.append(
            {
                "block": block,
                "fd_time": fd_time,
                "dK2_comm": _residual(first["by_block"]["dK2_comm"][block], fd["dK2"]),
                "LC_K3": _residual(first["by_block"]["LC_K3"][block], fd["dK3_first"]["LC_K3_ABC"]),
                "LB_K3": _residual(first["by_block"]["LB_K3"][block], fd["dK3_first"]["LB_K3_ABC"]),
                "d2K3": _residual(d2["by_block"][block], fd["d2K3"]),
            }
        )
    return rows


def run_velocity_and_closure_check():
    gens, f, U, K = _synthetic_setup(20260756)
    data = _data_from_U(U, K, gens, f)
    terms = _terms_from_data(data)
    first = analytic_first_derivatives_KJJSJ(sector_data=data)
    d2 = analytic_d2K3_KJJSJ(sector_data=data)
    fd = _fd_for_blocks(U, K, gens, f, eps_first=2.0e-5, eps_second=5.0e-4)
    derivatives = {
        "dK2": first["dK2_comm"],
        "dK3_first": {"LC_K3_ABC": first["LC_K3"], "LB_K3_ABC": first["LB_K3"]},
        "d2K3": d2,
    }
    fd_derivatives = {
        "dK2": fd["dK2"],
        "dK3_first": fd["dK3_first"],
        "d2K3": fd["d2K3"],
    }
    density = evaluate_test_density(U, gens, "dipole_trace")
    analytic_velocity = evaluate_velocity_from_terms(
        terms, density.score, density.hessian_score, **derivatives
    )["velocity"]
    fd_velocity = evaluate_velocity_from_terms(
        terms, density.score, density.hessian_score, **fd_derivatives
    )["velocity"]

    density_builder = lambda V: evaluate_test_density(V, gens, "dipole_trace")
    terms_builder = lambda V: _terms_from_data(_data_from_U(V, K, gens, f))
    derivatives_builder = lambda V: {
        "dK2": analytic_first_derivatives_KJJSJ(sector_data=_data_from_U(V, K, gens, f))[
            "dK2_comm"
        ],
        "dK3_first": {
            "LC_K3_ABC": analytic_first_derivatives_KJJSJ(sector_data=_data_from_U(V, K, gens, f))[
                "LC_K3"
            ],
            "LB_K3_ABC": analytic_first_derivatives_KJJSJ(sector_data=_data_from_U(V, K, gens, f))[
                "LB_K3"
            ],
        },
        "d2K3": analytic_d2K3_KJJSJ(sector_data=_data_from_U(V, K, gens, f)),
    }
    direct = evaluate_direct_density_operator(
        U,
        density_builder,
        terms,
        gens=gens,
        fd_eps=1.0e-3,
        terms_builder=terms_builder,
        active_outer_indices=(0,),
    )
    current = evaluate_current_divergence(
        U,
        density_builder,
        density.score,
        density.hessian_score,
        terms,
        derivatives,
        gens=gens,
        fd_eps=1.0e-3,
        terms_builder=terms_builder,
        derivatives_builder=derivatives_builder,
        density_builder=density_builder,
        active_outer_indices=(0,),
    )
    no_hessian = evaluate_current_divergence(
        U,
        density_builder,
        density.score,
        density.hessian_score,
        terms,
        derivatives,
        gens=gens,
        fd_eps=1.0e-3,
        terms_builder=terms_builder,
        derivatives_builder=derivatives_builder,
        density_builder=density_builder,
        active_outer_indices=(0,),
        omit_hessian_score=True,
    )
    return {
        "velocity": _residual(analytic_velocity, fd_velocity),
        "direct": direct.value,
        "current": current.value,
        "closure_abs": float(abs(direct.value - current.value)),
        "closure_rel": float(abs(direct.value - current.value) / (abs(direct.value) + abs(current.value) + 1e-30)),
        "omit_hessian_abs": float(abs(direct.value - no_hessian.value)),
    }


def run_physical_zero_smoke():
    rng = np.random.default_rng(20260759)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U = np.stack([random_su3(rng) for _ in range(2)])
    coords = np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)
    policy = KJSJIntegrationPolicy(
        quadrature_weights=np.ones(2) / 2,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
        description="KJJSJ analytic cubic expected-zero physical smoke policy",
    )
    config = PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=2.0e-5,
        fd_eps_second=5.0e-4,
    )
    analytic = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="analytic",
        sector_filter=("KJJSJ",),
    )
    fd = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="finite_difference",
        sector_filter=("KJJSJ",),
    )
    return {
        "dK2_comm": _residual(analytic.dK2, fd.dK2),
        "LC_K3": _residual(analytic.LC_K3, fd.LC_K3),
        "LB_K3": _residual(analytic.LB_K3, fd.LB_K3),
        "d2K3": _residual(analytic.d2K3, fd.d2K3),
        "metadata": analytic.metadata,
    }


def _fmt_res(res):
    return (
        f"{res['max_abs']:.3e} | {res['rel']:.3e} | {res['real_max']:.3e} | "
        f"{res['imag_max']:.3e} | {res['analytic_norm']:.3e} | {res['fd_norm']:.3e} | "
        f"{res['expected_zero']}"
    )


def write_report(synthetic_rows, block_rows, velocity_closure, physical_zero):
    report = ROOT / "reports" / "nlo_current" / "kjjsj_analytic_cubic_validation_report.md"
    lines = [
        "# KJJSJ Analytic Cubic Derivative Validation Report",
        "",
        "## Scope",
        "",
        "This report validates only the non-production analytic coefficient derivatives",
        "for the `KJJSJ` cubic sector. `KJJSSJ` remains explicitly pending. The",
        "finite-difference backend is the reference oracle.",
        "",
        "## Implemented Formulas",
        "",
        "- `LC_K3`: `(LC_K3)^{AB} = L_C K3^{ABC}` from the same canonicalized KJJSJ normal-form tensors.",
        "- `LB_K3`: `(LB_K3)^{AC} = L_B K3^{ABC}` with ordered combined-index contractions.",
        "- `d2K3`: `d2K3^A = L_B L_C K3^{ABC}` using ordered second adjoint derivatives.",
        "- `dK2_comm`: derivative of the KJJSJ quadratic commutator correction produced by canonicalization.",
        "- `dK1_comm`: classified from canonical linear terms; not assumed zero.",
        "",
        "The physical adapter supplies KLM-normalized cubic coefficients. The analytic",
        "derivative code does not apply the `(-i)` normalization again.",
        "",
        "## Full Synthetic Nonzero FD Step Scan",
        "",
        "| seed | eps first | eps second | quantity | max abs | relative | real max | imag max | analytic norm | FD norm | expected zero | analytic s | FD s | speedup |",
        "|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|",
    ]
    for row in synthetic_rows:
        for quantity, residual in row["residuals"].items():
            lines.append(
                f"| {row['seed']} | {row['eps_first']:.0e} | {row['eps_second']:.0e} | `{quantity}` | "
                f"{_fmt_res(residual)} | {row['analytic_time']:.3e} | {row['fd_time']:.3e} | {row['speedup']:.2f} |"
            )
    lines.extend(
        [
            "",
            "## Per-Block Residuals",
            "",
            "The block scan uses seed `20260756`, `eps_first=2e-5`, and",
            "`eps_second=5e-4`.",
            "",
            "| block | quantity | max abs | relative | real max | imag max | analytic norm | FD norm | expected zero |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in block_rows:
        for quantity in ("dK2_comm", "LC_K3", "LB_K3", "d2K3"):
            lines.append(f"| `{row['block']}` | `{quantity}` | {_fmt_res(row[quantity])} |")
    lines.extend(
        [
            "",
            "## Velocity and Closure",
            "",
            "| diagnostic | value |",
            "|---|---:|",
            f"| velocity max residual vs FD | {velocity_closure['velocity']['max_abs']:.3e} |",
            f"| velocity relative residual vs FD | {velocity_closure['velocity']['rel']:.3e} |",
            f"| direct density operator | `{velocity_closure['direct']}` |",
            f"| current divergence | `{velocity_closure['current']}` |",
            f"| closure absolute residual | {velocity_closure['closure_abs']:.3e} |",
            f"| closure relative residual | {velocity_closure['closure_rel']:.3e} |",
            f"| omit-Hessian absolute residual | {velocity_closure['omit_hessian_abs']:.3e} |",
            "",
            "## Physical Two-Site Expected-Zero Smoke Check",
            "",
            "For the two-site physical coordinate setup, the diagnostic unbarred physical",
            "`KJJSJ` array is structurally zero. This smoke check exercises the physical",
            "adapter/backend path but does not replace the synthetic nonzero oracle scan.",
            "",
            "| quantity | max abs | relative | real max | imag max | analytic norm | FD norm | expected zero |",
            "|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for quantity, residual in physical_zero.items():
        if quantity == "metadata":
            continue
        lines.append(f"| `{quantity}` | {_fmt_res(residual)} |")
    lines.extend(
        [
            "",
            "## Backend Status",
            "",
            "- `backend=\"analytic\"` uses no global FD fallback for `KJJSJ`.",
            "- `KJJSSJ` remains pending and must raise under `backend=\"analytic\"`.",
            "- No complex-to-real cast is required by the KJJSJ analytic path; tested",
            "  KLM-normalized derivative arrays are real to numerical precision.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def main():
    synthetic_rows = run_synthetic_scan()
    block_rows = run_block_scan()
    velocity_closure = run_velocity_and_closure_check()
    physical_zero = run_physical_zero_smoke()
    write_report(synthetic_rows, block_rows, velocity_closure, physical_zero)
    best = min(row["residuals"]["d2K3"]["max_abs"] for row in synthetic_rows)
    print(f"KJJSJ analytic cubic validation report written; best d2K3 max residual {best:.3e}")


if __name__ == "__main__":
    main()
