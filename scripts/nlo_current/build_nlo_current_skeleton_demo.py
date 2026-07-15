from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.coefficient_derivatives import (  # noqa: E402
    compute_all_coefficient_derivatives_fd,
    velocity_from_coeff_derivative_backend,
)
from nlo_current.nlo_current_skeleton import NLOCurrentTerms, assemble_nlo_current_terms  # noqa: E402
from nlo_current.nlo_velocity_evaluator import (  # noqa: E402
    evaluate_velocity_from_terms,
    evaluate_velocity_score_only,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.synthetic_kernels import synthetic_kernels_all  # noqa: E402


SEED = 20260704
NSITE = 2
DERIVATIVE_DEMO_NSITE = 1
DERIVATIVE_EPS_FIRST = 2e-5
DERIVATIVE_EPS_SECOND = 7e-4


def _markdown_norm_table(sectors: dict) -> str:
    lines = ["| sector | K1 norm | K2 norm | K3 norm | notes |", "|---|---:|---:|---:|---|"]
    for name in ("KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ"):
        info = sectors[name]
        norms = info["norms"]
        notes = []
        if info.get("missing"):
            notes.append("missing")
        if "include_commutators" in info:
            notes.append(f"commutators={info['include_commutators']}")
            notes.append(f"K2_comm_terms={info.get('quadratic_comm_terms', 0)}")
            notes.append(f"K1_comm_terms={info.get('linear_comm_terms', 0)}")
        if info.get("ordered_generic_second_order"):
            notes.append("ordered K2")
        lines.append(
            "| {name} | {K1:.16e} | {K2:.16e} | {K3:.16e} | {notes} |".format(
                name=name,
                K1=norms["K1"],
                K2=norms["K2"],
                K3=norms["K3"],
                notes=", ".join(notes) if notes else "none",
            )
        )
    return "\n".join(lines)


def _S_builder(gens):
    return lambda U_fund: np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])


def _coefficient_derivative_demo():
    rng = np.random.default_rng(SEED + 1)
    gens = su3_generators_fundamental()
    U_fund = np.stack([random_su3(rng) for _ in range(DERIVATIVE_DEMO_NSITE)])
    S_builder = _S_builder(gens)
    dim = DERIVATIVE_DEMO_NSITE * gens.shape[0]
    Q_coeff = rng.normal(size=(3, 3)) + 1.0j * rng.normal(size=(3, 3))
    M2 = 0.2 * rng.normal(size=(dim, dim))
    M3 = 0.04 * rng.normal(size=(dim, dim, dim))

    def phi(U):
        return float(np.real(np.trace(Q_coeff @ U[0])))

    def K2_callback(U, _S):
        return phi(U) * M2

    def K3_callback(U, _S):
        return phi(U) * M3

    derivatives = compute_all_coefficient_derivatives_fd(
        K2_callback,
        K3_callback,
        U_fund,
        S_builder,
        gens,
        eps_first=DERIVATIVE_EPS_FIRST,
        eps_second=DERIVATIVE_EPS_SECOND,
    )
    terms = NLOCurrentTerms(
        K1=np.zeros(dim),
        K2=K2_callback(U_fund, S_builder(U_fund)),
        K3=K3_callback(U_fund, S_builder(U_fund)),
        metadata={},
    )
    score = rng.normal(size=dim)
    raw_hessian = rng.normal(size=(dim, dim))
    hessian_score = 0.5 * (raw_hessian + raw_hessian.T)
    with_derivatives = velocity_from_coeff_derivative_backend(
        terms,
        score,
        hessian_score,
        derivatives,
    )
    without_derivatives = evaluate_velocity_from_terms(terms, score, hessian_score)
    return {
        "dim": dim,
        "dK2_norm": float(np.linalg.norm(derivatives["dK2"])),
        "LC_K3_norm": float(np.linalg.norm(derivatives["dK3_first"]["LC_K3_ABC"])),
        "LB_K3_norm": float(np.linalg.norm(derivatives["dK3_first"]["LB_K3_ABC"])),
        "d2K3_norm": float(np.linalg.norm(derivatives["d2K3"])),
        "velocity_with_norm": float(np.linalg.norm(with_derivatives["velocity"])),
        "velocity_without_norm": float(np.linalg.norm(without_derivatives["velocity"])),
        "velocity_difference_norm": float(
            np.linalg.norm(with_derivatives["velocity"] - without_derivatives["velocity"])
        ),
        "warnings_with": len(with_derivatives["diagnostics"]["warnings"]),
        "warnings_without": len(without_derivatives["diagnostics"]["warnings"]),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(NSITE)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    kernels = synthetic_kernels_all(NSITE, rng)

    terms = assemble_nlo_current_terms(
        U_fund,
        S_adj,
        kernels,
        gens,
        f,
        include_commutators=True,
    )
    terms.validate_shapes()

    score = rng.normal(size=terms.dim)
    raw_hessian = rng.normal(size=(terms.dim, terms.dim))
    hessian_score = 0.5 * (raw_hessian + raw_hessian.T)
    velocity_result = evaluate_velocity_from_terms(terms, score, hessian_score)
    score_only_velocity = evaluate_velocity_score_only(terms, score)
    derivative_demo = _coefficient_derivative_demo()

    total_norms = terms.norms()
    k2_asymmetry = float(np.linalg.norm(terms.K2 - terms.K2.T))
    report = ROOT / "reports" / "nlo_current" / "nlo_current_skeleton_demo_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# NLO Current Skeleton Demo Report\n\n"
        "This report was generated by `scripts/nlo_current/build_nlo_current_skeleton_demo.py`.\n\n"
        "## Run Configuration\n\n"
        f"- random_seed: {SEED}\n"
        f"- nsite: {NSITE}\n"
        f"- combined_dimension_D: {terms.dim}\n"
        f"- commutators_included: {terms.metadata['commutators']['included']}\n"
        "- kernel_source: synthetic tiny-lattice kernels only\n"
        "- production_status: non-production diagnostic skeleton only\n\n"
        "## Sector Norms\n\n"
        f"{_markdown_norm_table(terms.metadata['sectors'])}\n\n"
        "## Total Norms\n\n"
        f"- K1_norm: {total_norms['K1']:.16e}\n"
        f"- K2_norm: {total_norms['K2']:.16e}\n"
        f"- K3_norm: {total_norms['K3']:.16e}\n"
        f"- ordered_K2_asymmetry_norm: {k2_asymmetry:.16e}\n\n"
        "## Velocity Diagnostics\n\n"
        f"- velocity_norm: {velocity_result['diagnostics']['velocity_norm']:.16e}\n"
        f"- score_only_velocity_norm: {np.linalg.norm(score_only_velocity):.16e}\n"
        "- coefficient_derivatives_omitted: "
        f"{velocity_result['diagnostics']['coefficient_derivatives_omitted']}\n"
        "- derivative_warnings:\n"
        + "".join(f"  - {warning}\n" for warning in velocity_result["diagnostics"]["warnings"])
        + "\n## Coefficient-Derivative Tiny Demo\n\n"
        f"- derivative_demo_nsite: {DERIVATIVE_DEMO_NSITE}\n"
        f"- derivative_demo_D: {derivative_demo['dim']}\n"
        f"- eps_first: {DERIVATIVE_EPS_FIRST:.16e}\n"
        f"- eps_second: {DERIVATIVE_EPS_SECOND:.16e}\n"
        f"- dK2_norm: {derivative_demo['dK2_norm']:.16e}\n"
        f"- LC_K3_norm: {derivative_demo['LC_K3_norm']:.16e}\n"
        f"- LB_K3_norm: {derivative_demo['LB_K3_norm']:.16e}\n"
        f"- d2K3_norm: {derivative_demo['d2K3_norm']:.16e}\n"
        f"- derivative_enabled_velocity_norm: {derivative_demo['velocity_with_norm']:.16e}\n"
        f"- derivative_omitted_velocity_norm: {derivative_demo['velocity_without_norm']:.16e}\n"
        f"- velocity_difference_norm: {derivative_demo['velocity_difference_norm']:.16e}\n"
        f"- derivative_enabled_warning_count: {derivative_demo['warnings_with']}\n"
        f"- derivative_omitted_warning_count: {derivative_demo['warnings_without']}\n"
        "- finite_difference_note: feasible only for tiny diagnostic lattices\n"
        + "\n## Explicit Warnings\n\n"
        "- Physical NLO kernels were not used.\n"
        "- Coefficient derivative arrays were omitted in the main skeleton smoke path.\n"
        "- The derivative-enabled section uses synthetic callbacks on a tiny dense lattice.\n"
        "- Score and Hessian-score inputs were synthetic random arrays.\n"
        "- This is not a production NLO evolution implementation.\n",
        encoding="utf-8",
    )

    print(f"wrote {report}")
    print(f"D={terms.dim}")
    print(f"K norms={total_norms}")
    print(f"velocity_norm={velocity_result['diagnostics']['velocity_norm']:.16e}")
    print(f"coefficient_derivative_warnings={len(velocity_result['diagnostics']['warnings'])}")
    print(f"derivative_demo_velocity_difference={derivative_demo['velocity_difference_norm']:.16e}")


if __name__ == "__main__":
    main()
