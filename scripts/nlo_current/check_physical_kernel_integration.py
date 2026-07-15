from __future__ import annotations

from pathlib import Path
import sys
import warnings

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.nlo_current_skeleton import assemble_nlo_current_terms  # noqa: E402
from nlo_current.physical_kernel_adapter import (  # noqa: E402
    finite_kernel_stats,
    physical_kernels_for_skeleton,
)
from nlo_current.physical_kernels import (  # noqa: E402
    KJSJIntegrationPolicy,
    KJJSSJ_unbarred_value,
    KJJSJ_unbarred_value,
    KJSSJ_unbarred_value,
    KJSJ_unbarred_value,
    Kqbarq_unbarred_value,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)


REPORT = ROOT / "reports" / "nlo_current" / "physical_kernel_integration_report.md"
SEED = 20260711
PARAMS = {
    "Nc": 3,
    "nf": 2,
    "alpha_s": 0.3,
    "singularity_policy": "eps",
    "eps": 1.0e-6,
}
COORDS = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.2],
        [0.3, 1.1],
        [1.4, 1.3],
    ],
    dtype=float,
)
POINT_COORDS = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.2],
        [0.3, 1.1],
        [1.4, 1.3],
        [2.1, 0.7],
    ],
    dtype=float,
)


def _policy(ncoords: int, *, mu: float = 1.3) -> KJSJIntegrationPolicy:
    return KJSJIntegrationPolicy(
        quadrature_weights=np.ones(ncoords) / ncoords,
        mu=mu,
        exclude_coincident_labels=("x", "y", "z"),
        principal_value="none",
        subtraction="diagonal_zero",
        finite_volume_boundary="finite_coordinate_sum",
        description="diagnostic equal-weight finite coordinate sum",
    )


def _abs(value) -> float:
    return float(abs(value))


def symmetry_residuals() -> dict[str, float]:
    p = {k: v for k, v in PARAMS.items() if k != "singularity_policy" and k != "eps"}
    kjsj_policy = _policy(POINT_COORDS.shape[0])
    return {
        "KJSJ_xy": _abs(
            KJSJ_unbarred_value(
                POINT_COORDS,
                0,
                1,
                2,
                singularity_policy=PARAMS["singularity_policy"],
                eps=PARAMS["eps"],
                integration_policy=kjsj_policy,
                **p,
            )
            - KJSJ_unbarred_value(
                POINT_COORDS,
                1,
                0,
                2,
                singularity_policy=PARAMS["singularity_policy"],
                eps=PARAMS["eps"],
                integration_policy=kjsj_policy,
                **p,
            )
        ),
        "Kqbarq_xy": _abs(
            Kqbarq_unbarred_value(POINT_COORDS, 0, 1, 2, 3, **p)
            - Kqbarq_unbarred_value(POINT_COORDS, 1, 0, 2, 3, **p)
        ),
        "Kqbarq_zzp": _abs(
            Kqbarq_unbarred_value(POINT_COORDS, 0, 1, 2, 3, **p)
            - Kqbarq_unbarred_value(POINT_COORDS, 0, 1, 3, 2, **p)
        ),
        "KJSSJ_xy": _abs(
            KJSSJ_unbarred_value(POINT_COORDS, 0, 1, 2, 3, **p)
            - KJSSJ_unbarred_value(POINT_COORDS, 1, 0, 2, 3, **p)
        ),
        "KJSSJ_zzp": _abs(
            KJSSJ_unbarred_value(POINT_COORDS, 0, 1, 2, 3, **p)
            - KJSSJ_unbarred_value(POINT_COORDS, 0, 1, 3, 2, **p)
        ),
        "KJJSJ_xy_antisym": _abs(
            KJJSJ_unbarred_value(POINT_COORDS, 0, 1, 2, 3, **p)
            + KJJSJ_unbarred_value(POINT_COORDS, 0, 2, 1, 3, **p)
        ),
        "KJJSSJ_simultaneous_antisym": _abs(
            KJJSSJ_unbarred_value(POINT_COORDS, 0, 1, 2, 3, 4, **p)
            + KJJSSJ_unbarred_value(POINT_COORDS, 0, 2, 1, 4, 3, **p)
        ),
    }


def skeleton_metadata_check(kernels: dict) -> dict[str, object]:
    rng = np.random.default_rng(SEED)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(COORDS.shape[0])])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f, metadata_only=True)
    return {
        "metadata_only_pass": True,
        "sector_keys": sorted(terms.metadata["sectors"].keys()),
        "warnings": list(terms.metadata["warnings"]),
        "total_norms": terms.metadata["total_norms"],
    }


def skeleton_full_assembly_check() -> dict[str, object]:
    coords = COORDS[:3]
    rng = np.random.default_rng(SEED + 1)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(coords.shape[0])])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    kernels = physical_kernels_for_skeleton(coords, integration_policy=_policy(coords.shape[0]), **PARAMS)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f, metadata_only=False)
    return {
        "full_assembly_pass": True,
        "coordinate_count": int(coords.shape[0]),
        "term_shapes": {
            "K1": tuple(terms.K1.shape),
            "K2": tuple(terms.K2.shape),
            "K3": tuple(terms.K3.shape),
        },
        "norms": terms.norms(),
        "warnings": list(terms.metadata["warnings"]),
        "python_warnings": sorted(set(str(item.message) for item in caught)),
        "python_warning_count": len(caught),
    }


def write_report() -> dict[str, object]:
    policy = _policy(COORDS.shape[0])
    kernels = physical_kernels_for_skeleton(COORDS, integration_policy=policy, **PARAMS)
    stats = finite_kernel_stats(kernels)
    residuals = symmetry_residuals()
    skeleton = skeleton_metadata_check(kernels)
    full_assembly = skeleton_full_assembly_check()
    metadata = kernels["metadata"]
    cubic_convention = metadata.get("cubic_kernel_convention", {})

    lines = [
        "# Physical Kernel Integration Report",
        "",
        "Generated by `scripts/nlo_current/check_physical_kernel_integration.py`.",
        "",
        "## Scope",
        "",
        "This is a non-production dense coordinate diagnostic. It does not implement",
        "production evolution, score/Hessian-score model training, or physical",
        "positivity checks.",
        "",
        "## Source",
        "",
        "- formula notes: `docs/nlo_current/KLM_physical_kernel_formula_notes.md`",
        "- primary source: `references/WORKNLO.tex`",
        "",
        "## Configuration",
        "",
        f"- random_seed: `{SEED}`",
        f"- coordinates: `{COORDS.tolist()}`",
        f"- point_symmetry_coordinates: `{POINT_COORDS.tolist()}`",
        f"- Nc: `{PARAMS['Nc']}`",
        f"- nf: `{PARAMS['nf']}`",
        f"- alpha_s: `{PARAMS['alpha_s']}`",
        f"- singularity_policy: `{PARAMS['singularity_policy']}`",
        f"- eps: `{PARAMS['eps']}`",
        f"- KJSJ_policy: `{metadata['kjsj_integration_policy']}`",
        "",
        "## Implementation Status",
        "",
        f"- implemented_kernels: `{metadata['implemented_kernels']}`",
        f"- pending_kernels: `{metadata['pending_kernels']}`",
            "- KJSJ_status: implemented with an explicitly supplied diagnostic finite-sum policy.",
            "- KJSJ_policy_caveat: the finite sum is not a production regulator or",
            "  regulator-independent physical prescription.",
            f"- cubic_adapter_output: `{cubic_convention.get('adapter_output')}`",
            f"- cubic_normalization: `{cubic_convention.get('normalization')}`",
            "",
        "## Dense Kernel Statistics",
        "",
        "| kernel | shape | finite count | nonfinite count | finite norm |",
        "|---|---|---:|---:|---:|",
    ]
    for key in sorted(stats):
        item = stats[key]
        lines.append(
            f"| `{key}` | `{item['shape']}` | `{item['finite_count']}` | "
            f"`{item['nonfinite_count']}` | `{item['finite_norm']:.16e}` |"
        )
    lines.extend(
        [
            "",
            "## Symmetry Residuals",
            "",
            "| diagnostic | residual |",
            "|---|---:|",
        ]
    )
    for key, value in residuals.items():
        lines.append(f"| `{key}` | `{value:.16e}` |")
    lines.extend(
        [
            "",
            "## Skeleton Adapter Check",
            "",
            f"- metadata_only_pass: `{skeleton['metadata_only_pass']}`",
            f"- sector_keys: `{skeleton['sector_keys']}`",
            f"- warnings: `{skeleton['warnings']}`",
            f"- total_norms: `{skeleton['total_norms']}`",
            "",
            "The adapter output can be passed to `assemble_nlo_current_terms(...)` in",
            "metadata-only mode.",
            "",
            "## Full Dense Assembly Check",
            "",
            f"- full_assembly_pass: `{full_assembly['full_assembly_pass']}`",
            f"- coordinate_count: `{full_assembly['coordinate_count']}`",
            f"- term_shapes: `{full_assembly['term_shapes']}`",
            f"- norms: `{full_assembly['norms']}`",
            f"- warnings: `{full_assembly['warnings']}`",
            f"- python_warnings: `{full_assembly['python_warnings']}`",
            f"- python_warning_count: `{full_assembly['python_warning_count']}`",
            "",
            "The full dense assembly check uses the smallest valid non-degenerate",
            "coordinate subset with the same explicit non-production `eps` and",
            "`KJSJIntegrationPolicy` diagnostics.",
            "The physical adapter supplies KLM-normalized cubic coefficients, so raw",
            "physical cubic imaginary factors are not cast away during assembly.",
            "",
            "## Optional Dipole Recheck",
            "",
            "Run separately by `scripts/nlo_current/full_dipole_validation_physical_kernels.py`.",
            "",
            "## Warnings",
            "",
            "- `singularity_policy='eps'` is used as an explicit non-production",
            "  finite-value diagnostic regulator for dense arrays.",
            "- No barred/nonsinglet kernels are implemented in this diagnostic.",
            "- Physical-kernel positivity checks are future work.",
            "- The Pawula toy diagnostic does not prove physical NLO JIMWLK positivity",
            "  or non-positivity.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines))
    return {"stats": stats, "residuals": residuals, "skeleton": skeleton, "full_assembly": full_assembly}


def main() -> None:
    result = write_report()
    print(f"wrote {REPORT}")
    print("implemented kernels:", ", ".join(sorted(result["stats"].keys())))
    print("symmetry residuals:")
    for key, value in result["residuals"].items():
        print(f"  {key}: {value:.6e}")


if __name__ == "__main__":
    main()
