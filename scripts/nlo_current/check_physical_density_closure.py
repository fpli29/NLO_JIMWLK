from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.physical_density_closure import compare_physical_density_closure  # noqa: E402
from nlo_current.physical_kernels import KJSJIntegrationPolicy  # noqa: E402
from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig  # noqa: E402
from nlo_current.su3_adjoint import random_su3, structure_constants, su3_generators_fundamental  # noqa: E402


def _policy(nsite: int) -> KJSJIntegrationPolicy:
    return KJSJIntegrationPolicy(
        quadrature_weights=np.ones(nsite) / nsite,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
        description="density closure diagnostic finite z-prime policy",
    )


def _config(fd_eps_first: float, fd_eps_second: float) -> PhysicalNLOCurrentConfig:
    return PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=fd_eps_first,
        fd_eps_second=fd_eps_second,
    )


def _serial(value):
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def run_density_closure_scan() -> dict:
    rng = np.random.default_rng(20260724)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]], dtype=float)
    U = np.stack([random_su3(rng) for _ in range(2)])
    steps = [2.0e-3, 1.0e-3, 5.0e-4, 2.5e-4]
    densities = ["single_link_trace", "dipole_trace"]
    rows = []
    for density in densities:
        for step in steps:
            start = time.perf_counter()
            result = compare_physical_density_closure(
                U,
                density,
                coords=coords,
                gens=gens,
                f=f,
                physical_policy=_policy(2),
                config=_config(fd_eps_first=max(2.0e-5, step / 20.0), fd_eps_second=max(5.0e-4, step)),
                fd_eps=step,
                sector_mask=("KJSSJ", "Kqbarq"),
                active_outer_indices=(0,),
                include_by_sector=(density == "dipole_trace" and step == steps[1]),
            )
            rows.append(
                {
                    "density": density,
                    "fd_eps": step,
                    "direct": result.direct_value,
                    "current": result.current_value,
                    "abs_residual": result.abs_residual,
                    "rel_residual": result.rel_residual,
                    "real_residual": float((result.direct_value - result.current_value).real),
                    "imag_residual": float((result.direct_value - result.current_value).imag),
                    "sector_residuals": result.sector_residuals,
                    "wall_time_s": time.perf_counter() - start,
                    "active_outer_indices": [0],
                }
            )

    toggle_base = compare_physical_density_closure(
        U,
        "dipole_trace",
        coords=coords,
        gens=gens,
        f=f,
        physical_policy=_policy(2),
        config=_config(2.0e-5, 5.0e-4),
        fd_eps=1.0e-3,
        sector_mask=("KJSSJ", "Kqbarq"),
        active_outer_indices=(0,),
    )
    toggles = {
        "omit_coefficient_derivatives": compare_physical_density_closure(
            U,
            "dipole_trace",
            coords=coords,
            gens=gens,
            f=f,
            physical_policy=_policy(2),
            config=_config(2.0e-5, 5.0e-4),
            fd_eps=1.0e-3,
            sector_mask=("KJSSJ", "Kqbarq"),
            active_outer_indices=(0,),
            omit_coefficient_derivatives=True,
        ).abs_residual,
        "omit_hessian_score": compare_physical_density_closure(
            U,
            "dipole_trace",
            coords=coords,
            gens=gens,
            f=f,
            physical_policy=_policy(2),
            config=_config(2.0e-5, 5.0e-4),
            fd_eps=1.0e-3,
            sector_mask=("KJSSJ", "Kqbarq"),
            active_outer_indices=(0,),
            omit_hessian_score=True,
        ).abs_residual,
        "omit_commutators": compare_physical_density_closure(
            U,
            "dipole_trace",
            coords=coords,
            gens=gens,
            f=f,
            physical_policy=_policy(2),
            config=_config(2.0e-5, 5.0e-4),
            fd_eps=1.0e-3,
            sector_mask=("KJSSJ", "Kqbarq"),
            active_outer_indices=(0,),
            omit_commutators=True,
        ).abs_residual,
        "remove_cubic_normalization": compare_physical_density_closure(
            U,
            "dipole_trace",
            coords=coords,
            gens=gens,
            f=f,
            physical_policy=_policy(2),
            config=_config(2.0e-5, 5.0e-4),
            fd_eps=1.0e-3,
            sector_mask=("KJSSJ", "Kqbarq"),
            active_outer_indices=(0,),
            remove_cubic_normalization=True,
        ).abs_residual,
        "baseline_abs_residual": toggle_base.abs_residual,
    }
    return {
        "coordinate_set": coords.tolist(),
        "wilson_seed": 20260724,
        "densities": densities,
        "steps": steps,
        "rows": rows,
        "toggles": toggles,
        "caveat": (
            "The scan uses active_outer_indices=(0,) to keep the dense finite-"
            "difference diagnostic small. The default physical scan is restricted "
            "to KJSSJ and Kqbarq; cubic toggle coverage is in the unit tests and "
            "report. This is not production evolution."
        ),
    }


def main() -> None:
    data = run_density_closure_scan()
    print(json.dumps(data, indent=2, default=_serial, sort_keys=True))


if __name__ == "__main__":
    main()
