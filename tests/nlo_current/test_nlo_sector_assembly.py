from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.nlo_current_skeleton import assemble_nlo_current_terms, combined_dim  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.synthetic_kernels import synthetic_kernels_all  # noqa: E402


def _assembly_inputs(seed: int = 42001, nsite: int = 2):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(nsite)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    kernels = synthetic_kernels_all(nsite, rng)
    return rng, gens, f, U_fund, S_adj, kernels


def test_assemble_all_sectors_metadata_records_sector_map_and_commutators() -> None:
    _, gens, f, U_fund, S_adj, kernels = _assembly_inputs()

    terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f, include_commutators=True)

    sectors = terms.metadata["sectors"]
    for key in ("KJSJ", "KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ"):
        assert key in sectors
        assert sectors[key]["norms"]["K1"] >= 0.0
        assert sectors[key]["norms"]["K2"] >= 0.0
        assert sectors[key]["norms"]["K3"] >= 0.0
    assert "commutators" in terms.metadata
    assert terms.metadata["commutators"]["included"] is True
    assert sectors["KJJSJ"]["include_commutators"] is True
    assert sectors["KJJSSJ"]["include_commutators"] is True


def test_assembled_output_shapes_are_dense_combined_site_color_shapes() -> None:
    _, gens, f, U_fund, S_adj, kernels = _assembly_inputs(seed=42002)
    terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f)
    dim = combined_dim(nsite=2)

    assert terms.K1.shape == (dim,)
    assert terms.K2.shape == (dim, dim)
    assert terms.K3.shape == (dim, dim, dim)
    terms.validate_shapes()


def test_commutator_toggle_changes_lower_order_corrections_but_not_cubic_part() -> None:
    _, gens, f, U_fund, S_adj, kernels = _assembly_inputs(seed=42003)

    with_comm = assemble_nlo_current_terms(
        U_fund, S_adj, kernels, gens, f, include_commutators=True
    )
    without_comm = assemble_nlo_current_terms(
        U_fund, S_adj, kernels, gens, f, include_commutators=False
    )

    assert with_comm.metadata["commutators"]["included"] is True
    assert without_comm.metadata["commutators"]["included"] is False
    assert with_comm.metadata["sectors"]["KJJSJ"]["include_commutators"] is True
    assert without_comm.metadata["sectors"]["KJJSJ"]["include_commutators"] is False
    np.testing.assert_allclose(with_comm.K3, without_comm.K3, atol=1e-12, rtol=1e-12)
    lower_order_delta = np.linalg.norm(with_comm.K1 - without_comm.K1) + np.linalg.norm(
        with_comm.K2 - without_comm.K2
    )
    assert lower_order_delta > 1e-10


def test_missing_kernel_handling_records_zero_missing_sector_terms() -> None:
    _, gens, f, U_fund, S_adj, kernels = _assembly_inputs(seed=42004)
    subset = {"KJSJ": kernels["KJSJ"]}

    terms = assemble_nlo_current_terms(U_fund, S_adj, subset, gens, f)

    for key in ("KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ"):
        assert terms.metadata["sectors"][key]["missing"] is True
        assert terms.metadata["sectors"][key]["norms"] == {"K1": 0.0, "K2": 0.0, "K3": 0.0}
    assert all(
        f"missing kernel: {key}" in terms.metadata["warnings"]
        for key in ("KJSSJ", "Kqbarq", "KJJSJ", "KJJSSJ")
    )
    assert np.linalg.norm(terms.K2) > 0.0
    assert np.linalg.norm(terms.K1) == 0.0
    assert np.linalg.norm(terms.K3) == 0.0


def test_generic_ordered_K2_is_not_forcibly_symmetrized() -> None:
    _, gens, f, U_fund, S_adj, kernels = _assembly_inputs(seed=42005)

    terms = assemble_nlo_current_terms(U_fund, S_adj, kernels, gens, f)
    asymmetry = np.linalg.norm(terms.K2 - terms.K2.T)

    assert asymmetry > 1e-10
