from __future__ import annotations

import sys
import warnings
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.analytic_cubic_derivatives import (  # noqa: E402
    KJJSJ_BLOCKS,
    KJJSJSectorData,
    analytic_LB_K3_KJJSJ,
    analytic_LC_K3_KJJSJ,
    analytic_d2K3_KJJSJ,
    analytic_dK1_comm_KJJSJ,
    analytic_dK2_comm_KJJSJ,
    analytic_first_derivatives_KJJSJ,
    kjjsj_terms_from_blocks,
)
from nlo_current.coefficient_derivatives import compute_all_coefficient_derivatives_fd  # noqa: E402
from nlo_current.nlo_current_skeleton import NLOCurrentTerms  # noqa: E402
from nlo_current.nlo_velocity_evaluator import evaluate_velocity_from_terms  # noqa: E402
from nlo_current.physical_current_divergence import evaluate_current_divergence  # noqa: E402
from nlo_current.physical_density_closure import raw_cubic_terms_transform  # noqa: E402
from nlo_current.physical_density_operator import evaluate_direct_density_operator  # noqa: E402
from nlo_current.physical_kernels import KJSJIntegrationPolicy  # noqa: E402
from nlo_current.physical_coefficient_derivatives import compute_physical_coefficient_derivatives  # noqa: E402
from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig  # noqa: E402
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    left_perturb,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.test_densities import evaluate_test_density  # noqa: E402


def _S_builder(gens):
    return lambda U: np.stack([adjoint_from_fundamental(V, gens) for V in U])


def _terms_from_data(data: KJJSJSectorData, *, include_commutators=True, blocks=None) -> NLOCurrentTerms:
    K1, K2, K3 = kjjsj_terms_from_blocks(data, blocks=blocks)
    if not include_commutators:
        K1 = np.zeros_like(K1)
        K2 = np.zeros_like(K2)
    return NLOCurrentTerms(K1=K1, K2=K2, K3=K3, metadata={"sectors": {"KJJSJ": {}}})


def _perturb(U, flat_index, gens, eps):
    n_color = gens.shape[0]
    site, color = divmod(int(flat_index), n_color)
    out = np.array(U, copy=True)
    out[site] = left_perturb(out[site], color, eps, gens)
    return out


def _data_from_U(U, K, gens, f):
    return KJJSJSectorData(S_adj=_S_builder(gens)(U), KJJSJ=K, f=f)


@lru_cache(maxsize=1)
def _oracle():
    rng = np.random.default_rng(20260756)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U = np.stack([random_su3(rng) for _ in range(2)])
    K = 0.03 * rng.normal(size=(2, 2, 2, 2))
    data = _data_from_U(U, K, gens, f)
    first = analytic_first_derivatives_KJJSJ(sector_data=data, return_diagnostics=True)
    d2 = analytic_d2K3_KJJSJ(sector_data=data, return_diagnostics=True)
    fd = compute_all_coefficient_derivatives_fd(
        lambda V, S: _terms_from_data(KJJSJSectorData(S, K, f)).K2,
        lambda V, S: _terms_from_data(KJJSJSectorData(S, K, f)).K3,
        U,
        _S_builder(gens),
        gens,
        eps_first=2.0e-5,
        eps_second=5.0e-4,
    )
    return {"gens": gens, "f": f, "U": U, "K": K, "data": data, "first": first, "d2": d2, "fd": fd}


def _assert_close(name, analytic, fd, atol, rtol):
    diff = np.asarray(analytic) - np.asarray(fd)
    assert np.max(np.abs(diff)) < atol, name
    assert np.linalg.norm(diff) / (np.linalg.norm(fd) + 1.0e-30) < rtol, name


def test_kjjsj_full_LC_LB_d2_and_comm_derivatives_match_fd_oracle() -> None:
    o = _oracle()
    _assert_close("dK2_comm", o["first"]["dK2_comm"], o["fd"]["dK2"], 2.0e-10, 2.0e-8)
    _assert_close("LC_K3", o["first"]["LC_K3"], o["fd"]["dK3_first"]["LC_K3_ABC"], 2.0e-10, 2.0e-8)
    _assert_close("LB_K3", o["first"]["LB_K3"], o["fd"]["dK3_first"]["LB_K3_ABC"], 2.0e-10, 2.0e-8)
    _assert_close("d2K3", o["d2"]["total"], o["fd"]["d2K3"], 2.0e-7, 2.0e-5)


def _fd_first_entry(o, block, deriv_index, entry, eps=2.0e-5):
    gens, f, U, K = o["gens"], o["f"], o["U"], o["K"]

    def val(V):
        data = _data_from_U(V, K, gens, f)
        return _terms_from_data(data, blocks=(block,)).K3[entry]

    return (val(_perturb(U, deriv_index, gens, eps)) - val(_perturb(U, deriv_index, gens, -eps))) / (
        2.0 * eps
    )


def _fd_second_entry(o, block, first_index, second_index, entry, eps=5.0e-4):
    gens, U = o["gens"], o["U"]

    def inner(V):
        plus = _perturb(V, second_index, gens, eps)
        minus = _perturb(V, second_index, gens, -eps)
        return (
            _terms_from_data(_data_from_U(plus, o["K"], gens, o["f"]), blocks=(block,)).K3[entry]
            - _terms_from_data(_data_from_U(minus, o["K"], gens, o["f"]), blocks=(block,)).K3[entry]
        ) / (2.0 * eps)

    return (inner(_perturb(U, first_index, gens, eps)) - inner(_perturb(U, first_index, gens, -eps))) / (
        2.0 * eps
    )


def _max_entry(arr):
    idx = np.unravel_index(np.argmax(np.abs(arr)), arr.shape)
    return idx, arr[idx]


def test_kjjsj_llr_lrr_virtual_first_derivative_blocks_have_fd_spot_checks() -> None:
    o = _oracle()
    data = o["data"]
    for block in KJJSJ_BLOCKS:
        if block in {"virtual_LLL", "virtual_RRR"}:
            np.testing.assert_allclose(o["first"]["by_block"]["LC_K3"][block], 0.0, atol=1.0e-14)
            np.testing.assert_allclose(o["first"]["by_block"]["LB_K3"][block], 0.0, atol=1.0e-14)
            continue
        found = False
        for deriv_index in range(data.S_adj.shape[0] * data.S_adj.shape[1]):
            _, _, K3_d = kjjsj_terms_from_blocks(data, blocks=(block,), derivative="first", first_index=deriv_index)
            entry, value = _max_entry(K3_d)
            if abs(value) > 1.0e-8:
                np.testing.assert_allclose(
                    value,
                    _fd_first_entry(o, block, deriv_index, entry),
                    atol=2.0e-9,
                    rtol=2.0e-7,
                )
                found = True
                break
        assert found, block


def test_kjjsj_llr_lrr_virtual_ordered_d2_blocks_have_fd_spot_checks() -> None:
    o = _oracle()
    data = o["data"]
    dim = data.S_adj.shape[0] * data.S_adj.shape[1]
    for block in KJJSJ_BLOCKS:
        if block in {"virtual_LLL", "virtual_RRR"}:
            np.testing.assert_allclose(o["d2"]["by_block"][block], 0.0, atol=1.0e-14)
            continue
        found = False
        for first_index in range(dim):
            for second_index in range(dim):
                _, _, K3_dd = kjjsj_terms_from_blocks(
                    data,
                    blocks=(block,),
                    derivative="second",
                    first_index=first_index,
                    second_index=second_index,
                )
                entry, value = _max_entry(K3_dd)
                if abs(value) > 1.0e-8:
                    np.testing.assert_allclose(
                        value,
                        _fd_second_entry(o, block, first_index, second_index, entry),
                        atol=2.0e-7,
                        rtol=2.0e-5,
                    )
                    found = True
                    break
            if found:
                break
        assert found, block


def test_same_site_order_is_not_reversed_and_distinct_sites_commute_where_expected() -> None:
    o = _oracle()
    data = o["data"]
    _, _, K3_01 = kjjsj_terms_from_blocks(data, derivative="second", first_index=0, second_index=1)
    _, _, K3_10 = kjjsj_terms_from_blocks(data, derivative="second", first_index=1, second_index=0)
    assert np.linalg.norm(K3_01 - K3_10) > 1.0e-8

    _, _, K3_distinct_0_8 = kjjsj_terms_from_blocks(data, derivative="second", first_index=0, second_index=8)
    _, _, K3_distinct_8_0 = kjjsj_terms_from_blocks(data, derivative="second", first_index=8, second_index=0)
    np.testing.assert_allclose(K3_distinct_0_8, K3_distinct_8_0, atol=1.0e-12, rtol=1.0e-12)


def test_kjjsj_commutator_status_and_cubic_normalization_character() -> None:
    o = _oracle()
    data = o["data"]
    dK2 = analytic_dK2_comm_KJJSJ(sector_data=data, return_diagnostics=True)
    k1_status = analytic_dK1_comm_KJJSJ(sector_data=data)
    assert np.linalg.norm(dK2["total"]) > 0.0
    assert k1_status["status"] in {"structurally_zero", "nonzero"}
    assert np.max(np.abs(np.imag(o["first"]["LC_K3"]))) < 1.0e-14


def test_backend_analytic_kjjsj_no_fd_fallback_and_kjjssj_still_pending(monkeypatch) -> None:
    import nlo_current.physical_coefficient_derivatives as pcd

    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    rng = np.random.default_rng(20260757)
    U = np.stack([random_su3(rng) for _ in range(2)])
    coords = np.array([[0.0, 0.0], [1.0, 0.2]])
    policy = KJSJIntegrationPolicy(
        quadrature_weights=np.ones(2) / 2,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
    )
    config = PhysicalNLOCurrentConfig(nf=2, alpha_s=0.3, singularity_policy="eps", eps=1.0e-6)
    monkeypatch.setattr(
        pcd,
        "compute_all_coefficient_derivatives_fd",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("global FD was called")),
    )
    result = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="analytic",
        sector_filter=("KJJSJ",),
    )
    assert result.metadata["fallback_used"] is False
    assert result.metadata["implemented_analytic_sectors"] == ["KJJSJ"]
    with np.testing.assert_raises(NotImplementedError):
        compute_physical_coefficient_derivatives(
            U,
            coords,
            gens,
            f,
            integration_policy=policy,
            config=config,
            backend="analytic",
            sector_filter=("KJJSSJ",),
        )


def test_kjjsj_velocity_and_sparse_density_closure_with_analytic_derivatives() -> None:
    o = _oracle()
    data = o["data"]
    terms = _terms_from_data(data)
    derivatives = {
        "dK2": o["first"]["dK2_comm"],
        "dK3_first": {"LC_K3_ABC": o["first"]["LC_K3"], "LB_K3_ABC": o["first"]["LB_K3"]},
        "d2K3": o["d2"]["total"],
    }
    fd_derivatives = {
        "dK2": o["fd"]["dK2"],
        "dK3_first": o["fd"]["dK3_first"],
        "d2K3": o["fd"]["d2K3"],
    }
    density = evaluate_test_density(o["U"], o["gens"], "dipole_trace")
    analytic_velocity = evaluate_velocity_from_terms(
        terms, density.score, density.hessian_score, **derivatives
    )["velocity"]
    fd_velocity = evaluate_velocity_from_terms(
        terms, density.score, density.hessian_score, **fd_derivatives
    )["velocity"]
    np.testing.assert_allclose(analytic_velocity, fd_velocity, atol=2.0e-8, rtol=2.0e-6)

    density_builder = lambda V: evaluate_test_density(V, o["gens"], "dipole_trace")
    terms_builder = lambda V: _terms_from_data(_data_from_U(V, o["K"], o["gens"], o["f"]))
    derivatives_builder = lambda V: {
        "dK2": analytic_first_derivatives_KJJSJ(sector_data=_data_from_U(V, o["K"], o["gens"], o["f"]))[
            "dK2_comm"
        ],
        "dK3_first": {
            "LC_K3_ABC": analytic_first_derivatives_KJJSJ(
                sector_data=_data_from_U(V, o["K"], o["gens"], o["f"])
            )["LC_K3"],
            "LB_K3_ABC": analytic_first_derivatives_KJJSJ(
                sector_data=_data_from_U(V, o["K"], o["gens"], o["f"])
            )["LB_K3"],
        },
        "d2K3": analytic_d2K3_KJJSJ(sector_data=_data_from_U(V, o["K"], o["gens"], o["f"])),
    }
    direct = evaluate_direct_density_operator(
        o["U"],
        density_builder,
        terms,
        gens=o["gens"],
        fd_eps=1.0e-3,
        terms_builder=terms_builder,
        active_outer_indices=(0,),
    )
    current = evaluate_current_divergence(
        o["U"],
        density_builder,
        density.score,
        density.hessian_score,
        terms,
        derivatives,
        gens=o["gens"],
        fd_eps=1.0e-3,
        terms_builder=terms_builder,
        derivatives_builder=derivatives_builder,
        density_builder=density_builder,
        active_outer_indices=(0,),
    )
    assert abs(direct.value - current.value) < 5.0e-7

    no_hessian = evaluate_current_divergence(
        o["U"],
        density_builder,
        density.score,
        density.hessian_score,
        terms,
        derivatives,
        gens=o["gens"],
        fd_eps=1.0e-3,
        terms_builder=terms_builder,
        derivatives_builder=derivatives_builder,
        density_builder=density_builder,
        active_outer_indices=(0,),
        omit_hessian_score=True,
    )
    assert abs(direct.value - no_hessian.value) > abs(direct.value - current.value) + 1.0e-7


def test_kjjsj_no_complex_warning_and_raw_cubic_toggle_remains_meaningful() -> None:
    o = _oracle()
    terms = _terms_from_data(o["data"])
    density_builder = lambda V: evaluate_test_density(V, o["gens"], "dipole_trace")
    density = density_builder(o["U"])
    derivatives = {
        "dK2": o["first"]["dK2_comm"],
        "dK3_first": {"LC_K3_ABC": o["first"]["LC_K3"], "LB_K3_ABC": o["first"]["LB_K3"]},
        "d2K3": o["d2"]["total"],
    }
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        raw = evaluate_current_divergence(
            o["U"],
            density_builder,
            density.score,
            density.hessian_score,
            terms,
            derivatives,
            gens=o["gens"],
            fd_eps=1.0e-3,
            terms_builder=lambda V: _terms_from_data(_data_from_U(V, o["K"], o["gens"], o["f"])),
            derivatives_builder=lambda _V: derivatives,
            density_builder=density_builder,
            active_outer_indices=(0,),
            terms_transform=raw_cubic_terms_transform,
        )
    assert not any("ComplexWarning" in item.category.__name__ for item in caught)
    assert abs(raw.value.imag) > 1.0e-10
