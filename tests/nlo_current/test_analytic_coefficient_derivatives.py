from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.analytic_two_generator_derivatives import (  # noqa: E402
    analytic_dK2_KJSJ,
    analytic_dK2_KJSSJ,
    analytic_dK2_Kqbarq,
    analytic_dK2_Kqbarq_subtraction,
    analytic_dK2_Kqbarq_trace,
    kqbarq_subtraction_A_from_kernel,
    kqbarq_trace_A_from_kernel,
)
from nlo_current.coefficient_derivatives import fd_left_derivative_array  # noqa: E402
from nlo_current.nlo_current_skeleton import (  # noqa: E402
    assemble_kjsj_terms,
    assemble_kjssj_terms,
    assemble_kqbarq_terms,
    flatten_index,
    unflatten_index,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.two_generator_terms import kjssj_C_left_from_A  # noqa: E402


def _setup(seed=20260745, nsite=2):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U = np.stack([random_su3(rng) for _ in range(nsite)])
    S = np.stack([adjoint_from_fundamental(V, gens) for V in U])
    return rng, gens, f, U, S


def _S_builder(gens):
    return lambda U: np.stack([adjoint_from_fundamental(V, gens) for V in U])


def _flatten(C):
    nsite, n_color = C.shape[0], C.shape[1]
    out = np.zeros((nsite * n_color, nsite * n_color), dtype=np.result_type(C, complex))
    for x in range(nsite):
        for a in range(n_color):
            for y in range(nsite):
                for b in range(n_color):
                    out[flatten_index(x, a, n_color), flatten_index(y, b, n_color)] = C[x, a, y, b]
    return np.real_if_close(out)


def _fd_dK2_complex(K2_callback, U, gens, eps=1.0e-6):
    S_builder = _S_builder(gens)
    base = np.asarray(K2_callback(U, S_builder(U)))
    dim = base.shape[0]
    out = np.zeros(dim, dtype=np.result_type(base, complex))
    for b in range(dim):
        site, color = unflatten_index(b, gens.shape[0])
        derivative = fd_left_derivative_array(K2_callback, U, S_builder, site, color, gens, eps)
        out += derivative[:, b]
    return np.real_if_close(out)


def _random_kernel(shape, rng, scale=0.05):
    return scale * rng.normal(size=shape)


def test_analytic_dK2_KJSJ_matches_fd_oracle() -> None:
    rng, gens, f, U, S = _setup(seed=20260746)
    KJSJ = _random_kernel((2, 2, 2), rng)

    def callback(V, S_adj):
        return assemble_kjsj_terms(V, S_adj, KJSJ).K2

    analytic = analytic_dK2_KJSJ(S, KJSJ, f)
    fd = _fd_dK2_complex(callback, U, gens)
    np.testing.assert_allclose(analytic, fd, atol=2.0e-9, rtol=2.0e-7)


def test_analytic_dK2_KJSSJ_matches_fd_oracle() -> None:
    rng, gens, f, U, S = _setup(seed=20260747)
    KJSSJ = _random_kernel((2, 2, 2, 2), rng)

    def callback(V, S_adj):
        return assemble_kjssj_terms(V, S_adj, KJSSJ, f).K2

    analytic = analytic_dK2_KJSSJ(S, KJSSJ, f)
    fd = _fd_dK2_complex(callback, U, gens)
    np.testing.assert_allclose(analytic, fd, atol=2.0e-9, rtol=2.0e-7)


def test_analytic_dK2_Kqbarq_trace_subtraction_and_full_match_fd_oracle() -> None:
    rng, gens, f, U, S = _setup(seed=20260748)
    Kqbarq = _random_kernel((2, 2, 2, 2), rng)

    def trace_callback(V, S_adj):
        A = kqbarq_trace_A_from_kernel(V, Kqbarq, gens)
        return _flatten(kjssj_C_left_from_A(A, S_adj))

    def subtraction_callback(_V, S_adj):
        A = kqbarq_subtraction_A_from_kernel(S_adj, Kqbarq)
        return _flatten(kjssj_C_left_from_A(A, S_adj))

    def full_callback(V, S_adj):
        return assemble_kqbarq_terms(V, S_adj, Kqbarq, gens).K2

    trace_analytic = analytic_dK2_Kqbarq_trace(U, S, Kqbarq, gens, f)
    trace_fd = _fd_dK2_complex(trace_callback, U, gens)
    subtraction_analytic = analytic_dK2_Kqbarq_subtraction(U, S, Kqbarq, gens, f)
    subtraction_fd = _fd_dK2_complex(subtraction_callback, U, gens)
    full_analytic = analytic_dK2_Kqbarq(U, S, Kqbarq, gens, f)
    full_fd = _fd_dK2_complex(full_callback, U, gens)

    np.testing.assert_allclose(trace_analytic, trace_fd, atol=3.0e-9, rtol=3.0e-7)
    np.testing.assert_allclose(subtraction_analytic, subtraction_fd, atol=3.0e-9, rtol=3.0e-7)
    np.testing.assert_allclose(full_analytic, full_fd, atol=3.0e-9, rtol=3.0e-7)
    np.testing.assert_allclose(full_analytic, trace_analytic + subtraction_analytic, atol=1.0e-12)


def test_analytic_backend_supports_kjjsj_but_not_kjjssj() -> None:
    from nlo_current.physical_coefficient_derivatives import compute_physical_coefficient_derivatives
    from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig

    rng, gens, f, U, _ = _setup(seed=20260749)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]])
    from nlo_current.physical_kernels import KJSJIntegrationPolicy

    policy = KJSJIntegrationPolicy(
        quadrature_weights=np.ones(2) / 2,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
    )
    config = PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=2.0e-5,
        fd_eps_second=5.0e-4,
    )
    kjjsj = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="analytic",
        sector_filter=("KJJSJ",),
    )
    assert kjjsj.metadata["fallback_used"] is False
    assert kjjsj.metadata["implemented_analytic_sectors"] == ["KJJSJ"]
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


def test_structured_physical_analytic_backend_matches_fd_for_two_generator_sectors(monkeypatch) -> None:
    from nlo_current.physical_coefficient_derivatives import compute_physical_coefficient_derivatives
    from nlo_current.physical_kernels import KJSJIntegrationPolicy
    from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig

    _, gens, f, U, _ = _setup(seed=20260750)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]])
    policy = KJSJIntegrationPolicy(
        quadrature_weights=np.ones(2) / 2,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
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
        sector_filter=("KJSJ", "KJSSJ", "Kqbarq"),
    )
    fd = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="finite_difference",
        sector_filter=("KJSJ", "KJSSJ", "Kqbarq"),
    )
    np.testing.assert_allclose(analytic.dK2, fd.dK2, atol=2.0e-8, rtol=2.0e-6)
    np.testing.assert_allclose(analytic.LC_K3, 0.0, atol=1.0e-14)
    np.testing.assert_allclose(analytic.d2K3, 0.0, atol=1.0e-14)
    assert analytic.metadata["fallback_used"] is False

    import nlo_current.physical_coefficient_derivatives as pcd

    monkeypatch.setattr(
        pcd,
        "compute_dK2_fd",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("global FD was called")),
    )
    no_fd = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="analytic",
        sector_filter=("KJSSJ",),
    )
    assert np.linalg.norm(no_fd.dK2) > 0.0


def test_physical_velocity_and_closure_accept_analytic_two_generator_backend() -> None:
    from nlo_current.physical_density_closure import compare_physical_density_closure
    from nlo_current.physical_kernels import KJSJIntegrationPolicy
    from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig, evaluate_physical_nlo_velocity
    from nlo_current.test_densities import evaluate_test_density

    _, gens, f, U, _ = _setup(seed=20260751)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]])
    policy = KJSJIntegrationPolicy(
        quadrature_weights=np.ones(2) / 2,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
    )
    config = PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=2.0e-5,
        fd_eps_second=5.0e-4,
    )
    density = evaluate_test_density(U, gens, "dipole_trace")
    analytic_velocity = evaluate_physical_nlo_velocity(
        U,
        coords,
        gens,
        f,
        density.score,
        density.hessian_score,
        integration_policy=policy,
        config=config,
        sector_filter=("KJSSJ", "Kqbarq"),
        derivative_backend="analytic",
    )
    fd_velocity = evaluate_physical_nlo_velocity(
        U,
        coords,
        gens,
        f,
        density.score,
        density.hessian_score,
        integration_policy=policy,
        config=config,
        sector_filter=("KJSSJ", "Kqbarq"),
        derivative_backend="finite_difference",
    )
    np.testing.assert_allclose(
        analytic_velocity["velocity"],
        fd_velocity["velocity"],
        atol=2.0e-8,
        rtol=2.0e-6,
    )

    closure = compare_physical_density_closure(
        U,
        "dipole_trace",
        coords=coords,
        gens=gens,
        f=f,
        physical_policy=policy,
        config=config,
        derivative_backend="analytic",
        sector_mask=("KJSSJ", "Kqbarq"),
        fd_eps=1.0e-3,
        active_outer_indices=(0,),
    )
    assert closure.abs_residual < 5.0e-8
    assert closure.rel_residual < 5.0e-5


def test_hybrid_local_fd_is_explicitly_labeled_for_cubic_sectors() -> None:
    from nlo_current.physical_coefficient_derivatives import compute_physical_coefficient_derivatives
    from nlo_current.physical_kernels import KJSJIntegrationPolicy
    from nlo_current.physical_nlo_current import PhysicalNLOCurrentConfig

    _, gens, f, U, _ = _setup(seed=20260752)
    coords = np.array([[0.0, 0.0], [1.0, 0.2]])
    policy = KJSJIntegrationPolicy(
        quadrature_weights=np.ones(2) / 2,
        mu=1.3,
        exclude_coincident_labels=("x", "y", "z"),
    )
    config = PhysicalNLOCurrentConfig(
        nf=2,
        alpha_s=0.3,
        singularity_policy="eps",
        eps=1.0e-6,
        fd_eps_first=2.0e-5,
        fd_eps_second=5.0e-4,
    )
    hybrid = compute_physical_coefficient_derivatives(
        U,
        coords,
        gens,
        f,
        integration_policy=policy,
        config=config,
        backend="hybrid_local_fd",
        sector_filter=("KJJSJ",),
    )
    assert hybrid.backend == "hybrid_local_fd"
    assert hybrid.metadata["fallback_used"] is True
    assert hybrid.metadata["fd_components"] == ["KJJSJ"]
