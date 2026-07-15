from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.dipole_appendix_targets import (  # noqa: E402
    appendix_target_available,
    klm_normalized_cubic_direct_action,
    tilde_K_JJSSJ,
    target_KJJSJ_appendix,
    target_KJJSJ_appendix_real,
    target_KJJSJ_appendix_virtual,
    target_KJJSSJ_appendix,
    target_KJJSSJ_appendix_real,
    target_KJJSSJ_appendix_virtual,
    target_KJSSJ_appendix,
    target_KJSJ_appendix,
    target_Kqbarq_appendix,
    target_Kqbarq_trace_current_appendix,
)
from nlo_current.dipole_hamiltonian_action import (  # noqa: E402
    _action_LR_from_A,
    action_KJJSJ_direct,
    action_KJJSSJ_direct,
    action_KJSSJ_direct,
    action_KJSJ_direct,
    action_Kqbarq_direct,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.synthetic_kernels import (  # noqa: E402
    synthetic_kjssj_kernel,
    synthetic_kjsj_kernel,
    synthetic_kqbarq_kernel,
)
from nlo_current.three_generator_terms import (  # noqa: E402
    synthetic_kjjsj_kernel,
    synthetic_kjjssj_kernel,
)
from nlo_current.two_generator_terms import qbarq_trace_block  # noqa: E402


NC = 3


def _setup(seed: int = 62001, nsite: int = 3):
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    U_fund = np.stack([random_su3(rng) for _ in range(nsite)])
    S_adj = np.stack([adjoint_from_fundamental(U, gens) for U in U_fund])
    return rng, gens, f, U_fund, S_adj


def _appendix_compatible_kjsj(nsite, rng):
    kernel = synthetic_kjsj_kernel(nsite, rng)
    for x in range(nsite):
        kernel[x, x, :] = 0.0
    return 0.5 * (kernel + np.swapaxes(kernel, 0, 1))


def _appendix_compatible_kjssj(nsite, rng):
    kernel = synthetic_kjssj_kernel(nsite, rng)
    for x in range(nsite):
        kernel[x, x, :, :] = 0.0
    kernel = 0.5 * (kernel + np.swapaxes(kernel, 0, 1))
    return 0.5 * (kernel + np.swapaxes(kernel, 2, 3))


def _trace_word(U_fund, factors):
    product = np.eye(np.asarray(U_fund).shape[-1], dtype=complex)
    for dagger, site in factors:
        U_site = U_fund[site]
        product = product @ (U_site.conj().T if dagger else U_site)
    return np.trace(product)


def _kjssj_combined_bracket(U_fund, u, v, z, zp):
    s_uv = np.trace(U_fund[u].conj().T @ U_fund[v]) / NC
    s_u_zp = np.trace(U_fund[u].conj().T @ U_fund[zp]) / NC
    s_zp_z = np.trace(U_fund[zp].conj().T @ U_fund[z]) / NC
    s_z_v = np.trace(U_fund[z].conj().T @ U_fund[v]) / NC
    s_u_z = np.trace(U_fund[u].conj().T @ U_fund[z]) / NC
    return (
        (NC**3) * s_u_zp * s_zp_z * s_z_v
        - _trace_word(U_fund, ((False, v), (True, z), (False, zp), (True, u), (False, z), (True, zp)))
        - (NC**3) * s_u_z * s_z_v
        + NC * s_uv
    )


def _kjssj_subsection_partial_target(U_fund, KJSSJ, u, v):
    total = 0.0j
    for z in range(np.asarray(U_fund).shape[0]):
        for zp in range(np.asarray(U_fund).shape[0]):
            partial_bracket = (
                (NC**2)
                * (np.trace(U_fund[u].conj().T @ U_fund[zp]) / NC)
                * (np.trace(U_fund[zp].conj().T @ U_fund[z]) / NC)
                * (np.trace(U_fund[z].conj().T @ U_fund[v]) / NC)
                - _trace_word(U_fund, ((True, u), (False, z), (True, zp), (False, v), (True, z), (False, zp)))
                / NC
            )
            total += -KJSSJ[u, v, z, zp] * partial_bracket
    return total


def _kjjssj_pure_eight_real_target(U_fund, KJJSSJ, u, v):
    total = 0.0j
    for z in range(np.asarray(U_fund).shape[0]):
        for zp in range(np.asarray(U_fund).shape[0]):
            combo = (
                KJJSSJ[u, u, u, z, zp]
                - KJJSSJ[u, v, u, z, zp]
                + KJJSSJ[u, v, v, z, zp]
                - KJJSSJ[u, u, v, z, zp]
                + KJJSSJ[v, u, v, z, zp]
                - KJJSSJ[v, u, u, z, zp]
                + KJJSSJ[v, v, u, z, zp]
                - KJJSSJ[v, v, v, z, zp]
            )
            triple = (
                (NC**3)
                * (np.trace(U_fund[z].conj().T @ U_fund[v]) / NC)
                * (np.trace(U_fund[zp].conj().T @ U_fund[z]) / NC)
                * (np.trace(U_fund[u].conj().T @ U_fund[zp]) / NC)
            )
            total += 0.5j * combo * triple / NC
    return total


def _kjjssj_tilde_real_target(U_fund, KJJSSJ, u, v):
    total = 0.0j
    for z in range(np.asarray(U_fund).shape[0]):
        for zp in range(np.asarray(U_fund).shape[0]):
            triple = (
                (NC**3)
                * (np.trace(U_fund[z].conj().T @ U_fund[v]) / NC)
                * (np.trace(U_fund[zp].conj().T @ U_fund[z]) / NC)
                * (np.trace(U_fund[u].conj().T @ U_fund[zp]) / NC)
            )
            trace = _trace_word(U_fund, ((False, v), (True, z), (False, zp), (True, u), (False, z), (True, zp)))
            total += tilde_K_JJSSJ(KJJSSJ, u, v, z, zp) * (triple - trace) / NC
    return total


def _qbarq_trace_current_A(U_fund, Kqbarq, gens):
    nsite = np.asarray(U_fund).shape[0]
    trace_blocks = np.empty((nsite, nsite, gens.shape[0], gens.shape[0]), dtype=complex)
    for z in range(nsite):
        for zp in range(nsite):
            trace_blocks[z, zp] = qbarq_trace_block(U_fund[z], U_fund[zp], gens)
    return np.einsum("xyuv,uvab->xyab", Kqbarq, trace_blocks, optimize=True)


def test_zero_kernels_give_zero_direct_action_for_all_sectors() -> None:
    _, gens, f, U_fund, S_adj = _setup(seed=62002, nsite=2)
    u, v = 0, 1

    np.testing.assert_allclose(action_KJSJ_direct(U_fund, S_adj, np.zeros((2, 2, 2)), u, v, gens), 0.0)
    np.testing.assert_allclose(
        action_KJSSJ_direct(U_fund, S_adj, np.zeros((2, 2, 2, 2)), u, v, f, gens),
        0.0,
    )
    np.testing.assert_allclose(
        action_Kqbarq_direct(U_fund, S_adj, np.zeros((2, 2, 2, 2)), u, v, gens),
        0.0,
    )
    np.testing.assert_allclose(
        action_KJJSJ_direct(U_fund, S_adj, np.zeros((2, 2, 2, 2)), u, v, f, gens),
        0.0,
    )
    np.testing.assert_allclose(
        action_KJJSSJ_direct(U_fund, S_adj, np.zeros((2, 2, 2, 2, 2)), u, v, f, gens),
        0.0,
    )


def test_direct_actions_are_linear_in_each_sector_kernel() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62003, nsite=2)
    u, v = 0, 1
    sector_cases = [
        (
            action_KJSJ_direct,
            synthetic_kjsj_kernel(2, rng),
            synthetic_kjsj_kernel(2, rng),
            (U_fund, S_adj, u, v, gens),
        ),
        (
            action_KJSSJ_direct,
            synthetic_kjssj_kernel(2, rng),
            synthetic_kjssj_kernel(2, rng),
            (U_fund, S_adj, u, v, f, gens),
        ),
        (
            action_Kqbarq_direct,
            synthetic_kqbarq_kernel(2, rng),
            synthetic_kqbarq_kernel(2, rng),
            (U_fund, S_adj, u, v, gens),
        ),
        (
            action_KJJSJ_direct,
            synthetic_kjjsj_kernel(2, rng, xy_symmetry="unconstrained"),
            synthetic_kjjsj_kernel(2, rng, xy_symmetry="unconstrained"),
            (U_fund, S_adj, u, v, f, gens),
        ),
        (
            action_KJJSSJ_direct,
            synthetic_kjjssj_kernel(2, rng, klm_antisym=False),
            synthetic_kjjssj_kernel(2, rng, klm_antisym=False),
            (U_fund, S_adj, u, v, f, gens),
        ),
    ]

    for action, kernel_a, kernel_b, trailing_args in sector_cases:
        lhs = action(trailing_args[0], trailing_args[1], kernel_a + kernel_b, *trailing_args[2:])
        rhs = action(trailing_args[0], trailing_args[1], kernel_a, *trailing_args[2:]) + action(
            trailing_args[0],
            trailing_args[1],
            kernel_b,
            *trailing_args[2:],
        )
        np.testing.assert_allclose(lhs, rhs, atol=1e-12, rtol=1e-12)


def test_KJSJ_direct_action_matches_appendix_A_target() -> None:
    rng, gens, _, U_fund, S_adj = _setup(seed=62004, nsite=3)
    KJSJ = _appendix_compatible_kjsj(3, rng)

    assert appendix_target_available("KJSJ") is True
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJSJ_direct(U_fund, S_adj, KJSJ, u, v, gens)
        target = target_KJSJ_appendix(U_fund, KJSJ, u, v)
        np.testing.assert_allclose(direct, target, atol=2e-13, rtol=2e-12)


def test_unavailable_appendix_targets_are_marked_pending() -> None:
    _, _, _, _, _ = _setup(seed=62005, nsite=2)

    assert appendix_target_available("KJSSJ") is True
    assert appendix_target_available("Kqbarq") is True
    assert appendix_target_available("KJJSJ") is True
    assert appendix_target_available("KJJSSJ") is True


def test_KJSSJ_appendix_target_zero_kernel_and_linearity() -> None:
    rng, _, _, U_fund, _ = _setup(seed=62016, nsite=3)
    zero = np.zeros((3, 3, 3, 3))
    kernel_a = _appendix_compatible_kjssj(3, rng)
    kernel_b = _appendix_compatible_kjssj(3, rng)

    np.testing.assert_allclose(target_KJSSJ_appendix(U_fund, zero, 0, 1), 0.0, atol=1e-13, rtol=1e-13)
    for u, v in ((0, 1), (1, 2), (0, 2)):
        lhs = target_KJSSJ_appendix(U_fund, kernel_a + kernel_b, u, v)
        rhs = target_KJSSJ_appendix(U_fund, kernel_a, u, v) + target_KJSSJ_appendix(U_fund, kernel_b, u, v)
        np.testing.assert_allclose(lhs, rhs, atol=2e-13, rtol=2e-12)


def test_KJSSJ_full_target_matches_appendix_A_direct_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62017, nsite=3)
    KJSSJ = _appendix_compatible_kjssj(3, rng)

    assert appendix_target_available("KJSSJ") is True
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJSSJ_direct(U_fund, S_adj, KJSSJ, u, v, f, gens)
        target = target_KJSSJ_appendix(U_fund, KJSSJ, u, v)
        np.testing.assert_allclose(direct, target, atol=4e-13, rtol=4e-12)


def test_KJSSJ_subsection_partial_target_fails_as_full_sector_target() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62018, nsite=3)
    KJSSJ = _appendix_compatible_kjssj(3, rng)

    failures = []
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJSSJ_direct(U_fund, S_adj, KJSSJ, u, v, f, gens)
        partial = _kjssj_subsection_partial_target(U_fund, KJSSJ, u, v)
        failures.append(abs(direct - partial))

    assert max(failures) > 1e-8


def test_KJSSJ_isolated_target_excludes_tildeK_contamination() -> None:
    rng, _, _, U_fund, _ = _setup(seed=62019, nsite=3)
    KJSSJ = _appendix_compatible_kjssj(3, rng)
    fake_tilde = rng.normal(size=(3, 3))

    contaminated_deltas = []
    for u, v in ((0, 1), (1, 2), (0, 2)):
        isolated = target_KJSSJ_appendix(U_fund, KJSSJ, u, v)
        contaminated = isolated
        for z in range(3):
            for zp in range(3):
                contaminated += fake_tilde[z, zp] * _kjssj_combined_bracket(U_fund, u, v, z, zp) / NC
        contaminated_deltas.append(abs(contaminated - isolated))

    assert max(contaminated_deltas) > 1e-8


def test_Kqbarq_subsection_target_matches_trace_current_component_only() -> None:
    rng, gens, _, U_fund, _ = _setup(seed=62008, nsite=3)
    Kqbarq = synthetic_kqbarq_kernel(3, rng)
    A_trace = _qbarq_trace_current_A(U_fund, Kqbarq, gens)

    mismatch_to_full_seen = False
    for u, v in ((0, 1), (1, 2), (0, 2)):
        trace_current = _action_LR_from_A(U_fund, A_trace, u, v, gens)
        subsection_target = target_Kqbarq_trace_current_appendix(U_fund, Kqbarq, u, v)
        full_target = target_Kqbarq_appendix(U_fund, Kqbarq, u, v)

        np.testing.assert_allclose(trace_current, subsection_target, atol=3e-13, rtol=3e-12)
        mismatch_to_full_seen |= abs(subsection_target - full_target) > 1e-10

    assert mismatch_to_full_seen


def test_Kqbarq_full_target_matches_appendix_A_direct_action() -> None:
    rng, gens, _, U_fund, S_adj = _setup(seed=62009, nsite=3)
    Kqbarq = synthetic_kqbarq_kernel(3, rng)

    assert appendix_target_available("Kqbarq") is True
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_Kqbarq_direct(U_fund, S_adj, Kqbarq, u, v, gens)
        target = target_Kqbarq_appendix(U_fund, Kqbarq, u, v)
        np.testing.assert_allclose(direct, target, atol=4e-13, rtol=4e-12)


def test_Kqbarq_trace_product_target_matches_four_site_direct_action() -> None:
    rng, gens, _, U_fund, S_adj = _setup(seed=62028, nsite=4)
    raw = rng.normal(size=(4, 4, 4, 4))
    Kqbarq = 0.5 * (raw + np.swapaxes(raw, 0, 1))
    Kqbarq = 0.5 * (Kqbarq + np.swapaxes(Kqbarq, 2, 3))

    for u, v in ((0, 1), (0, 2), (1, 3)):
        A_trace = _qbarq_trace_current_A(U_fund, Kqbarq, gens)
        trace_direct = _action_LR_from_A(U_fund, A_trace, u, v, gens)
        trace_target = target_Kqbarq_trace_current_appendix(U_fund, Kqbarq, u, v, gens=gens)
        full_direct = action_Kqbarq_direct(U_fund, S_adj, Kqbarq, u, v, gens)
        full_target = target_Kqbarq_appendix(U_fund, Kqbarq, u, v, gens=gens)

        np.testing.assert_allclose(trace_direct, trace_target, atol=2e-13, rtol=2e-12)
        np.testing.assert_allclose(full_direct, full_target, atol=2e-13, rtol=2e-12)


def test_Kqbarq_subsection_target_fails_as_full_sector_target() -> None:
    rng, gens, _, U_fund, S_adj = _setup(seed=62010, nsite=3)
    Kqbarq = np.zeros((3, 3, 3, 3))
    diagonal_weights = rng.normal(size=(3, 3, 3))
    for x in range(3):
        for y in range(3):
            for z in range(3):
                Kqbarq[x, y, z, z] = diagonal_weights[x, y, z]

    full_direct = action_Kqbarq_direct(U_fund, S_adj, Kqbarq, 0, 1, gens)
    full_target = target_Kqbarq_appendix(U_fund, Kqbarq, 0, 1)
    subsection_target = target_Kqbarq_trace_current_appendix(U_fund, Kqbarq, 0, 1)

    np.testing.assert_allclose(full_direct, full_target, atol=3e-13, rtol=3e-13)
    np.testing.assert_allclose(full_direct, 0.0, atol=3e-13, rtol=3e-13)
    assert abs(subsection_target - full_direct) > 1e-10


def test_KJJSJ_appendix_target_zero_kernel_and_linearity() -> None:
    rng, _, _, U_fund, _ = _setup(seed=62011, nsite=3)
    zero = np.zeros((3, 3, 3, 3))
    kernel_a = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")
    kernel_b = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")

    np.testing.assert_allclose(target_KJJSJ_appendix(U_fund, zero, 0, 1), 0.0, atol=1e-13, rtol=1e-13)
    for u, v in ((0, 1), (1, 2), (0, 2)):
        lhs = target_KJJSJ_appendix(U_fund, kernel_a + kernel_b, u, v)
        rhs = target_KJJSJ_appendix(U_fund, kernel_a, u, v) + target_KJJSJ_appendix(U_fund, kernel_b, u, v)
        np.testing.assert_allclose(lhs, rhs, atol=2e-13, rtol=2e-12)


def test_KJJSJ_real_target_matches_calibrated_direct_real_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62012, nsite=3)
    KJJSJ = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")

    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct_real = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens, virtual_scale=0.0)
        target_real = target_KJJSJ_appendix_real(U_fund, KJJSJ, u, v)
        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct_real),
            target_real,
            atol=4e-13,
            rtol=4e-12,
        )


def test_KJJSJ_virtual_target_matches_calibrated_direct_virtual_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62013, nsite=3)
    KJJSJ = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")

    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct_real = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens, virtual_scale=0.0)
        direct_full = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens)
        target_virtual = target_KJJSJ_appendix_virtual(U_fund, KJJSJ, u, v)
        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct_full - direct_real),
            target_virtual,
            atol=4e-13,
            rtol=4e-12,
        )


def test_KJJSJ_full_target_matches_calibrated_direct_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62014, nsite=3)
    KJJSJ = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")

    assert appendix_target_available("KJJSJ") is True
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens)
        target = target_KJJSJ_appendix(U_fund, KJJSJ, u, v)
        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct),
            target,
            atol=4e-13,
            rtol=4e-12,
        )


def test_KJJSJ_raw_uncalibrated_direct_action_fails_generic_appendix_target() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62015, nsite=3)
    KJJSJ = synthetic_kjjsj_kernel(3, rng, xy_symmetry="antisymmetric")

    failures = []
    calibrated_residuals = []
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, u, v, f, gens)
        target = target_KJJSJ_appendix(U_fund, KJJSJ, u, v)
        failures.append(abs(direct - target))
        calibrated_residuals.append(abs(klm_normalized_cubic_direct_action(direct) - target))

    assert max(failures) > 1e-8
    assert max(calibrated_residuals) < 4e-13


def test_KJJSSJ_appendix_target_zero_kernel_and_linearity() -> None:
    rng, _, _, U_fund, _ = _setup(seed=62020, nsite=3)
    zero = np.zeros((3, 3, 3, 3, 3))
    kernel_a = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)
    kernel_b = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    np.testing.assert_allclose(target_KJJSSJ_appendix(U_fund, zero, 0, 1), 0.0, atol=1e-13, rtol=1e-13)
    for u, v in ((0, 1), (1, 2), (0, 2)):
        lhs = target_KJJSSJ_appendix(U_fund, kernel_a + kernel_b, u, v)
        rhs = target_KJJSSJ_appendix(U_fund, kernel_a, u, v) + target_KJJSSJ_appendix(U_fund, kernel_b, u, v)
        np.testing.assert_allclose(lhs, rhs, atol=3e-13, rtol=3e-12)


def test_KJJSSJ_real_target_matches_calibrated_direct_real_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62021, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct_real = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens, virtual_scale=0.0)
        target_real = target_KJJSSJ_appendix_real(U_fund, KJJSSJ, u, v)
        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct_real),
            target_real,
            atol=5e-13,
            rtol=5e-12,
        )


def test_KJJSSJ_virtual_target_matches_calibrated_direct_virtual_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62022, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct_real = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens, virtual_scale=0.0)
        direct_full = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens)
        target_virtual = target_KJJSSJ_appendix_virtual(U_fund, KJJSSJ, u, v)
        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct_full - direct_real),
            target_virtual,
            atol=5e-13,
            rtol=5e-12,
        )


def test_KJJSSJ_full_target_matches_calibrated_direct_action() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62023, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    assert appendix_target_available("KJJSSJ") is True
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens)
        target = target_KJJSSJ_appendix(U_fund, KJJSSJ, u, v)
        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct),
            target,
            atol=5e-13,
            rtol=5e-12,
        )


def test_KJJSSJ_raw_uncalibrated_direct_action_fails_generic_appendix_target() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62024, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    failures = []
    calibrated_residuals = []
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens)
        target = target_KJJSSJ_appendix(U_fund, KJJSSJ, u, v)
        failures.append(abs(direct - target))
        calibrated_residuals.append(abs(klm_normalized_cubic_direct_action(direct) - target))

    assert max(failures) > 1e-8
    assert max(calibrated_residuals) < 5e-13


def test_KJJSSJ_virtual_one_third_factor_is_required_for_appendix_target() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62025, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    changes = []
    for u, v in ((0, 1), (1, 2), (0, 2)):
        direct_real = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens, virtual_scale=0.0)
        direct_full = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens)
        direct_no_factor = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, u, v, f, gens, virtual_scale=3.0)
        target_virtual = target_KJJSSJ_appendix_virtual(U_fund, KJJSSJ, u, v)
        calibrated_virtual_no_factor = klm_normalized_cubic_direct_action(direct_no_factor - direct_real)

        np.testing.assert_allclose(
            klm_normalized_cubic_direct_action(direct_full - direct_real),
            target_virtual,
            atol=5e-13,
            rtol=5e-12,
        )
        changes.append(abs(calibrated_virtual_no_factor - target_virtual))

    assert max(changes) > 1e-8


def test_KJJSSJ_tildeK_is_nonzero_for_generic_synthetic_kernel() -> None:
    rng, _, _, _, _ = _setup(seed=62026, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    values = [
        tilde_K_JJSSJ(KJJSSJ, x, y, z, zp)
        for x in range(3)
        for y in range(3)
        for z in range(3)
        for zp in range(3)
    ]
    assert max(abs(value) for value in values) > 1e-10


def test_KJJSSJ_real_target_exercises_pure_eight_and_tilde_terms() -> None:
    rng, _, _, U_fund, _ = _setup(seed=62027, nsite=3)
    KJJSSJ = synthetic_kjjssj_kernel(3, rng, klm_antisym=True)

    eight_values = []
    tilde_values = []
    for u, v in ((0, 1), (1, 2), (0, 2)):
        eight = _kjjssj_pure_eight_real_target(U_fund, KJJSSJ, u, v)
        tilde = _kjjssj_tilde_real_target(U_fund, KJJSSJ, u, v)
        target = target_KJJSSJ_appendix_real(U_fund, KJJSSJ, u, v)
        np.testing.assert_allclose(eight + tilde, target, atol=3e-13, rtol=3e-12)
        eight_values.append(abs(eight))
        tilde_values.append(abs(tilde))

    assert max(eight_values) > 1e-10
    assert max(tilde_values) > 1e-10


def test_qbarq_zprime_equals_z_subtraction_gives_zero_action() -> None:
    rng, gens, _, U_fund, S_adj = _setup(seed=62006, nsite=3)
    Kqbarq = np.zeros((3, 3, 3, 3))
    diagonal_weights = rng.normal(size=(3, 3, 3))
    for x in range(3):
        for y in range(3):
            for z in range(3):
                Kqbarq[x, y, z, z] = diagonal_weights[x, y, z]

    action = action_Kqbarq_direct(U_fund, S_adj, Kqbarq, 0, 1, gens)

    np.testing.assert_allclose(action, 0.0, atol=2e-13, rtol=2e-13)


def test_cubic_virtual_one_third_factor_sensitivity() -> None:
    rng, gens, f, U_fund, S_adj = _setup(seed=62007, nsite=2)
    KJJSJ = synthetic_kjjsj_kernel(2, rng, xy_symmetry="unconstrained")
    KJJSSJ = synthetic_kjjssj_kernel(2, rng, klm_antisym=False)

    kjjsj_default = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, 0, 1, f, gens)
    kjjsj_no_factor = action_KJJSJ_direct(U_fund, S_adj, KJJSJ, 0, 1, f, gens, virtual_scale=3.0)
    kjjssj_default = action_KJJSSJ_direct(U_fund, S_adj, KJJSSJ, 0, 1, f, gens)
    kjjssj_no_factor = action_KJJSSJ_direct(
        U_fund,
        S_adj,
        KJJSSJ,
        0,
        1,
        f,
        gens,
        virtual_scale=3.0,
    )

    assert abs(kjjsj_default - kjjsj_no_factor) > 1e-10
    assert abs(kjjssj_default - kjjssj_no_factor) > 1e-10
