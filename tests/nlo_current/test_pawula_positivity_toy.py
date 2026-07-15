from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = ROOT / "scripts" / "nlo_current" / "check_pawula_positivity_toy.py"
SPEC = importlib.util.spec_from_file_location("pawula_toy", SCRIPT_PATH)
assert SPEC is not None
pawula_toy = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = pawula_toy
SPEC.loader.exec_module(pawula_toy)


def _case(name: str):
    for case in pawula_toy.default_case_definitions():
        if case.name == name:
            return case
    raise AssertionError(f"missing case {name}")


def _finite_values(obj) -> list[float]:
    values = []
    if isinstance(obj, dict):
        for value in obj.values():
            values.extend(_finite_values(value))
    elif isinstance(obj, list):
        for value in obj:
            values.extend(_finite_values(value))
    elif isinstance(obj, bool):
        pass
    elif isinstance(obj, (int, float, np.floating)):
        values.append(float(obj))
    return values


def test_periodic_derivative_matrices_conserve_constants() -> None:
    _, _, d1, d2, d3 = pawula_toy.periodic_derivative_matrices(64)
    ones = np.ones(64)

    np.testing.assert_allclose(d1 @ ones, 0.0, atol=1e-12)
    np.testing.assert_allclose(d2 @ ones, 0.0, atol=1e-12)
    np.testing.assert_allclose(d3 @ ones, 0.0, atol=1e-10)
    np.testing.assert_allclose(ones @ d1, 0.0, atol=1e-12)
    np.testing.assert_allclose(ones @ d2, 0.0, atol=1e-12)
    np.testing.assert_allclose(ones @ d3, 0.0, atol=1e-10)


def test_normalization_drift_is_near_zero_for_periodic_cases() -> None:
    for case in pawula_toy.default_case_definitions():
        result = pawula_toy.run_case(case)
        assert abs(result["normalization_drift"]) < 1e-12


def test_diffusion_case_has_no_short_step_negativity() -> None:
    result = pawula_toy.run_case(_case("lo_like_diffusion"))

    assert result["min_w_after"] >= -1e-14
    assert result["negative_mass_after"] <= 1e-14
    assert result["offdiagonal_signs"]["negative_offdiagonal_count"] == 0


def test_pure_k3_case_shows_positivity_warning() -> None:
    result = pawula_toy.run_case(_case("pure_third_order"))

    assert result["positivity_warning"] is True
    assert result["negative_mass_after"] > 0.0
    assert result["positive_maximum_principle"]["negative_rhs_near_zero_count"] > 0
    assert result["offdiagonal_signs"]["negative_offdiagonal_count"] > 0


def test_diagnostics_report_finite_values() -> None:
    for result in pawula_toy.run_all_cases():
        values = _finite_values(result)
        assert values
        assert np.all(np.isfinite(values))
