from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
REPORT_PATH = ROOT / "reports" / "nlo_current" / "pawula_positivity_diagnostic_report.md"
GRID_SIZE = 96
LENGTH = 2.0 * np.pi
DY_VALUES = (2.5e-4, 5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3)


@dataclass(frozen=True)
class ToyCase:
    name: str
    description: str
    density_name: str
    dy: float
    k1_description: str
    k2_description: str
    k3_description: str

    def density(self, theta: np.ndarray) -> np.ndarray:
        if self.density_name == "smooth_positive":
            return np.exp(1.4 * np.cos(theta - 0.3)) + 0.08
        if self.density_name == "near_zero_notch":
            return (1.0 + np.cos(theta - 0.2)) ** 2 + 0.05 * (
                1.0 + np.sin(3.0 * theta)
            ) ** 2 + 1.0e-8
        if self.density_name == "variable_positive":
            return (
                np.exp(1.1 * np.cos(theta + 0.2))
                * (1.0 + 0.2 * np.sin(2.0 * theta - 0.1))
                + 0.05
            )
        raise ValueError(f"unknown density_name: {self.density_name}")

    def coefficients(self, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        zeros = np.zeros_like(theta)
        if self.name == "lo_like_diffusion":
            return zeros, 0.08 * np.ones_like(theta), zeros
        if self.name == "pure_third_order":
            return zeros, zeros, 10.0 * np.ones_like(theta)
        if self.name == "mixed_nlo_like":
            return zeros, 0.08 * np.ones_like(theta), 0.15 * np.ones_like(theta)
        if self.name == "variable_coefficients":
            k1 = 0.01 * np.sin(theta)
            k2 = 0.06 * (1.0 + 0.35 * np.cos(theta - 0.1))
            k3 = 0.25 * (1.0 + 0.25 * np.sin(2.0 * theta))
            return k1, k2, k3
        raise ValueError(f"unknown case name: {self.name}")


def periodic_derivative_matrices(
    n_grid: int = GRID_SIZE,
    length: float = LENGTH,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray, np.ndarray]:
    """Return periodic grid, spacing, and central derivative matrices."""

    h = length / n_grid
    theta = h * np.arange(n_grid)
    d1 = np.zeros((n_grid, n_grid), dtype=float)
    d2 = np.zeros((n_grid, n_grid), dtype=float)
    for i in range(n_grid):
        d1[i, (i + 1) % n_grid] = 1.0 / (2.0 * h)
        d1[i, (i - 1) % n_grid] = -1.0 / (2.0 * h)
        d2[i, (i + 1) % n_grid] = 1.0 / h**2
        d2[i, i] = -2.0 / h**2
        d2[i, (i - 1) % n_grid] = 1.0 / h**2
    d3 = d1 @ d2
    return theta, h, d1, d2, d3


def normalize_density(w: np.ndarray, h: float) -> np.ndarray:
    clipped = np.maximum(np.asarray(w, dtype=float), 0.0)
    mass = h * float(np.sum(clipped))
    if mass <= 0.0:
        raise ValueError("density has zero mass after clipping")
    return clipped / mass


def build_generator_matrix(
    k1: np.ndarray,
    k2: np.ndarray,
    k3: np.ndarray,
    d1: np.ndarray,
    d2: np.ndarray,
    d3: np.ndarray,
) -> np.ndarray:
    return (
        -d1 @ np.diag(k1)
        + 0.5 * d2 @ np.diag(k2)
        - (1.0 / 6.0) * d3 @ np.diag(k3)
    )


def negative_mass(w: np.ndarray, h: float) -> float:
    return h * float(np.sum(np.maximum(-w, 0.0)))


def near_zero_mask(w: np.ndarray, fraction_of_max: float = 0.05) -> np.ndarray:
    mask = w <= fraction_of_max * float(np.max(w))
    if not np.any(mask):
        mask[int(np.argmin(w))] = True
    return mask


def positive_maximum_principle_diagnostic(
    generator: np.ndarray,
    w: np.ndarray,
) -> dict[str, float | int | bool]:
    rhs = generator @ w
    mask = near_zero_mask(w)
    rhs_near = rhs[mask]
    return {
        "near_zero_count": int(np.sum(mask)),
        "min_w_near_zero": float(np.min(w[mask])),
        "min_rhs_near_zero": float(np.min(rhs_near)),
        "negative_rhs_near_zero_count": int(np.sum(rhs_near < -1.0e-12)),
        "warning": bool(np.any(rhs_near < -1.0e-12)),
    }


def offdiagonal_sign_diagnostic(generator: np.ndarray) -> dict[str, float | int | bool]:
    n = generator.shape[0]
    offdiag = generator[~np.eye(n, dtype=bool)]
    negative = offdiag < -1.0e-14
    positive = offdiag > 1.0e-14
    return {
        "offdiagonal_count": int(offdiag.size),
        "negative_offdiagonal_count": int(np.sum(negative)),
        "positive_offdiagonal_count": int(np.sum(positive)),
        "negative_offdiagonal_fraction": float(np.sum(negative) / offdiag.size),
        "min_offdiagonal": float(np.min(offdiag)),
        "max_offdiagonal": float(np.max(offdiag)),
        "warning": bool(np.any(negative)),
    }


def default_case_definitions() -> list[ToyCase]:
    return [
        ToyCase(
            name="lo_like_diffusion",
            description="LO-like second-order diffusion with K3=0 and K2>0.",
            density_name="smooth_positive",
            dy=1.0e-3,
            k1_description="0",
            k2_description="0.08",
            k3_description="0",
        ),
        ToyCase(
            name="pure_third_order",
            description="Pure third-order generalized FP term with K1=K2=0.",
            density_name="near_zero_notch",
            dy=1.0e-3,
            k1_description="0",
            k2_description="0",
            k3_description="10",
        ),
        ToyCase(
            name="mixed_nlo_like",
            description="Positive diffusion plus a small third-order term.",
            density_name="near_zero_notch",
            dy=1.0e-3,
            k1_description="0",
            k2_description="0.08",
            k3_description="0.15",
        ),
        ToyCase(
            name="variable_coefficients",
            description="Synthetic variable K1, positive variable K2, and nonzero variable K3.",
            density_name="variable_positive",
            dy=1.0e-3,
            k1_description="0.01 sin(theta)",
            k2_description="0.06 [1 + 0.35 cos(theta - 0.1)]",
            k3_description="0.25 [1 + 0.25 sin(2 theta)]",
        ),
    ]


def run_case(
    case: ToyCase,
    n_grid: int = GRID_SIZE,
    dy_values: tuple[float, ...] = DY_VALUES,
) -> dict[str, object]:
    theta, h, d1, d2, d3 = periodic_derivative_matrices(n_grid)
    w0 = normalize_density(case.density(theta), h)
    k1, k2, k3 = case.coefficients(theta)
    generator = build_generator_matrix(k1, k2, k3, d1, d2, d3)
    rhs = generator @ w0
    w_after = w0 + case.dy * rhs
    mass_before = h * float(np.sum(w0))
    mass_after = h * float(np.sum(w_after))
    scaling = []
    for dy in dy_values:
        stepped = w0 + dy * rhs
        scaling.append(
            {
                "dY": float(dy),
                "min_w_after": float(np.min(stepped)),
                "negative_mass": negative_mass(stepped, h),
            }
        )

    max_principle = positive_maximum_principle_diagnostic(generator, w0)
    offdiag = offdiagonal_sign_diagnostic(generator)
    warning = (
        negative_mass(w_after, h) > 1.0e-15
        or bool(max_principle["warning"])
        or bool(offdiag["warning"])
    )
    return {
        "name": case.name,
        "description": case.description,
        "density_name": case.density_name,
        "n_grid": int(n_grid),
        "spacing": float(h),
        "dY": float(case.dy),
        "k1_description": case.k1_description,
        "k2_description": case.k2_description,
        "k3_description": case.k3_description,
        "mass_before": mass_before,
        "mass_after": mass_after,
        "normalization_drift": mass_after - mass_before,
        "min_w_initial": float(np.min(w0)),
        "min_w_after": float(np.min(w_after)),
        "negative_mass_after": negative_mass(w_after, h),
        "negative_mass_scaling": scaling,
        "positive_maximum_principle": max_principle,
        "offdiagonal_signs": offdiag,
        "positivity_warning": bool(warning),
    }


def run_all_cases(n_grid: int = GRID_SIZE) -> list[dict[str, object]]:
    return [run_case(case, n_grid=n_grid) for case in default_case_definitions()]


def _format_scaling(entries: list[dict[str, float]]) -> str:
    rows = ["| dY | min(W_after) | negative mass |", "|---:|---:|---:|"]
    for entry in entries:
        rows.append(
            f"| `{entry['dY']:.4e}` | `{entry['min_w_after']:.16e}` | "
            f"`{entry['negative_mass']:.16e}` |"
        )
    return "\n".join(rows)


def write_report(results: list[dict[str, object]], report_path: Path = REPORT_PATH) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    n_grid = int(results[0]["n_grid"])
    spacing = float(results[0]["spacing"])
    lines = [
        "# Pawula Positivity Toy Diagnostic Report",
        "",
        "This report was generated by `scripts/nlo_current/check_pawula_positivity_toy.py`.",
        "",
        "## Scope",
        "",
        "This is a non-production one-dimensional periodic-grid diagnostic. It does not",
        "implement physical kernels, production evolution, or score/Hessian-score model",
        "training.",
        "",
        "## Grid And Derivative Scheme",
        "",
        f"- grid_size: `{n_grid}`",
        "- domain: `[0, 2 pi)`",
        f"- spacing: `{spacing:.16e}`",
        "- derivative_scheme: second-order central periodic finite-difference matrices",
        "- third_derivative_matrix: `D3 = D1 @ D2`",
        "- time_stepper: one explicit Euler diagnostic step",
        "",
        "## Case Summary",
        "",
        "| case | density | K1 | K2 | K3 | dY | min(W_after) | negative mass | norm drift | warning |",
        "|---|---|---|---|---|---:|---:|---:|---:|---|",
    ]
    for result in results:
        lines.append(
            f"| `{result['name']}` | `{result['density_name']}` | "
            f"`{result['k1_description']}` | `{result['k2_description']}` | "
            f"`{result['k3_description']}` | `{result['dY']:.4e}` | "
            f"`{result['min_w_after']:.16e}` | "
            f"`{result['negative_mass_after']:.16e}` | "
            f"`{result['normalization_drift']:.16e}` | "
            f"`{result['positivity_warning']}` |"
        )
    lines.extend(["", "## Negative Mass Scaling", ""])
    for result in results:
        lines.extend(
            [
                f"### `{result['name']}`",
                "",
                _format_scaling(result["negative_mass_scaling"]),  # type: ignore[arg-type]
                "",
            ]
        )
    lines.extend(["## Positive Maximum-Principle Diagnostic", ""])
    lines.extend(
        [
            "| case | near-zero count | min W near zero | min RHS near zero | negative RHS count | warning |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for result in results:
        diag = result["positive_maximum_principle"]
        assert isinstance(diag, dict)
        lines.append(
            f"| `{result['name']}` | `{diag['near_zero_count']}` | "
            f"`{diag['min_w_near_zero']:.16e}` | "
            f"`{diag['min_rhs_near_zero']:.16e}` | "
            f"`{diag['negative_rhs_near_zero_count']}` | `{diag['warning']}` |"
        )
    lines.extend(["", "## Off-Diagonal Sign Diagnostic", ""])
    lines.extend(
        [
            "| case | negative off-diagonal | positive off-diagonal | negative fraction | min offdiag | max offdiag | warning |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for result in results:
        diag = result["offdiagonal_signs"]
        assert isinstance(diag, dict)
        lines.append(
            f"| `{result['name']}` | `{diag['negative_offdiagonal_count']}` | "
            f"`{diag['positive_offdiagonal_count']}` | "
            f"`{diag['negative_offdiagonal_fraction']:.16e}` | "
            f"`{diag['min_offdiagonal']:.16e}` | "
            f"`{diag['max_offdiagonal']:.16e}` | `{diag['warning']}` |"
        )
    lines.extend(
        [
            "",
            "## Conclusion",
            "",
            "The LO-like diffusion case has no negative off-diagonal generator entries and",
            "shows no short-step negative mass in this setup. The pure third-order case",
            "shows negative off-diagonal entries, a positive maximum-principle warning,",
            "and negative mass for a short Euler step on the constructed near-zero",
            "density. Mixed and variable-coefficient cases retain warning signs from",
            "the third-order term even when the selected short step stays nonnegative.",
            "",
            "This demonstrates positivity risk of finite third-order generalized FP operators. "
            "It does not prove physical NLO JIMWLK positivity or non-positivity.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines))


def main() -> None:
    results = run_all_cases()
    write_report(results)
    for result in results:
        print(
            f"{result['name']}: min_after={result['min_w_after']:.6e}, "
            f"negative_mass={result['negative_mass_after']:.6e}, "
            f"warning={result['positivity_warning']}"
        )
    print(f"wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
