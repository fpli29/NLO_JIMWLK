#!/usr/bin/env python3
"""Coefficient-level diagnostics for cubic commutator corrections."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nlo_current.cubic_commutator_terms import (  # noqa: E402
    canonicalize_cubic_block_terms,
    cubic_block_terms_from_LLR,
    cubic_block_terms_from_LRR,
    cubic_block_terms_from_virtual_LLL,
    cubic_block_terms_from_virtual_RRR,
    llr_right_to_left,
    lrr_right_to_left,
    virtual_rrr_right_to_left,
)
from nlo_current.su3_adjoint import (  # noqa: E402
    adjoint_from_fundamental,
    random_su3,
    structure_constants,
    su3_generators_fundamental,
)
from nlo_current.three_generator_terms import (  # noqa: E402
    kjjsj_A_LLR_from_kernel,
    kjjsj_B_LRR_from_kernel,
    kjjsj_V_virtual_from_kernel,
    kjjssj_A_LLR_from_kernel,
    kjjssj_B_LRR_from_kernel,
    kjjssj_V_virtual_from_kernel,
    synthetic_kjjsj_kernel,
    synthetic_kjjssj_kernel,
)


def _combine(*term_dicts):
    out = {}
    for terms in term_dicts:
        for word, coeff in terms.items():
            out[word] = out.get(word, 0.0) + float(coeff)
            if abs(out[word]) < 1e-14:
                del out[word]
    return out


def _norm(terms):
    return float(np.sqrt(sum(float(coeff) ** 2 for coeff in terms.values())))


def _build_jjssj_terms(S_adj, kernel, f):
    A = kjjsj_A_LLR_from_kernel(S_adj, kernel, f)
    B = kjjsj_B_LRR_from_kernel(S_adj, kernel, f)
    V = kjjsj_V_virtual_from_kernel(kernel, f)
    C_llr = llr_right_to_left(A, S_adj)
    C_lrr = lrr_right_to_left(B, S_adj)
    C_rrr = virtual_rrr_right_to_left(V, S_adj, sign=-1.0)
    return _combine(
        cubic_block_terms_from_LLR(C_llr),
        cubic_block_terms_from_LRR(C_lrr),
        cubic_block_terms_from_virtual_LLL(V),
        cubic_block_terms_from_virtual_RRR(C_rrr),
    )


def _build_jjssj_double_s_terms(S_adj, kernel, f):
    A = kjjssj_A_LLR_from_kernel(S_adj, kernel, f)
    B = kjjssj_B_LRR_from_kernel(S_adj, kernel, f)
    V = kjjssj_V_virtual_from_kernel(kernel, f)
    C_llr = llr_right_to_left(A, S_adj)
    C_lrr = lrr_right_to_left(B, S_adj)
    C_rrr = virtual_rrr_right_to_left(V, S_adj, sign=-1.0)
    return _combine(
        cubic_block_terms_from_LLR(C_llr),
        cubic_block_terms_from_LRR(C_lrr),
        cubic_block_terms_from_virtual_LLL(V),
        cubic_block_terms_from_virtual_RRR(C_rrr),
    )


def _diagnose(name, raw_terms, f):
    cubic, quadratic, linear = canonicalize_cubic_block_terms(raw_terms, f)
    return {
        "name": name,
        "raw_count": len(raw_terms),
        "cubic_count": len(cubic),
        "quadratic_count": len(quadratic),
        "linear_count": len(linear),
        "raw_norm": _norm(raw_terms),
        "cubic_norm": _norm(cubic),
        "quadratic_norm": _norm(quadratic),
        "linear_norm": _norm(linear),
        "comm_nonzero": _norm(quadratic) > 1e-12 or _norm(linear) > 1e-12,
    }


def write_report(results, seed, nsite):
    report = ROOT / "reports" / "nlo_current" / "cubic_commutator_corrections_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Cubic Commutator Corrections Report",
        "",
        f"random_seed: {seed}",
        f"nsite: {nsite}",
        "",
        "Raw left-basis cubic words were canonicalized using",
        "`L_x^a L_x^b = L_x^b L_x^a + f^{abc} L_x^c`.",
        "",
    ]
    for result in results:
        lines.extend(
            [
                f"## {result['name']}",
                "",
                f"raw_cubic_terms: {result['raw_count']}",
                f"canonical_cubic_terms: {result['cubic_count']}",
                f"quadratic_commutator_terms: {result['quadratic_count']}",
                f"linear_commutator_terms: {result['linear_count']}",
                "",
                "| class | coefficient norm |",
                "|---|---:|",
                f"| raw cubic | {result['raw_norm']:.16e} |",
                f"| canonical cubic | {result['cubic_norm']:.16e} |",
                f"| quadratic commutator | {result['quadratic_norm']:.16e} |",
                f"| linear commutator | {result['linear_norm']:.16e} |",
                "",
                f"commutator_corrections_nonzero: {result['comm_nonzero']}",
                "",
            ]
        )
    any_nonzero = any(result["comm_nonzero"] for result in results)
    lines.extend(
        [
            "## Conclusion",
            "",
            f"any_commutator_corrections_nonzero: {any_nonzero}",
            "",
        ]
    )
    if any_nonzero:
        lines.append(
            "Commutator-induced lower-order terms are present in these coincident-sector diagnostics."
        )
    else:
        lines.append(
            "The synthetic diagnostic produced numerically zero lower-order terms; this is not a proof."
        )
    lines.append("")
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main():
    seed = 20260811
    nsite = 2
    rng = np.random.default_rng(seed)
    gens = su3_generators_fundamental()
    f = structure_constants(gens)
    S_adj = np.stack([adjoint_from_fundamental(random_su3(rng), gens) for _ in range(nsite)])

    kjjsj_kernel = synthetic_kjjsj_kernel(nsite, rng, xy_symmetry="antisymmetric")
    kjjssj_kernel = synthetic_kjjssj_kernel(nsite, rng, klm_antisym=True)

    results = [
        _diagnose("K_JJSJ", _build_jjssj_terms(S_adj, kjjsj_kernel, f), f),
        _diagnose("K_JJSSJ", _build_jjssj_double_s_terms(S_adj, kjjssj_kernel, f), f),
    ]
    report = write_report(results, seed, nsite)
    for result in results:
        print(
            f"{result['name']}: raw={result['raw_count']} cubic={result['cubic_count']} "
            f"quad_norm={result['quadratic_norm']:.16e} "
            f"lin_norm={result['linear_norm']:.16e}"
        )
    print(f"wrote {report}")


if __name__ == "__main__":
    main()

