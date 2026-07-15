# Derivation Summary Workflow Start Status

Workflow: `codex_nlo_current_derivation_summary_workflow.md`

Mode: no-git mode. `git rev-parse --is-inside-work-tree` reported that this
workspace is not a Git repository.

Scope: documentation only. This workflow does not implement physical kernels,
does not modify production evolution code, and does not train score or
Hessian-score models.

## Required Artifacts

All Phase 0 documentation and report artifacts listed in the workflow were
present at start:

- derivation docs under `docs/nlo_current/`: present
- validation and diagnostic reports under `reports/nlo_current/`: present
- `reports/nlo_current/file_manifest.md`: present

## Test Status

Command run before summarizing:

```bash
python3 -m pytest tests/nlo_current -q
```

Result:

```text
84 passed in 9.42s
```

## Appendix A Dipole Target Status

The latest `reports/nlo_current/full_dipole_validation_report.md` marks all
five Appendix A dipole sectors as passed:

| sector | status | max residual |
|---|---|---:|
| `KJSJ` | passed | `1.6172982426630268e-16` |
| `KJSSJ` | passed | `1.3877787807814457e-16` |
| `Kqbarq` | passed | `5.6325111571547335e-17` |
| `KJJSJ` | passed | `1.1801832636420706e-15` |
| `KJJSSJ` | passed | `1.9087382418229414e-15` |

## Notes

No source-code edits are required by this workflow. The planned edits are
limited to Markdown summaries and the no-git file manifest.

