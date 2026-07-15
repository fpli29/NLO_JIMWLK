# Physical Kernel Workflow Start Status

## Execution Mode

- no-git mode: yes
- workflow scope: non-production physical KLM kernel diagnostics only
- production evolution code changes: not planned
- score/Hessian-score model training: not planned

## Source Availability

- `references/WORKNLO.tex`: available
- `references/1405.0418v2.pdf`: available for cross-check if needed

## Appendix A Baseline

`docs/nlo_current/dipole_validation_status.md` marks all five Appendix A dipole
targets as passed:

| sector | status |
|---|---|
| \(K_{JSJ}\) | passed |
| \(K_{JSSJ}\) | passed |
| \(K_{q\bar q}\) | passed |
| \(K_{JJSJ}\) | passed |
| \(K_{JJSSJ}\) | passed |

## Baseline Tests

Command:

```bash
python3 -m pytest tests/nlo_current -q
```

Result:

```text
89 passed in 10.19s
```

The count differs from the older workflow expectation of approximately
`84 passed` because the non-production Pawula/positivity diagnostic tests are
now included.

## Positivity Caveat

The existing Pawula/positivity diagnostic remains a caveat for this workflow.
Physical-kernel positivity checks are future work; the toy diagnostic does not
prove physical NLO JIMWLK positivity or non-positivity.
