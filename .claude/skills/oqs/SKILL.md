---
name: oqs
description: Development guide for contributing to OpenQuantumSystems.jl, the Julia open quantum systems thesis repo. Use when working on issues, PRs, or any improvements for detrin/OpenQuantumSystems.jl.
user-invokable: true
---

# oqs — OpenQuantumSystems.jl Development Guide

OpenQuantumSystems.jl is a Julia numerical framework for simulating open quantum systems, focused on quantum biology. Built as a master thesis project.

---

## Repository setup

- **Repo:** `https://github.com/detrin/OpenQuantumSystems.jl`
- **gh CLI:** authenticated as `detrin`
- **Language:** Julia
- **Main branch:** `master`
- **Dev branch:** `devel` — commit fixes directly here, PRs merge devel → master

---

## Environment

- **Julia path:** `~/.juliaup/bin/julia` (not in PATH — always use full path)
- **Run tests:** `~/.juliaup/bin/julia --project=. -e "using Pkg; Pkg.test()"`
- **Run specific test:** `~/.juliaup/bin/julia --project=. test/test_<module>.jl`
- **Attribution:** Disabled — NO co-author lines in commits. Git hooks in `.githooks/` enforce this.
- **Git hooks:** Set with `git config core.hooksPath .githooks`

---

## Branching strategy

We are the repo owners — no fork needed. Work directly on `detrin/OpenQuantumSystems.jl`.

- **Commit directly to `devel`** for all fixes and features
- **Do NOT create feature branches** — commit straight to devel
- One commit per logical change, issue ID at the start: `#NNN Description`
- Never commit directly to `master`
- Open PRs from `devel` → `master` to merge batches of work

---

## Subagent workflow

When dispatching subagents to work on issues:

- Subagents must commit directly to `devel`, NOT create feature branches
- Prompt subagents with the OQS skill context
- Each subagent must run tests before exiting
- NO co-author lines in commits

---

## Commit message format

```
#NNN Short description of what and why

Optional longer explanation if the change is non-obvious.
```

Issue ID always first. Every commit must have an issue ID.

---

## Testing

```bash
# Run full test suite
~/.juliaup/bin/julia --project=. -e "using Pkg; Pkg.test()"

# Run specific test file
~/.juliaup/bin/julia --project=. test/test_<module>.jl
```

Test files live in `/test/` and mirror source files (e.g., `src/trace.jl` → `test/test_trace.jl`).

---

## ISSUES.md tracking

All issues are tracked in `ISSUES.md` at the repo root with three states:

- **BACKLOG** — Not started
- **IN PROGRESS** — Committed in `devel` branch
- **DONE** — Merged into `master`

When working on an issue, update BOTH:
1. The `**Status:**` line in the issue's section
2. The Status column in the summary table

---

## Key files

| File | Purpose |
|---|---|
| `src/OpenQuantumSystems.jl` | Main module, exports |
| `src/core.jl` | Shared constants and helpers (_SAFE_DIV_TOL, ELECTRONIC_GROUND/EXCITED) |
| `src/timeevolution_base.jl` | Integration helpers (_setup_delayed_integration) |
| `src/master_ansatz.jl` | QME ansatz solvers |
| `src/master_iterative.jl` | Iterative QME solvers |
| `src/evolution.jl` | Evolution operators |
| `src/trace.jl` | Bath tracing operations |
| `src/aggregateOperators.jl` | Hamiltonian and operator construction |
| `test/runtests.jl` | Test suite entry point |
| `Project.toml` | Julia project manifest |
| `ISSUES.md` | Issue tracker with status (BACKLOG/IN PROGRESS/DONE) |

---

## Known issues and technical debt

Active tracking issue: **https://github.com/detrin/OpenQuantumSystems.jl/issues/49**

Issues #50–#65 resolved in v0.3.0–v0.4.0. Remaining:

- `test/test_memory_kernel.jl` commented out in `runtests.jl`
- Docs say "under construction" — tutorial is broken (#67)
- `OpenQuantumSystems.mul!` in schroedinger.jl — intentionally kept (from QuantumOpticsBase)

---

## Issue and PR conventions

- Open a GitHub issue for every fix before or alongside the change
- Document what was found, why it matters, and how it was fixed
- Keep PRs focused — one concern per PR
- Update the tracking issue (#49) when new issues or PRs are opened
- Check off items in #49 only after the PR is merged into `master`
- **Never close an issue until its PR is merged** — issues are closed manually after merge, not via `Closes` keywords
- Use `Addresses #NNN` in PR bodies instead of `Closes #NNN`

---

## What we do NOT touch

- Scientific correctness of quantum mechanics implementations
- Algorithm logic unless a bug is clearly identified
- Any change that requires deep physics knowledge without explicit instruction
