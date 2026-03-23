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

---

## Branching strategy

We are the repo owners — no fork needed. Work directly on `detrin/OpenQuantumSystems.jl`.

- Feature branches off `master` with descriptive names
- One commit per logical change, issue ID at the start: `#NNN Description`
- Never commit directly to `master`
- Open PRs against `master` on the same repo

---

## PR workflow

1. Branch off `master` with a descriptive name
2. One commit per logical change: `#NNN Short description`
3. Test locally before pushing — never push a known red state
4. Push branch and open PR against `master`
5. Link issues in PR body with `Closes #NNN`
6. Document findings in the relevant GitHub issues before and after fixing

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
julia --project=. -e "using Pkg; Pkg.test()"

# Run specific test file
julia --project=. test/test_<module>.jl
```

Test files live in `/test/` and mirror source files (e.g., `src/trace.jl` → `test/test_trace.jl`).

---

## Key files

| File | Purpose |
|---|---|
| `src/OpenQuantumSystems.jl` | Main module, exports |
| `src/master_ansatz.jl` | QME ansatz solvers |
| `src/master_iterative.jl` | Iterative QME solvers |
| `src/evolution.jl` | Evolution operators |
| `src/trace.jl` | Bath tracing operations |
| `src/aggregateOperators.jl` | Hamiltonian and operator construction |
| `test/runtests.jl` | Test suite entry point |
| `Project.toml` | Julia project manifest |

---

## Known issues and technical debt

Active tracking issue: **https://github.com/detrin/OpenQuantumSystems.jl/issues/49**

Issues #50–#58 were resolved in v0.3.0 (PR #59). Remaining debt:

- `test/test_memory_kernel.jl` commented out in `runtests.jl`
- Docs say "under construction" — only one tutorial exists
- 13 unresolved TODO comments in source

---

## Issue and PR conventions

- Open a GitHub issue for every fix before or alongside the change
- Document what was found, why it matters, and how it was fixed
- Keep PRs focused — one concern per PR
- Update the tracking issue (#49) when new issues or PRs are opened
- Check off items in #49 only after the PR is merged into `master`
- **Never close an issue until its PR is merged**

---

## What we do NOT touch

- Scientific correctness of quantum mechanics implementations
- Algorithm logic unless a bug is clearly identified
- Any change that requires deep physics knowledge without explicit instruction
