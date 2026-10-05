# Iterative QME: verification and benchmarks

Scripts behind the numbers in `docs/src/theory/memory_kernel.md` and
`docs/src/theory/iterative.md` (issue #105). They are not part of the test suite;
the fast versions of these checks live in `test/test_iterative_qme.jl`.

## Setup

```bash
cd benchmark/iterative_qme
julia --project=. -e 'using Pkg; Pkg.develop(path="../.."); Pkg.instantiate()'
```

`HierarchicalEOM` (HEOM reference) and `Plots` are only needed here, not by the package.
The scripts that need neither can also run from the package environment
(`julia --project=../.. script.jl`).

## Scripts

| script | what it checks | runtime |
|---|---|---|
| `validate_building_blocks.jl` | correlation functions, second- and fourth-order kernels and Born/Redfield solvers against the already validated finite-basis code | ~5 min |
| `check_derivation.jl` | thesis RBP equation vs exact dŵ/dt (missing quotient term), odd-order terms vanish, size of thesis vs projected correction | ~2 min |
| `order_test.jl` | disconnected-term identity; error vs Huang-Rhys factor against exact finite-bath dynamics (table in `memory_kernel.md`) | ~10 min |
| `heom_benchmark.jl` | Drude-Lorentz (infinite) bath, all methods vs HEOM for λ = 2.5…40 cm⁻¹; writes `heom_benchmark.jls` | ~1–2 h |
| `plot_heom.jl OUTDIR` | plots from `heom_benchmark.jls` (populations at λ = 10, 40 and error vs λ) | ~1 min |

Use `julia --project=. -e 'include("heom_benchmark.jl")' > heom.log 2>&1` and note that
Julia buffers stdout when redirected; the results table is printed at the end.
Memory: the earlier run with 4 Matsubara terms per bath and tier 10 was killed
(17 GB machine); with `NMATS = 1` tier 10 has 4004 ADO dimensions and is fine.
