```@meta
CurrentModule=OpenQuantumSystems
```

# OpenQuantumSystems.jl

**OpenQuantumSystems.jl** is a numerical framework written in Julia for simulating
open quantum systems, with a focus on quantum biology and molecular aggregates in
a finite basis. It is inspired by
[quantarhei](https://github.com/tmancal74/quantarhei) and
[QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

## Features

- **Aggregate construction** — Build Hamiltonians for aggregates of molecules with
  vibrational modes modelled as displaced harmonic oscillators.
- **Multiple solvers** — Liouville-von Neumann, exact QME, Redfield, ansatz-based QME,
  iterative QME, and a unified `solve()` entry point.
- **Förster theory** — Absorption/emission spectra, spectral overlap, and Förster
  excitation energy transfer rates.
- **Modified Redfield theory** — Exciton basis transformation, lineshape functions,
  and population transfer rates beyond the weak-coupling limit.
- **Dipole couplings** — Compute electronic couplings from molecular positions and
  transition dipole vectors.
- **Bath tracing** — Trace over bath degrees of freedom using Franck-Condon factors.
- **Initial states** — Thermal states and ultrafast laser excitation.
- **Post-processing** — Picture conversion, validation, scoring, and comparison tools.
- **Convenience constructors** — `setup_dimer`, `setup_trimer`, `setup_linear_chain`
  for quick prototyping.

## Installation

```julia
pkg> add OpenQuantumSystems
```

or

```julia
using Pkg; Pkg.add("OpenQuantumSystems")
```

## Quick Start

```julia
using OpenQuantumSystems

# Set up a dimer with one vibrational mode per molecule
agg = setup_dimer(E1=12500.0, E2=12700.0, J=100.0)

# Initial thermal state with laser excitation
W0 = thermal_state_composite(300.0, [0.0, 0.5, 0.5], agg)

# Time grid
tspan = get_tspan(0.0, 0.5, 500)

# Solve with QME ansatz
result = solve(W0, tspan, agg; method=:qme_ansatz)
```

See the [Solver Selection Guide](@ref) for choosing the right solver, and the
[Dimer examples](@ref) tutorial for worked examples.
