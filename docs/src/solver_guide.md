```@meta
CurrentModule=OpenQuantumSystems
```

# Solver Selection Guide

OpenQuantumSystems.jl provides several solver families for propagating the density matrix
of molecular aggregates coupled to vibrational baths. This page helps you choose the
right solver for your problem.

## Decision Tree

```
Start here
  |
  |-- Do you need the full system+bath density matrix W(t)?
  |     |
  |     YES --> Use LvN (Liouville-von Neumann) on the full space
  |     |         |
  |     |         |-- Want the interaction picture?  --> LvN_SI
  |     |         |-- Want the Schroedinger picture? --> LvN_SS
  |     |
  |     NO  --> You need the reduced electronic density matrix rho(t)
  |               |
  |               |-- Is the system small enough for exact propagation
  |               |   (few molecules, few vibrational levels)?
  |               |     |
  |               |     YES --> Exact evolution with bath trace
  |               |     |         |
  |               |     |         |-- Interaction picture? --> Evolution_sI_exact
  |               |     |         |-- Schroedinger picture? --> Evolution_sS_exact
  |               |     |
  |               |     NO  --> You need a master equation
  |               |               |
  |               |               |-- Do you need high accuracy and can
  |               |               |   afford the cost of numerical integration?
  |               |               |     |
  |               |               |     YES --> QME exact (second-order, no approximations)
  |               |               |     |         |
  |               |               |     |         |-- Reduced density matrix?
  |               |               |     |         |     --> QME_sI_exact / QME_sS_exact
  |               |               |     |         |-- Full density matrix?
  |               |               |     |               --> QME_SI_exact / QME_SS_exact
  |               |               |     |
  |               |               |     NO  --> Choose an approximate method
  |               |               |               |
  |               |               |               |-- Is the system-bath coupling weak
  |               |               |               |   and Markovian-like?
  |               |               |               |     |
  |               |               |               |     YES --> QME_sI_Redfield
  |               |               |               |     |
  |               |               |               |     NO  --> Do you have a zeroth-order
  |               |               |               |             solution already?
  |               |               |               |               |
  |               |               |               |               YES --> QME_sI_iterative
  |               |               |               |               |
  |               |               |               |               NO  --> QME_sI_ansatz
  |               |               |               |                       (pick an ansatz
  |               |               |               |                        for the bath state)
```

## Solver Families

### 1. Liouville-von Neumann (LvN)

Solves the standard Liouville-von Neumann equation for the **full** system+bath density
matrix `W(t)`:

```math
\frac{d}{dt} W(t) = -i [H, W(t)]
```

These solvers propagate the entire Hilbert space without any bath tracing or
perturbative expansion. They serve as the exact reference when the full space
is computationally tractable.

| Function | Picture | Density matrix | Notes |
|----------|---------|---------------|-------|
| `LvN_sI` | Interaction | Reduced (traced) | Traces bath at each step |
| `LvN_sS` | Schroedinger | Reduced (traced) | Traces bath at each step |
| `LvN_SI` | Interaction | Full W(t) | No bath trace |
| `LvN_SS` | Schroedinger | Full W(t) | No bath trace |

### 2. Exact Evolution (non-ODE)

Direct matrix-exponential propagation using `U(t) = exp(-iHt)`. These compute
`W(t) = U(t) W(0) U^dag(t)` at each time point without solving a differential
equation.

| Function | Picture | Output | Notes |
|----------|---------|--------|-------|
| `Evolution_sI_exact` | Interaction | Reduced rho(t) | Reference for reduced dynamics |
| `Evolution_sS_exact` | Schroedinger | Reduced rho(t) | Reference for reduced dynamics |
| `Evolution_SI_exact` | Interaction | Full W(t) | Reference for full dynamics |
| `Evolution_SS_exact` | Schroedinger | Full W(t) | Reference for full dynamics |
| `evolution_exact` | Schroedinger | Operator array | General-purpose exact evolution |
| `evolution_approximate` | Schroedinger | Operator array | Equidistant step approximation |

### 3. QME Exact (Second-Order Quantum Master Equation)

Solves the second-order quantum master equation **without** further approximations.
The memory kernel is evaluated by numerical integration (QuadGK) at each time step,
making these solvers accurate but expensive.

```math
\frac{d}{dt} \rho^{(I)}(t) = -i [H_I^{(I)}(t), \rho(t_0)]
  - \int_0^t d\tau\; [H_I^{(I)}(t), [H_I^{(I)}(\tau), \rho^{(I)}(\tau)]]
```

| Function | Picture | Density matrix | Notes |
|----------|---------|---------------|-------|
| `QME_sI_exact` | Interaction | Reduced | Bath traced; delayed DE |
| `QME_sS_exact` | Schroedinger | Reduced | Bath traced; delayed DE |
| `QME_SI_exact` | Interaction | Full | Delayed DE on full space |
| `QME_SS_exact` | Schroedinger | Full | Delayed DE on full space |

### 4. QME Ansatz (Approximate Bath State)

Like QME exact, but the bath state at time `tau` inside the memory integral is
replaced by an ansatz. This avoids propagating the full `W(tau)` and dramatically
reduces the cost. Different ansatze offer different accuracy/cost trade-offs.

Use the unified interface:

```julia
QME_sI_ansatz(W0, tspan, agg; ansatz=:const_sch)
```

Available ansatze (passed as the `ansatz` keyword):

| Ansatz symbol | Description |
|---------------|-------------|
| `:test` | Uses exact bath state (for validation; same cost as exact) |
| `:const_int` | Bath frozen at t=0, interaction picture |
| `:const_sch` | Bath frozen at t=0, Schroedinger picture |
| `:linear_sch` | First-order linear expansion of bath, Schroedinger picture |
| `:linear2_sch` | Piecewise-linear bath evolution, Schroedinger picture |
| `:upart1_sch` | Block-diagonal unitary bath evolution, Schroedinger picture |
| `:upart1_int` | Block-diagonal unitary bath evolution, interaction picture |
| `:upart2_sch` | Full-block unitary bath evolution, Schroedinger picture |
| `:upart2_int` | Full-block unitary bath evolution, interaction picture |

### 5. Redfield

Second-order master equation with the Redfield approximation: the density matrix
inside the memory integral is evaluated at the current time `t` rather than at
the integration variable `tau`, and the correlation function depends on `t - tau`.

```julia
QME_sI_Redfield(W0, tspan, agg)
```

This is appropriate when the bath correlation time is short compared to the
system dynamics (weak-coupling, Markovian-like regime).

### 6. Iterative QME

A self-consistent iterative scheme. You first obtain a zeroth-order solution
(e.g., from QME ansatz), then feed it back to refine the bath state and
re-solve. This can systematically improve accuracy.

```julia
QME_sI_iterative(W0, rho_0_int_t, W_0_bath_t, tspan, agg)
```

Markov-level variants are available via the `method` keyword:

| Function / method | Description |
|-------------------|-------------|
| `QME_sI_iterative` (`:default`) | Full iterative correction |
| `QME_sI_iterative_markov0` (`:markov0`) | Zeroth-order Markov bath correction |
| `QME_sI_iterative_markov1` (`:markov1`) | First-order Markov bath correction |

## Comparison Table

| Solver | Picture | System | Comp. cost | Memory integral | Best for |
|--------|---------|--------|-----------|-----------------|----------|
| `LvN_SS` | Schroedinger | Full | Low | None | Small aggregates; exact reference |
| `LvN_SI` | Interaction | Full | Low | None | Small aggregates; interaction-picture reference |
| `LvN_sS` | Schroedinger | Reduced | Low | None | Quick reduced dynamics of small systems |
| `LvN_sI` | Interaction | Reduced | Low | None | Quick reduced dynamics of small systems |
| `Evolution_sI_exact` | Interaction | Reduced | Low | None | Non-ODE exact reference |
| `QME_sI_exact` | Interaction | Reduced | High | Exact (QuadGK) | Accurate benchmark; small-to-medium systems |
| `QME_sS_exact` | Schroedinger | Reduced | High | Exact (QuadGK) | Accurate benchmark in Schroedinger picture |
| `QME_SI_exact` | Interaction | Full | Very high | Exact (QuadGK) | Full-space benchmark |
| `QME_sI_ansatz` | Interaction | Reduced | Medium | Approximate bath | Production runs; tunable accuracy |
| `QME_sI_Redfield` | Interaction | Reduced | Medium | Redfield approx. | Weak coupling; Markovian baths |
| `QME_sI_iterative` | Interaction | Reduced | High | Iterative refinement | Systematic improvement over ansatz |

**Cost key:** Low = single matrix exponential per step; Medium = numerical integration
of memory kernel with approximate bath; High = numerical integration with exact or
iterative bath; Very high = full-space delayed differential equation.

## Common Usage Example

A typical workflow for a molecular dimer: set up the aggregate, choose a solver,
and convert results to the Schroedinger picture if needed.

```julia
using OpenQuantumSystems

# --- Define the aggregate ---
HR = 0.01
shift = sqrt(2.0 * HR)
mode = Mode(300.0, shift)
mols = [
    Molecule([mode], 5, [12500.0, 12750.0]),
    Molecule([mode], 5, [12500.0, 12800.0]),
]
agg = Aggregate(mols)
agg.coupling[2, 3] = 100.0
agg.coupling[3, 2] = 100.0
setup_aggregate!(agg)

# --- Initial state ---
T = 300.0  # temperature in K
W0 = thermal_state_composite(T, [0.0, 0.5, 0.5], agg.operators.Ham, agg.tools)

# --- Time grid ---
t_max = 0.5
tspan = get_tspan(0.0, t_max, 500)

# --- Solve with QME ansatz (good default for production) ---
tspan_out, rho_int_t = QME_sI_ansatz(W0, tspan, agg; ansatz=:const_sch)

# --- Convert from interaction picture to Schroedinger picture ---
rho_sch_t = interaction_pic_to_schroedinger_pic(rho_int_t, agg.operators.Ham_0, tspan)
```

## Tips

- **Start with `Evolution_sI_exact` or `LvN_sI`** to establish a reference solution
  for your system size, then switch to an approximate solver for production.
- **Interaction picture** (suffix `_sI` or `_SI`) is generally more numerically stable
  for oscillatory dynamics because the fast free-evolution phase is factored out.
- **Tolerance tuning:** The `int_reltol` and `int_abstol` parameters control the
  accuracy of the memory-kernel quadrature. Tighten them for benchmarks, loosen for
  speed.
- **Ansatz choice:** Start with `:const_sch` (cheapest nontrivial ansatz). If accuracy
  is insufficient, move to `:upart1_sch` or `:upart2_sch`, or switch to the iterative
  solver.
