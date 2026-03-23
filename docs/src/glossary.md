# Naming Conventions and Glossary

This page explains the abbreviations and naming patterns used throughout the
OpenQuantumSystems.jl API. Understanding these conventions makes it much easier
to find and use the right function.

## Density matrix variables: `W` vs `rho`

| Variable | Meaning |
|----------|---------|
| `rho`    | Reduced (system-only) density matrix, obtained after tracing out bath degrees of freedom. |
| `W`      | Full density matrix of the composite system (electronic + vibrational/bath). |
| `W_bath` | Bath-only density matrix (partial trace over electronic states). |

## System-size and picture prefixes

Function names encode **two choices** as a two-character prefix:

| Character | Position | Meaning |
|-----------|----------|---------|
| `s` (lowercase) | 1st | **Reduced system** -- the function works with or returns the traced-over (small) density matrix `rho`. |
| `S` (uppercase) | 1st | **Full system** -- the function works with or returns the composite density matrix `W`. |
| `I` (uppercase) | 2nd | **Interaction picture** -- operators and states are in the interaction picture. |
| `S` (uppercase) | 2nd | **Schroedinger picture** -- operators and states are in the Schroedinger picture. |

The four combinations are:

| Prefix | System size | Picture | Example |
|--------|-------------|---------|---------|
| `sI` | Reduced | Interaction | `QME_sI_exact`, `LvN_sI` |
| `sS` | Reduced | Schroedinger | `QME_sS_exact`, `LvN_sS` |
| `SI` | Full | Interaction | `QME_SI_exact`, `LvN_SI` |
| `SS` | Full | Schroedinger | `QME_SS_exact`, `LvN_SS` |

## Bath-state ansatz suffixes

Functions in `master_ansatz.jl` carry additional suffixes that specify how the
bath state is approximated during time integration:

| Suffix | Meaning |
|--------|---------|
| `const_int` | Constant bath state, interaction picture. |
| `const_sch` | Constant bath state, Schroedinger picture. |
| `linear_sch` | Linear interpolation of bath state, Schroedinger picture. |
| `linear2_sch` | Second-order linear interpolation, Schroedinger picture. |
| `upart1_sch` | Partial evolution operator (variant 1), Schroedinger-picture bath. |
| `upart1_int` | Partial evolution operator (variant 1), interaction-picture bath. |
| `upart2_sch` | Partial evolution operator (variant 2), Schroedinger-picture bath. |
| `upart2_int` | Partial evolution operator (variant 2), interaction-picture bath. |

## Abbreviation glossary

| Abbreviation | Meaning |
|--------------|---------|
| `QME` | Quantum Master Equation |
| `LvN` | Liouville-von Neumann (equation) |
| `FC` | Franck-Condon (factors) |
| `Ham` | Hamiltonian |
| `agg` | Aggregate (molecular aggregate) |
| `mol` | Molecule |
| `vib` | Vibrational (mode / index) |
| `el` | Electronic (state / index) |
| `ad` | Adjoint assembly -- constructs composite `W` from reduced `rho` and `W_bath` |
| `Nvib` | Number of vibrational levels |
| `intp` | Interpolation |
| `op` | Operator |

## Aggregate Hamiltonian components

`getAggHam*` functions build different parts of the aggregate Hamiltonian:

| Function | Returns |
|----------|---------|
| `getAggHamSystemSmall` | Electronic-only system Hamiltonian (small basis). |
| `getAggHamSystemBig` | System Hamiltonian in the full electronic+vibrational basis. |
| `getAggHamBathSmall` | Bath (vibrational) Hamiltonian in the small basis. |
| `getAggHamBathBig` | Bath Hamiltonian in the full basis. |
| `getAggHamSystemBath` | System-bath coupling Hamiltonian. |
| `getAggHamInteraction` | Interaction Hamiltonian (off-diagonal electronic couplings). |
| `getAggHamiltonian` | Total Hamiltonian of the aggregate. |

## Function name decoding examples

Here are some fully decoded function names:

- **`QME_sI_exact`** -- Quantum Master Equation, reduced system, interaction picture, exact (non-iterative) solver.
- **`QME_sI_ansatz_upart2_sch`** -- QME, reduced system, interaction picture, ansatz method with partial-U variant 2 and Schroedinger-picture bath state.
- **`QME_sI_iterative_markov0`** -- QME, reduced system, interaction picture, iterative solver with zeroth-order Markov approximation.
- **`QME_sI_Redfield`** -- QME, reduced system, interaction picture, Redfield (Markovian) approximation.
- **`LvN_SS`** -- Liouville-von Neumann equation, full system, Schroedinger picture.
- **`Evolution_SI_exact`** -- Time evolution, full system, interaction picture, exact propagation.
- **`K_aabb_W_bath_intp`** -- Rate constant with `aa`/`bb` electronic index pattern, computed from full `W` bath state with interpolation.
- **`M_abcd_W_bath_intp`** -- Memory kernel element with `abcd` electronic index pattern, computed from full `W` bath state with interpolation.
- **`getAggHamBathBig`** -- Get aggregate bath Hamiltonian in the full (big) basis.
- **`getFCproduct`** -- Get Franck-Condon factor product.
