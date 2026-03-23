# OpenQuantumSystems.jl - Open Issues

Issue tracker for code quality, refactoring, and technical debt.
Tracking issue: [#49](https://github.com/detrin/OpenQuantumSystems.jl/issues/49)

**Status definitions:**
- **BACKLOG** — Not started
- **IN PROGRESS** — Committed in `devel` branch
- **DONE** — Merged into `master`

---

## Near-Future Features

### #77 - Implement Förster theory for excitation energy transfer rates
**Status:** BACKLOG
**Severity:** Feature (near-future roadmap)

Förster theory provides excitation energy transfer (EET) rates in the weak electronic coupling regime. It is one of the two near-future items on the roadmap.

**What exists:**
- Site energies, electronic couplings, and Franck-Condon factors
- Discrete LHO bath modes with frequencies and Huang-Rhys factors (`Mode`)
- Rate constant infrastructure (`rate_constant.jl`)
- Correlation functions (`correlation_function` in trace.jl)

**What needs to be implemented:**
1. **Spectral density function** — derive from discrete LHO modes (frequencies + shifts)
2. **Line shape functions** — absorption and emission lineshapes for each molecule, computed from the spectral density and temperature
3. **Förster rate expression** — `k_DA = (2π/ℏ) |J_DA|² ∫ f_D(ω) a_A(ω) dω`, where `f_D` is donor fluorescence lineshape, `a_A` is acceptor absorption lineshape
4. **Förster rate matrix** — compute pairwise rates for all molecules in an aggregate

**Proposed files:**
- `src/spectral_density.jl` — spectral density and line shape functions (shared with #78)
- `src/forster.jl` — Förster rate calculation
- `test/test_forster.jl` — tests

**Reference:** Thesis sections on Förster theory; May & Kühn "Charge and Energy Transfer Dynamics in Molecular Systems"

---

### #78 - Implement modified Redfield theory
**Status:** BACKLOG
**Severity:** Feature (near-future roadmap)

Modified Redfield theory improves on standard Redfield by treating diagonal (population) fluctuations non-perturbatively, while keeping off-diagonal (coherence) elements perturbative. This gives better results for intermediate coupling regimes where standard Redfield breaks down.

**What exists:**
- Standard Redfield solver (`QME_sI_Redfield` in `redfield.jl`)
- Aggregate Hamiltonians in site and exciton basis
- Evolution operators and interaction picture transformations
- Spectral density infrastructure (to be added in #77)

**What needs to be implemented:**
1. **Spectral density in exciton basis** — transform site-basis spectral densities using exciton coefficients: `C_mn(ω) = Σ_k |c_mk|² |c_nk|² J_k(ω)`
2. **Modified Redfield tensor** — population transfer rates use the full lineshape (non-perturbative), coherence dephasing uses standard Redfield form
3. **Modified Redfield rate expression** — `R_mn = 2 Re ∫₀^∞ dt [g̈_mn(t) + (ġ_mm(t) - ġ_nn(t)) ġ_mn(t)] × exp(-...)` where `g(t)` is the lineshape function
4. **Solver integration** — either a standalone `QME_sI_modified_redfield` solver or an option in the existing Redfield solver

**Proposed files:**
- `src/spectral_density.jl` — shared with #77 (spectral density, lineshape functions, `g(t)` and derivatives)
- `src/modified_redfield.jl` — modified Redfield tensor and solver
- `test/test_modified_redfield.jl` — tests, compare against standard Redfield in weak-coupling limit

**Blocked by:** #77 (spectral density infrastructure)

**Reference:** Thesis sections on modified Redfield; Zhang et al. J. Chem. Phys. 108, 7763 (1998)

---

### #79 - Compute transition dipole moments from molecular coordinates
**Status:** BACKLOG
**Severity:** Feature (near-future roadmap)

Currently, electronic couplings must be set manually in the coupling matrix. For realistic simulations, dipole-dipole couplings should be computed automatically from molecular geometry.

**What needs to be implemented:**
1. **Transition dipole moment storage** — extend `Molecule` or add a wrapper to store 3D position and transition dipole vector (position, orientation, magnitude)
2. **Dipole-dipole coupling** — `J_DA = (1/4πε₀) [μ_D · μ_A / R³ - 3(μ_D · R̂)(μ_A · R̂) / R³]`
3. **Auto-populate coupling matrix** — given molecular coordinates and dipole vectors, compute the full coupling matrix for an `AggregateCore`
4. **Convenience constructor** — `AggregateCore(molecules; positions, dipoles)` that computes couplings automatically

**Proposed files:**
- `src/dipole.jl` — dipole moment types, dipole-dipole coupling calculation
- `test/test_dipole.jl` — tests with known geometries (e.g., parallel, head-to-tail dimers)

---

## Nice-to-Have Features

### #80 - Mixed state decomposition via linear programming
**Status:** BACKLOG
**Severity:** Feature (nice-to-have)

Decompose a mixed density matrix into a convex combination of pure states using linear programming. Useful for physical interpretation of simulation results.

**What needs to be implemented:**
1. **LP formulation** — express ρ = Σ_i p_i |ψ_i⟩⟨ψ_i| as a convex optimization problem
2. **Solver integration** — use JuMP.jl + HiGHS.jl (free LP solver) as optional dependencies
3. **API** — `decompose(rho; max_terms=10)` returning weights and pure states

**New dependency:** JuMP.jl, HiGHS.jl (consider making this a package extension to keep core lightweight)

**Proposed files:**
- `ext/MixedStateDecompositionExt.jl` — package extension (loaded when JuMP is available)
- `test/test_decomposition.jl` — tests

---

### #81 - Hamiltonian loading and data storage
**Status:** BACKLOG
**Severity:** Feature (nice-to-have)

Serialize and deserialize aggregates, Hamiltonians, and simulation results to disk for reproducibility and data sharing.

**What needs to be implemented:**
1. **Save/load Aggregate** — serialize the full `Aggregate` (core, tools, operators) to JLD2 format
2. **Save/load SimulationResult** — store time series and metadata
3. **Import from external formats** — read Hamiltonians from CSV/JSON for interop with other packages
4. **API** — `save_aggregate(filename, agg)`, `load_aggregate(filename)`, `save_result(filename, result)`, `load_result(filename)`

**New dependency:** JLD2.jl (consider as package extension)

**Proposed files:**
- `src/io.jl` or `ext/IOExt.jl` — serialization functions
- `test/test_io.jl` — round-trip tests

---

### #82 - Anharmonic oscillators
**Status:** BACKLOG
**Severity:** Feature (nice-to-have, physics extension)

Replace or supplement LHO modes with anharmonic potentials (e.g., Morse oscillator). This enables modeling of more realistic vibrational modes where anharmonicity matters.

**What needs to be implemented:**
1. **Anharmonic mode type** — `MorseMode` or `AnharmonicMode` with Morse potential parameters (dissociation energy, width)
2. **Numerical Franck-Condon factors** — analytic FC factors only work for LHO; anharmonic modes need numerical overlap integrals of eigenstates
3. **Eigenstate computation** — solve 1D Schrödinger equation for the Morse/anharmonic potential (DVR method or analytical Morse eigenstates)
4. **Integration with Molecule/Aggregate** — `Mode` becomes abstract, `HarmonicMode` and `MorseMode` are subtypes

**Impact:** Changes the `Mode` struct fundamentally. Needs careful design to maintain backward compatibility.

**Proposed files:**
- `src/anharmonic.jl` — anharmonic mode types and numerical FC factors
- `test/test_anharmonic.jl` — tests, verify LHO limit recovers analytic FC factors

---

### #83 - Double excited electronic states
**Status:** BACKLOG
**Severity:** Feature (nice-to-have, physics extension)

Currently the system supports ground + singly excited electronic states (N+1 states for N molecules). Double excitations (two molecules excited simultaneously) are needed for nonlinear spectroscopy (2D spectra, pump-probe).

**What needs to be implemented:**
1. **Extended electronic state space** — from `molCount+1` to `molCount+1 + molCount*(molCount-1)/2` states (ground + single + double excitations)
2. **Double-excitation Hamiltonian blocks** — couplings between single and double excitation manifolds
3. **Extended aggregate construction** — update `AggregateCore`, `AggregateTools`, `AggregateOperators` to handle the larger basis
4. **Updated tracing** — `trace_bath` and related functions must handle the new block structure

**Impact:** Architectural change affecting aggregate construction, Hamiltonians, and tracing. The `electronicIndices` function and `indicesMap` need extension. Consider a flag: `Aggregate(...; max_excitations=1)` defaulting to current behavior.

**Proposed files:**
- Modifications to `src/aggregateCore.jl`, `src/aggregateTools.jl`, `src/aggregateOperators.jl`
- `test/test_double_excitation.jl` — tests

---

### #84 - quantarhei interoperability
**Status:** BACKLOG
**Severity:** Feature (nice-to-have, ecosystem)

Enable data exchange with [quantarhei](https://github.com/tmancal74/quantarhei), the Python package that inspired OpenQuantumSystems.jl. This allows users to leverage both ecosystems.

**What needs to be implemented:**
1. **File-based interop** — define a common JSON/HDF5 format for aggregates, Hamiltonians, and spectral densities that both packages can read/write
2. **Import quantarhei Aggregate** — read quantarhei's serialized aggregate data and construct an OpenQuantumSystems `Aggregate`
3. **Export to quantarhei format** — write OpenQuantumSystems data in quantarhei-compatible format

**Alternative approach:** Use PyCall.jl for direct Python interop, but file-based exchange is simpler and more portable.

**Proposed files:**
- `src/quantarhei_io.jl` or `ext/QuantarheiExt.jl` — import/export functions
- `test/test_quantarhei_io.jl` — round-trip tests with sample data

---

## Summary of Open Issues

| # | Title | Category | Priority | Status |
|---|-------|----------|----------|--------|
| #77 | Förster theory for EET rates | Feature | Near-future | BACKLOG |
| #78 | Modified Redfield theory | Feature | Near-future | BACKLOG |
| #79 | Dipole moments from coordinates | Feature | Near-future | BACKLOG |
| #80 | Mixed state decomposition via LP | Feature | Nice-to-have | BACKLOG |
| #81 | Hamiltonian loading/storage | Feature | Nice-to-have | BACKLOG |
| #82 | Anharmonic oscillators | Feature | Nice-to-have | BACKLOG |
| #83 | Double excited states | Feature | Nice-to-have | BACKLOG |
| #84 | quantarhei interoperability | Feature | Nice-to-have | BACKLOG |

## Priority Order

**Near-future (target for next release):**
1. **#79** — Dipole moments (standalone, no dependencies)
2. **#77** — Förster theory (provides spectral density infrastructure for #78)
3. **#78** — Modified Redfield (depends on #77)

**Nice-to-have (future releases):**
4. **#81** — Hamiltonian loading/storage (practical, low effort)
5. **#83** — Double excited states (enables nonlinear spectroscopy)
6. **#82** — Anharmonic oscillators (physics extension)
7. **#80** — Mixed state decomposition (niche use case)
8. **#84** — quantarhei interop (ecosystem)
