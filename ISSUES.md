# OpenQuantumSystems.jl - Open Issues

Issue tracker for code quality, refactoring, and technical debt.
Tracking issue: [#49](https://github.com/detrin/OpenQuantumSystems.jl/issues/49)

**Status definitions:**
- **BACKLOG** — Not started
- **IN PROGRESS** — Committed in `devel` branch
- **DONE** — Merged into `master`

---

## Near-Future Features

### #78 - Implement Förster theory for excitation energy transfer rates
**Status:** DONE
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
- `src/spectral_density.jl` — spectral density and line shape functions (shared with #79)
- `src/forster.jl` — Förster rate calculation
- `test/test_forster.jl` — tests

**Reference:** Thesis sections on Förster theory; May & Kühn "Charge and Energy Transfer Dynamics in Molecular Systems"

---

### #79 - Implement modified Redfield theory
**Status:** DONE
**Severity:** Feature (near-future roadmap)

Modified Redfield theory improves on standard Redfield by treating diagonal (population) fluctuations non-perturbatively, while keeping off-diagonal (coherence) elements perturbative. This gives better results for intermediate coupling regimes where standard Redfield breaks down.

**What exists:**
- Standard Redfield solver (`QME_sI_Redfield` in `redfield.jl`)
- Aggregate Hamiltonians in site and exciton basis
- Evolution operators and interaction picture transformations
- Spectral density infrastructure (to be added in #78)

**What needs to be implemented:**
1. **Spectral density in exciton basis** — transform site-basis spectral densities using exciton coefficients: `C_mn(ω) = Σ_k |c_mk|² |c_nk|² J_k(ω)`
2. **Modified Redfield tensor** — population transfer rates use the full lineshape (non-perturbative), coherence dephasing uses standard Redfield form
3. **Modified Redfield rate expression** — `R_mn = 2 Re ∫₀^∞ dt [g̈_mn(t) + (ġ_mm(t) - ġ_nn(t)) ġ_mn(t)] × exp(-...)` where `g(t)` is the lineshape function
4. **Solver integration** — either a standalone `QME_sI_modified_redfield` solver or an option in the existing Redfield solver

**Proposed files:**
- `src/spectral_density.jl` — shared with #78 (spectral density, lineshape functions, `g(t)` and derivatives)
- `src/modified_redfield.jl` — modified Redfield tensor and solver
- `test/test_modified_redfield.jl` — tests, compare against standard Redfield in weak-coupling limit

**Blocked by:** #78 (spectral density infrastructure)

**Reference:** Thesis sections on modified Redfield; Zhang et al. J. Chem. Phys. 108, 7763 (1998)

---

### #80 - Compute transition dipole moments from molecular coordinates
**Status:** DONE
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

## Bug Fixes

### #86 - Fix missing rho factor in default iterative QME bath correction
**Status:** IN PROGRESS
**Severity:** Bug

In `master_iterative.jl`, the default I ansatz (`method=:default`) had a line that overwrote the correctly computed `ad(rho_t2, W_bath_t2, ...)` with just `W_bath_t2`, dropping the ρ_{cd}(t₂) weighting required by thesis Eq. 3.28. The I.M1 and I.M2 variants were not affected (the overwrite was already commented out there).

**Fix:** Removed the erroneous `tmp1[:, :] = W_bath_t2[:, :]` on line 43 of `W_abcd_1_bath_core`.

---

### #87 - Fix docstring labels and type parameters in memory_kernel.jl and master_exact.jl
**Status:** IN PROGRESS
**Severity:** Bug (minor)

Two minor issues found during thesis verification:
1. All four `MemoryKernel_*_traced` functions in `memory_kernel.jl` had the same LaTeX label `M_1` instead of `M_1`..`M_4`.
2. `QME_sS_exact` in `master_exact.jl` declared unused type parameters `U` and `V`.

**Fix:** Corrected docstring labels to `M_2`, `M_3`, `M_4` and removed unused type parameters.

---

## Corrected Memory Kernel (Thesis Section 3.7)

### #88 - Zeroth-order memory kernel via correlation functions
**Status:** DONE
**Severity:** Feature (thesis core contribution)
**Blocked by:** None

Express the memory kernel M_abcd(t₁, t₂; ŵ⁰⁽ᴵ⁾) using first-order correlation functions C_nn(t) instead of explicit operator products. This enables use with infinite baths where explicit bath operators cannot be constructed.

**Thesis reference:** Section 3.7.1, Equations 3.53–3.55

**What needs to be implemented:**
1. **Analytic first-order correlation function** `C_nn(t₁, t₂)` — closed-form from Eq. 1.36 using Bose-Einstein distribution, mode frequencies, and Huang-Rhys factors
2. **Exciton transformation coefficients** `K^nm_ab,cd` from Eq. 3.62
3. **Multi-frequency phase factor** `Ω_ab,cd,ef,gh` from Eq. 3.61
4. **Zeroth-order memory kernel** `M_abcd(t₁, t₂; ŵ⁰)` from Eq. 3.55

**Proposed files:**
- `src/evolution/correlation.jl` — analytic correlation function C_nn(t₁,t₂) and helpers (K, Ω)
- `src/evolution/corrected_memory_kernel.jl` — corrected memory kernel
- `test/test_corrected_memory_kernel.jl` — tests, verify against existing `MemoryKernel_traced` for finite systems

---

### #89 - First-order corrected memory kernel
**Status:** IN PROGRESS
**Severity:** Feature (thesis core contribution)
**Blocked by:** #88

Implement the first-order correction to the memory kernel M_abcd(t₁, t₂; ŵ¹⁽ᴵ⁾) which accounts for the non-equilibrium evolution of the relative bath part. This is the central theoretical result of the thesis.

**Thesis reference:** Section 3.7.1, Equations 3.56–3.66

**What needs to be implemented:**
1. **General corrected memory kernel** `M_abcd(t₁, t₂; ŵ¹)` from Eq. 3.59 — 16 terms (M_1a..M_4d) involving double integrals over t₃, t₄ with second-order correlation functions decomposed into first-order ones
2. **Population-only version** `M_aabb(t₁, t₂; ŵ¹)` from Eq. 3.66 — simplified form for population transfer rates
3. **Markov approximation** ρ^{1,(I)}(t₄) ≈ ρ^{1,(I)}(t₂) for numerical feasibility (thesis page 75)
4. **QME-RDM iterative solver** using corrected memory kernel — the iterative scheme ŵ⁰ → ρ¹ → ρ² → ... (Eq. 3.65)

**Proposed files:**
- `src/evolution/corrected_memory_kernel.jl` — extend with first-order correction
- `test/test_corrected_memory_kernel.jl` — extend tests

---

### #90 - High-temperature limit of corrected memory kernel
**Status:** BACKLOG
**Severity:** Feature (thesis core contribution)
**Blocked by:** #89

In the high-temperature limit, correlation functions become real (Eq. 1.55), reducing the 16 terms of the corrected memory kernel to 5 unique groups (Eq. 3.75). This simplification makes the corrected kernel more practical for numerical use.

**Thesis reference:** Section 3.7.3, Equations 3.75 and surrounding discussion

**What needs to be implemented:**
1. **High-temperature correlation function** — real-valued C_nn(t₁,t₂) from Eq. 1.55
2. **Grouped memory kernel terms** — 5 unique second-order correlation function groups from Eq. 3.75
3. **Simplified corrected memory kernel** — exploit symmetry C_{n,m}(t₁,t₂) = C_{cd,ab}(t₂,t₁) (Eq. 1.56) to halve the computation

**Proposed files:**
- `src/evolution/corrected_memory_kernel.jl` — high-temperature specializations
- `test/test_corrected_memory_kernel.jl` — extend tests, verify high-T limit matches general case

---

## Nice-to-Have Features

### #82 - Hamiltonian loading and data storage
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

### #83 - Anharmonic oscillators
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

### #84 - Double excited electronic states
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

### #85 - quantarhei interoperability
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
| #78 | Förster theory for EET rates | Feature | Near-future | DONE |
| #79 | Modified Redfield theory | Feature | Near-future | DONE |
| #80 | Dipole moments from coordinates | Feature | Near-future | DONE |
| #82 | Hamiltonian loading/storage | Feature | Nice-to-have | BACKLOG |
| #83 | Anharmonic oscillators | Feature | Nice-to-have | BACKLOG |
| #84 | Double excited states | Feature | Nice-to-have | BACKLOG |
| #85 | quantarhei interoperability | Feature | Nice-to-have | BACKLOG |
| #86 | Fix missing rho factor in iterative QME | Bug | Critical | IN PROGRESS |
| #87 | Fix docstring labels and type params | Bug | Minor | IN PROGRESS |
| #88 | Zeroth-order memory kernel via correlation functions | Feature | Near-future | DONE |
| #89 | First-order corrected memory kernel | Feature | Near-future | IN PROGRESS |
| #90 | High-temperature limit of corrected memory kernel | Feature | Near-future | BACKLOG |

## Priority Order

**Near-future (target for next release):**
1. **#80** — Dipole moments (standalone, no dependencies)
2. **#78** — Förster theory (provides spectral density infrastructure for #79)
3. **#79** — Modified Redfield (depends on #78)

**Corrected memory kernel (thesis Section 3.7):**
4. **#88** — Zeroth-order memory kernel via correlation functions (foundation for #89)
5. **#89** — First-order corrected memory kernel (depends on #88)
6. **#90** — High-temperature limit (depends on #89)

**Nice-to-have (future releases):**
7. **#82** — Hamiltonian loading/storage (practical, low effort)
8. **#84** — Double excited states (enables nonlinear spectroscopy)
9. **#83** — Anharmonic oscillators (physics extension)
10. **#85** — quantarhei interop (ecosystem)
