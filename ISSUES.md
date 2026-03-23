# OpenQuantumSystems.jl - Open Issues

Issue tracker for code quality, refactoring, and technical debt.
Tracking issue: [#49](https://github.com/detrin/OpenQuantumSystems.jl/issues/49)

---

## Existing Issues

### #60 - Fix unsafe divisions in trace.jl and rate_constant.jl
**Status:** Open
**Severity:** Correctness risk

Issue #52 covers risky division in `master_iterative.jl:340`, but similar unsafe division patterns exist in other files:

- `trace.jl:297` - `W[...] / rho[el1_p, el2_p]` -- division without explicit zero check
- `trace.jl:305` - Same pattern -- denominator could be zero
- `rate_constant.jl:38` - `M_tr[a, a] / rho_t[a, a]` -- diagonal element as denominator
- `rate_constant.jl:90` - `M_tr[a, b] / rho_t[a, b]` -- off-diagonal element as denominator

These can silently produce `Inf` or `NaN` when density matrix elements are zero.

**Proposed fix:** Add explicit zero-checks or a safe division helper:
```julia
safe_div(a, b; tol=1e-10) = abs(b) < tol ? zero(a) : a / b
```

Related: #52

---

## Issues To Create

### Extract shared integration setup boilerplate into helper function
**Severity:** Moderate (DRY -- ~150 duplicated lines)

The integration setup pattern is duplicated across 4 files:
- `master_ansatz.jl:207-255`
- `master_exact.jl:31-78`
- `master_iterative.jl:99-196`
- `redfield.jl:4-52`

Common boilerplate repeated each time:
```julia
history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
rho0 = trace_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=...)
tmp1 = copy(W0.data)
tmp2 = copy(W0.data)
tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
```

**Proposed fix:** Extract a shared helper (e.g. `_setup_integration`) that returns the common objects, reducing each call site to 2-3 lines.

Related: #50, #51

---

### Refactor deeply nested getAggHamInteraction in aggregateOperators.jl
**Severity:** Moderate (readability/maintainability)

`getAggHamInteraction` in `aggregateOperators.jl:206-312` is 107 lines long with 5 levels of nesting:

```
if vib_basis == :ground_ground
    for mol_i = 1:aggCore.molCount
        if mol_i == elOrder1 - 1
            for mode_i = 1:length(mol.modes)
                if diff_vib == 1
                    # 5 levels deep
```

**Proposed fix:**
- Extract inner loops into named helper functions
- Use early returns / `continue` to reduce nesting
- Consider splitting the vib_basis branches into separate methods via dispatch

---

### Remove Union{T, Nothing} fields from Aggregate struct
**Severity:** Moderate (SOLID -- Liskov Substitution)

The `Aggregate` struct (`aggregate.jl:6-16`) uses `Union{T, Nothing}` for its fields:

```julia
mutable struct Aggregate
    core::Union{AggregateCore, Nothing}
    tools::Union{AggregateTools, Nothing}
    operators::Union{AggregateOperators, Nothing}
end
```

Every function that accepts `Aggregate` must check for `nothing` before accessing fields, adding defensive code throughout the codebase.

**Proposed fix:**
- Make fields non-nullable (require all at construction time), OR
- Use a builder pattern where partial aggregates are a different type
- Remove scattered `nothing` checks once fields are guaranteed present

---

### Replace hardcoded 2-level system assumptions in molecules.jl
**Severity:** Moderate (extensibility)

`molecules.jl` and `aggregateOperators.jl` hardcode magic numbers assuming a 2-level system (ground + excited):

- `molecules.jl:217-223` - Electronic state values `1` (ground) and `2` (excited) used directly
- `aggregateOperators.jl:40-41` - `E_agg = zeros(Float64, (2, aggCore.molCount))` -- hardcoded 2 levels
- `aggregateOperators.jl:78,82` - `E_agg[1, :]` and `E_agg[2, :]` -- magic indices

This prevents extension to 3-level or N-level systems.

**Proposed fix:**
- Define constants: `const GROUND = 1; const EXCITED = 2`
- Better: parameterize on number of electronic levels
- Use `mol.nlevels` instead of hardcoded `2`

---

### Remove hardcoded OpenQuantumSystems module self-references
**Severity:** Minor (code smell)

Several files explicitly reference the module name instead of using direct function calls:

- `aggregateTools.jl:101,102,168,172` - `OpenQuantumSystems.elIndOrder(elind1)`
- `trace.jl:27,31,129,133` - `OpenQuantumSystems.elIndOrder(...)`

These should be direct calls since they are within the same module.

**Proposed fix:**
- Replace `OpenQuantumSystems.elIndOrder(x)` with `elIndOrder(x)` throughout
- If scope issues exist, fix the underlying export/include order

---

### Standardize naming conventions across codebase
**Severity:** Minor (breaking change -- plan for major version)

The codebase mixes camelCase and snake_case inconsistently:

**camelCase (Julia non-standard):**
- `getAggHamSystemSmall`, `getMolStateEnergy`, `getMolStateFC`
- `elIndOrder`, `electronicIndices`, `vibrationalIndices`

**snake_case (Julia standard):**
- `trace_bath`, `vib_basis`, `trace_bath_ground_excited`
- `get_rho_bath`, `operator_recast`

**Proposed fix:**
- Rename exported functions to snake_case (e.g. `get_agg_ham_system_small`)
- Keep old names as deprecated aliases for one release cycle
- Update internal variable names to be consistent

---

## Summary of All Open Issues

| # | Title | Category | Severity |
|---|-------|----------|----------|
| #49 | Code quality tracking | Meta | -- |
| #50 | Refactor master_ansatz.jl: eliminate 8-way duplication | DRY | Critical |
| #51 | Refactor master_iterative.jl: consolidate duplicated solvers | DRY | Critical |
| #52 | Fix risky division in master_iterative.jl:340 | Correctness | Critical |
| #53 | Replace typeof(rho) <: Operator anti-pattern with dispatch | SOLID | Moderate |
| #54 | Replace bare parameter tuples with NamedTuple or struct | SOLID | Moderate |
| #55 | Replace magic symbol strings with @enum | SOLID | Moderate |
| #56 | Fix potential nothing return in trace.jl:71-79 | Correctness | Moderate |
| #57 | Improve type annotations throughout | Quality | Moderate |
| #58 | Clean up technical debt: TODOs, dead code, commented-out tests | Debt | Moderate |
| #60 | Fix unsafe divisions in trace.jl and rate_constant.jl | Correctness | Moderate |
| TBD | Extract shared integration setup boilerplate | DRY | Moderate |
| TBD | Refactor deeply nested aggregateOperators.jl | Quality | Moderate |
| TBD | Remove Union{T, Nothing} from Aggregate struct | SOLID | Moderate |
| TBD | Replace hardcoded 2-level system assumptions | Extensibility | Moderate |
| TBD | Remove hardcoded module self-references | Code smell | Minor |
| TBD | Standardize naming conventions | Convention | Minor |

## Priority Order (Internal Quality)

1. **#52, #60, #56** -- Correctness risks (divisions, nothing returns)
2. **#53** -- typeof anti-pattern (enables better dispatch)
3. **#50, #51** + integration boilerplate -- Biggest DRY wins (~230+ lines)
4. **#55** -- Magic symbols to enums (unlocks Open/Closed)
5. **#58** -- Tech debt cleanup (commented tests, dead code)
6. **#54** + Union{Nothing} + aggregateOperators nesting -- Structural quality
7. **#57** + 2-level assumptions -- Type safety and extensibility
8. Module self-references + naming conventions -- Polish

---

## User Experience & Developer Experience Issues

The 17 issues above fix **internal** code quality. The issues below address the **external** experience -- how physicists and developers actually interact with the library.

### UX-1: Tutorial uses old API and does not work
**Severity:** Critical (blocks all new users)

The only tutorial (`docs/src/tutorials/dimer.md`) uses an API that no longer exists:

```julia
# Tutorial shows (BROKEN):
aggInds, vibindices, aggIndLen, basis, FCFact, FCProd, Ham, Ham_0, Ham_I =
    setupAggregate(agg; verbose=false)
trace_bath(W.data, agg, FCProd, aggInds, vibindices)

# Current API:
agg = setupAggregate(aggCore)  # returns Aggregate
trace_bath(W, agg)
```

A new user following the docs will immediately get errors.

**Proposed fix:** Rewrite the dimer tutorial to use the current API. Add a quickstart section to the docs landing page.

---

### UX-2: No SimulationResult type -- inconsistent return values
**Severity:** High (API confusion, overload explosion)

Every solver returns something different:

| Function | Returns |
|----------|---------|
| `thermal_state(...)` | `Operator` |
| `ultrafast_laser_excitation(...)` | `(W0, rho0, W0_bath)` -- 3 values |
| `LvN_sI(...)` | `(tspan, rho_t)` |
| `QME_sI_iterative(...)` | `(tspan, rho_t, W_1_bath_t)` -- 3 values |
| `evolution_approximate(...)` | `(tspan, W_t)` |

And `rho_t` is sometimes `Vector{Operator}`, sometimes `Array{ComplexF64, 3}`, causing the 4-overload explosion in `scoring.jl` (12 methods for 3 functions).

**Proposed fix:** Introduce a `SimulationResult` type:

```julia
struct SimulationResult
    tspan::Vector{Float64}
    states::Vector{Operator}
    method::Symbol
    aggregate::Aggregate
end

populations(r::SimulationResult) = [real(diag(s.data)) for s in r.states]
Base.getindex(r::SimulationResult, i) = r.states[i]
```

This eliminates the `OperatorVector` / `OperatorVectorArray` duality and collapses the 4-overload pattern in scoring.jl to single methods.

---

### UX-3: Coupling matrix includes ground state index -- confusing for physicists
**Severity:** High (domain confusion)

The coupling matrix is `(molCount+1) x (molCount+1)` because index 1 = ground state:

```julia
aggCore = AggregateCore([mol1, mol2])
aggCore.coupling[2, 3] = 100   # Why 2,3? Because 1 = ground state!
aggCore.coupling[3, 2] = 100
```

A physicist expects a 2x2 coupling matrix for a dimer. The +1 offset is an implementation detail that leaks into the API.

**Proposed fix:** Accept a molecule-only coupling matrix and handle the ground state offset internally:

```julia
# User writes:
AggregateCore([mol1, mol2]; coupling=[0 100; 100 0])

# Library internally pads to (molCount+1) x (molCount+1)
```

---

### UX-4: Function names are cryptic abbreviations
**Severity:** Moderate (discoverability)

| Name | Decoded meaning |
|------|----------------|
| `QME_sI_ansatz_upart2_sch` | QME, reduced system, interaction picture, partial-U ansatz v2, Schrodinger bath |
| `LvN_sI` vs `LvN_SI` | lowercase `s` = reduced; uppercase `S` = full. One character changes the physics. |
| `ad` | Not adjoint -- inverse of trace_bath. 2 chars for a core operation. |
| `W0` vs `rho0` | `W` = full density matrix; `rho` = reduced. Undocumented convention. |

Compare with QuantumOptics.jl's `timeevolution.master` -- immediately clear.

**Proposed fix:**
- Document the `s/S` and `I/S` naming conventions prominently
- Rename `ad` to `reconstruct_from_bath` or `embed_bath` (breaking, defer to major version)
- Add a glossary to the docs: W, rho, s, S, I, FC, etc.

---

### UX-5: No unified solver entry point
**Severity:** Moderate (API discoverability)

Users must discover and remember 8+ top-level solver functions:
`LvN_sI`, `LvN_sS`, `LvN_SI`, `LvN_SS`, `QME_sI_exact`, `QME_sS_exact`, `QME_sI_ansatz`, `QME_sI_Redfield`, `QME_sI_iterative`, ...

**Proposed fix:** Add a unified entry point alongside the existing functions:

```julia
solve(W0, tspan, agg; method=:qme_ansatz, kwargs...) -> SimulationResult
```

The individual functions remain for advanced use, but `solve` provides discoverability.

---

### UX-6: No physical validation
**Severity:** Moderate (silent numerical errors)

If a solver produces a non-physical state (trace != 1, not positive semidefinite), the user gets no warning. Silent `Inf`/`NaN` propagation is possible.

**Proposed fix:** Add opt-in validation:

```julia
function validate(result::SimulationResult; tol=1e-6)
    for (i, rho) in enumerate(result.states)
        tr_val = real(tr(rho))
        abs(tr_val - 1.0) < tol || @warn "Trace = $tr_val at t=$(result.tspan[i])"
    end
end
```

---

### UX-7: Unit conversion is fragile and duplicated
**Severity:** Minor

`convert_units(E; from="1/cm", to="1/fs")` supports only one conversion path. Constants `c`, `pi` are redefined locally instead of using Julia's `pi`. The same conversion is duplicated in `tspan_cm_to_fs` in `postprocessing.jl`.

**Proposed fix:**
- Use Julia's `pi` constant
- Consolidate `tspan_cm_to_fs` into `convert_units`
- Consider `Unitful.jl` integration for a future version

---

### UX-8: No solver selection guide in docs
**Severity:** Minor (onboarding)

A new user has no guidance on which solver to use. Proposed decision tree for documentation:

```
Which solver should I use?
+-- Need exact dynamics? -> evolution_exact (expensive, reference)
+-- Small system, short time? -> LvN (Liouville-von Neumann)
+-- Need open system dynamics?
|   +-- Weak coupling? -> Redfield
|   +-- Need memory effects?
|   |   +-- Known bath ansatz? -> QME_ansatz (choose variant)
|   |   +-- Iterative correction? -> QME_iterative
|   +-- Just the kernel? -> MemoryKernel
+-- Rate constants only? -> M_aabb / K_aabb
```

---

### UX-9: Convenience constructors for common systems
**Severity:** Low (quality of life)

Common systems (dimer, trimer, linear chain) require 15+ lines of boilerplate. A convenience layer would lower the barrier to entry:

```julia
dimer = setup_dimer(
    E = [12750., 12800.],
    J = 100.0,
    mode = Mode(omega=300., hr_factor=0.01),
    nvib = 5
)
```

These are thin wrappers -- they don't change the core, just provide a friendlier entry point.

---

### UX-10: Adopt Julia-idiomatic patterns to replace ad-hoc dispatch and parameter passing
**Severity:** High (architectural -- enables all future extensibility)

The codebase uses OOP-influenced patterns (symbol dispatch tables, bare parameter tuples, Union{Nothing} builders) where Julia-native patterns would be simpler and more extensible.

#### Pattern 1: Type dispatch instead of Symbol dispatch

**Current** (closed for extension):
```julia
if vib_basis == :ground_ground
    trace_bath_ground_ground(...)
elseif vib_basis == :ground_excited
    trace_bath_ground_excited(...)
end
```

**Proposed** (open for extension):
```julia
struct GroundGround end
struct GroundExcited end
trace_bath(W, agg, ::GroundGround) = ...
trace_bath(W, agg, ::GroundExcited) = ...
# Adding a new basis = new type + new method. No existing code modified.
```

Applies to: `vib_basis` dispatch in trace.jl/aggregateTools.jl/aggregateOperators.jl, ansatz dispatch in master_ansatz.jl, method dispatch in master_iterative.jl.

#### Pattern 2: Abstract type hierarchy for ansatzes

**Current** (`_bath_state_fn` symbol lookup table):
```julia
function _bath_state_fn(ansatz::Symbol)
    ansatz == :test        && return _bath_state_test
    ansatz == :const_int   && return _bath_state_const_int
    # ... 9 branches
end
```

**Proposed** (dispatch on abstract type):
```julia
abstract type AbstractAnsatz end
struct TestAnsatz <: AbstractAnsatz end
struct ConstIntAnsatz <: AbstractAnsatz end

bath_state(::TestAnsatz, s, p, rho) = ...
bath_state(::ConstIntAnsatz, s, p, rho) = ...

# Generic kernel dispatches automatically:
kernel(t, s, ansatz::AbstractAnsatz, ...) = ...
```

New ansatz = new type + one method definition. No modification of existing code.

#### Pattern 3: Callable structs (functors) for solver state

**Current** (bare tuples, anonymous closures):
```julia
p = (aggCore=agg.core, aggTools=agg.tools, aggOperators=agg.operators,
     W0=W0, W0_bath=W0_bath, t_mk_bath_step=t_mk_bath_step, elementtype=eltype(W0))
dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz(
    t, rho, drho, history_fun, tmp1, tmp2, p, int_reltol, int_abstol, ansatz)
```

**Proposed** (callable struct):
```julia
struct QMESolver{A<:AbstractAnsatz}
    agg::Aggregate
    ansatz::A
    W0::Operator
    W0_bath::Operator
    tmp1::Matrix{ComplexF64}
    tmp2::Matrix{ComplexF64}
    int_reltol::Float64
    int_abstol::Float64
end

(s::QMESolver)(t, rho, drho, history_fun) = dQME_sI(t, rho, drho, history_fun, s)
```

Eliminates: anonymous closures capturing variables, bare tuple unpacking (`aggCore, aggTools, aggOperators, W0, elementtype = p`), and the 9-element tuples in rate_constant.jl. The solver **is** its parameters.

#### Pattern 4: Pipeline constructors instead of Union{Nothing} builder

**Current** (broken builder):
```julia
agg = Aggregate(aggCore, nothing, nothing)
# ... user must call setupAggregate! later or fields are nothing
```

**Proposed** (pipeline / smart constructors):
```julia
# Option A: pipeline
agg = AggregateCore(molecules, coupling) |> setup_aggregate

# Option B: single constructor (no partial state)
agg = Aggregate(molecules, coupling; vib_basis=GroundGround())
# Internally computes tools + operators. No nothing fields possible.
```

#### Pattern 5: Single canonical representation for time series

**Current** (dual representation, 4 overloads per function):
```julia
const OperatorVector = Vector{Operator{...}}
const OperatorVectorArray = Array{ComplexF64, 3}
# Every function needs: (OV,OV), (OV,OVA), (OVA,OV), (OVA,OVA)
```

**Proposed** (single type via SimulationResult from UX-2):
```julia
# All solvers return SimulationResult with Vector{Operator}
# One conversion function for interop:
to_array(r::SimulationResult) = stack(s.data for s in r.states)
```

#### Patterns to avoid in Julia

| OOP Pattern | Why it doesn't fit |
|-------------|-------------------|
| Composite | Aggregate isn't a tree -- flat struct is correct |
| Factory | Multiple dispatch constructors already do this |
| Singleton | Module-level constants suffice |
| Observer/Event | The existing `fout` callback is simpler |

#### Summary of pattern adoption

| Pattern | Replaces | Files affected | Blocked by |
|---------|----------|---------------|------------|
| Type dispatch | Symbol if/else chains | trace.jl, aggregateTools.jl, aggregateOperators.jl, master_ansatz.jl | #55 |
| Abstract ansatz types | `_bath_state_fn` lookup | master_ansatz.jl | #50 |
| Callable structs | Bare tuple `p`, closures | All solver files | #54 |
| Pipeline constructors | Union{Nothing} fields | aggregate.jl | Union{Nothing} issue |
| Single representation | OperatorVector/Array duality | postprocessing.jl, scoring.jl | UX-2 |

---

## UX Priority Order

| Priority | Issue | Impact |
|----------|-------|--------|
| 1 | **UX-1** Fix tutorial to use current API | Unblocks all new users |
| 2 | **UX-10** Adopt Julia-idiomatic patterns | Architectural foundation for everything below |
| 3 | **UX-2** Add `SimulationResult` type | Eliminates overload explosion, consistent API |
| 4 | **UX-3** Molecule-count coupling matrix | Removes biggest physicist confusion |
| 5 | **UX-5** Unified `solve()` entry point | Discoverable API, reduced cognitive load |
| 6 | **UX-4** Document naming conventions + glossary | Helps existing users read the code |
| 7 | **UX-6** Physical validation | Catches silent numerical errors |
| 8 | **UX-8** Solver selection guide in docs | Guides new users |
| 9 | **UX-7** Fix unit conversion | Small cleanup |
| 10 | **UX-9** Convenience constructors | Lowers barrier for common cases |
