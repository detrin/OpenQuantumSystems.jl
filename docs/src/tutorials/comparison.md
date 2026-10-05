```@meta
CurrentModule=OpenQuantumSystems
```

# Comparing Simulation Methods

This tutorial compares the different simulation methods available in
OpenQuantumSystems.jl on a standard molecular dimer benchmark. The system
parameters follow Section 3.4 of [Her23] and allow direct comparison of
population dynamics across approaches.

## Benchmark System

We use a homodimer with one vibrational mode per molecule:

| Parameter | Value |
|-----------|-------|
| Site energies | ``E_1 = 12\,500``, ``E_2 = 12\,700`` cm⁻¹ |
| Electronic coupling | ``J = 100`` cm⁻¹ |
| Vibrational frequency | ``\omega = 200`` cm⁻¹ |
| Temperature | ``T = 300`` K |
| Huang-Rhys factor | varied: ``S = 0.01, 0.05, 0.1`` |

The Huang-Rhys factor ``S`` controls the system-bath coupling strength.
Small ``S`` (weak coupling) favours Redfield-type approaches, while larger
``S`` tests their limits.

## Setup

```julia
using OpenQuantumSystems
using LinearAlgebra

function make_dimer(S; nvib=3, J=100.0, T=300.0)
    mode = Mode(; omega=200.0, hr_factor=S)
    mol1 = Molecule([mode], nvib, [0.0, 12500.0])
    mol2 = Molecule([mode], nvib, [0.0, 12700.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = J
    aggCore.coupling[3, 2] = J
    return aggCore
end
```

## Method 1: Exact Evolution (Reference)

The exact evolution propagates the full system+bath density matrix via
``W(t) = U(t) W(0) U^\dagger(t)`` and traces over the bath. This is the
reference against which all approximate methods are compared.

```julia
S = 0.05
nvib = 5
mode = Mode(; omega=200.0, hr_factor=S)
mol1 = Molecule([mode], nvib, [0.0, 12500.0])
mol2 = Molecule([mode], nvib, [0.0, 12700.0])
aggCore = AggregateCore([mol1, mol2])
aggCore.coupling[2, 3] = 100.0
aggCore.coupling[3, 2] = 100.0

agg = setup_aggregate(aggCore)
Ham = agg.operators.Ham

T = 300.0
W0 = thermal_state_composite(T, [0.0, 1.0, 0.0], agg)

t_max = 0.5
tspan = collect(0.0:0.001:t_max)
tspan_out, W_t = evolution_approximate(W0, tspan, Ham)

elLen = aggCore.molCount + 1
rho_exact = zeros(Float64, length(tspan), elLen)
for t_i in 1:length(tspan)
    rho_traced = trace_bath(W_t[t_i].data, agg)
    for n in 1:elLen
        rho_exact[t_i, n] = real(rho_traced[n, n])
    end
end
# rho_exact[:, 2] is site 1 population, rho_exact[:, 3] is site 2
```

The exact evolution is limited by the vibrational basis size: computation
scales as ``O(N_\text{vib}^{2N_\text{mol}})`` per time step.

## Method 2: QME with Redfield Approximation

The Redfield equation replaces the density matrix inside the memory integral
with its current-time value, yielding a time-local master equation. It works
well for weak coupling (small ``S``).

```julia
agg = setup_aggregate(aggCore)
W0 = thermal_state_composite(T, [0.0, 1.0, 0.0], agg)

# memory-kernel solvers: keep t_max short (cost grows faster than linearly)
tspan = collect(0.0:0.001:0.05)
tspan_out, rho_int_t = QME_sI_Redfield(W0, tspan, agg; reltol=1e-6, abstol=1e-8)

rho_redfield = zeros(Float64, length(tspan_out), elLen)
rho_sch_t = interaction_pic_to_schroedinger_pic(rho_int_t, tspan_out, agg)
for t_i in 1:length(tspan_out)
    for n in 1:elLen
        rho_redfield[t_i, n] = real(rho_sch_t[t_i, n, n])
    end
end
```

**When to use:** Weak coupling (``S \lesssim 0.05``), Markovian-like bath
(short correlation time relative to system dynamics).

## Method 3: QME with Bath Ansatz

The ansatz approach approximates the bath state inside the memory integral
with various levels of sophistication. The simplest ansatz (`:const_int`)
freezes the bath at its initial value in the interaction picture.

```julia
tspan_out, rho_int_t = QME_sI_ansatz(W0, tspan, agg; ansatz=:const_int,
    reltol=1e-6, abstol=1e-8)

rho_ansatz = zeros(Float64, length(tspan_out), elLen)
rho_sch_t = interaction_pic_to_schroedinger_pic(rho_int_t, tspan_out, agg)
for t_i in 1:length(tspan_out)
    for n in 1:elLen
        rho_ansatz[t_i, n] = real(rho_sch_t[t_i, n, n])
    end
end
```

**When to use:** General-purpose; the `:const_int` ansatz is the cheapest
nontrivial option. For better accuracy, try `:upart1_int` or `:upart2_int`.
The Schrödinger-picture variants (`:*_sch`) are considerably less accurate.

## Method 4: Modified Redfield Theory

Modified Redfield treats diagonal fluctuations non-perturbatively. It works
in the exciton basis and produces time-independent population transfer rates.

```julia
aggCore = make_dimer(0.05; nvib=3)

tspan = collect(0.0:0.001:0.5)
ts, pop_mr = modified_redfield_dynamics(aggCore, [1.0, 0.0], tspan; T=300.0)
# pop_mr[:, 1] = lower exciton population
# pop_mr[:, 2] = upper exciton population
```

!!! note
    Modified Redfield produces exciton-basis populations. To compare with
    exact (site-basis) results, transform using the exciton coefficients:
    ``\rho_k^{(\text{site})} = \sum_m |c_{mk}|^2 \rho_m^{(\text{exc})}``.

**When to use:** Intermediate coupling where standard Redfield breaks down
but full QME is too expensive. Provides rates, not full coherence dynamics.

## Method 5: Förster Theory

Förster theory computes EET rates from the spectral overlap of donor emission
and acceptor absorption. Valid in the weak electronic coupling limit.

```julia
aggCore = make_dimer(0.05; nvib=5)

K = forster_rate_matrix(aggCore; T=300.0, sigma=30.0)
# K[1,2] = rate from mol 2 to mol 1
# K[2,1] = rate from mol 1 to mol 2

# Propagate with Förster rates (simple forward Euler)
tspan = collect(0.0:0.001:0.5)
p = [1.0, 0.0]
pop_forster = zeros(Float64, length(tspan), 2)
pop_forster[1, :] .= p
for i in 2:length(tspan)
    dt = tspan[i] - tspan[i-1]
    dp = K * p
    p .= p .+ dt .* dp
    clamp!(p, 0.0, 1.0)
    s = sum(p); s > 0 && (p ./= s)
    pop_forster[i, :] .= p
end
```

**When to use:** Weak electronic coupling (``|J| \ll \lambda``), incoherent
energy transfer regime. Does not capture coherent oscillations.

## Method 6: Corrected Memory Kernel

The corrected memory kernel (thesis Section 3.7) improves on standard Redfield
by accounting for the bath's response to the system evolution, without
explicitly constructing the bath state. The zeroth-order kernel reproduces
standard Redfield-type rates; the first-order correction includes the
non-equilibrium bath relaxation.

```julia
aggCore = make_dimer(0.05; nvib=3)

# Zeroth-order corrected rates (Redfield-like)
r0 = corrected_rates_cf(aggCore; T=300.0, t_ref=0.1, order=0)

# First-order corrected rates (includes bath correction)
r1 = corrected_rates_cf(aggCore; T=300.0, t_ref=0.1, order=1,
    rtol=1e-3, atol=1e-5)

# Full iterative QME-RDM dynamics
tspan = collect(0.0:0.001:0.5)
t_out, pop_corr, all_iters = corrected_qme_rdm(aggCore, [1.0, 0.0], tspan;
    T=300.0, t_ref=0.1, max_iter=3, rtol=1e-3, atol=1e-5)

# all_iters[1] = zeroth-order populations
# all_iters[2] = first-order corrected populations
# all_iters[3] = second iteration (further refined)
```

!!! note
    The corrected kernel works in the exciton basis. The `t_ref` parameter
    sets the reference time for rate extraction. For short-time dynamics,
    use a small `t_ref`; for steady-state rates, use a larger value.

**When to use:** When you need systematic improvement over Redfield without
the cost of full QME. Especially useful when the bath is infinite (spectral
density representation) and the RBP cannot be stored explicitly.

## Results

All simulations start from the **upper exciton eigenstate** ``|\psi_2\rangle``,
which ensures all methods agree on the initial condition ``\rho_{22}^{\mathrm{exc}}(0) = 1``.
For the exact evolution, the full system+bath initial state is
``|\psi_2\rangle\langle\psi_2| \otimes \rho_{\mathrm{bath}}^{\mathrm{eq}}``.

### Weak coupling (``S = 0.01``)

![Comparison S=0.01](../assets/comparison_S0p01.png)

At weak coupling the upper exciton is nearly stationary: the **exact**
population oscillates between ``\approx 0.92`` and 1 (vibronic beating) and
slowly dephases.
- **Redfield** and the **QME ansatz** (`:const_int`) reproduce these
  oscillations, with a small phase drift and without the slow damping.
- **Modified Redfield** and the **corrected kernel** (both orders) predict
  no transfer at all.
- **Förster** decays spuriously to ``\approx 0.58``.

### Moderate coupling (``S = 0.05``)

![Comparison S=0.05](../assets/comparison_S0p05.png)

- The **exact** population drops to ``\approx 0.7`` within ``t \approx 0.03``
  and then oscillates, damped, around ``\approx 0.85``.
- **Redfield** and the **QME ansatz** follow the initial drop and the first
  minimum closely, but afterwards keep oscillating with undamped amplitude
  (Redfield between ``\approx 0.68`` and 1, the ansatz between
  ``\approx 0.58`` and 1).
- **Modified Redfield** relaxes to ``\approx 0.51`` and **Förster** to
  ``\approx 0.58``; neither shows the oscillations.
- The **corrected kernel** gives no transfer at zeroth order, and at first
  order a linear decay that crosses zero (unphysical).

### Intermediate coupling (``S = 0.1``)

![Comparison S=0.1](../assets/comparison_S0p1.png)

- The **exact** population drops to ``\approx 0.53`` and oscillates around
  ``\approx 0.75``.
- **Redfield** still captures the first drop and oscillates between
  ``\approx 0.5`` and 1; the **QME ansatz** overshoots strongly (down to
  ``\approx 0.24``).
- **Förster** relaxes to ``\approx 0.58``; **modified Redfield** and the
  zeroth-order corrected kernel show no transfer; the first-order corrected
  kernel drops below zero almost immediately.

### Discussion

1. **Rate-based methods only capture an envelope.** Modified Redfield,
   Förster and the corrected kernel solve ``d\mathbf{p}/dt = R\mathbf{p}``
   and cannot reproduce the coherent oscillations of the exact solution.
   Their predicted transfer is also qualitatively off here: either none, or
   relaxation to a fixed value regardless of ``S``.

2. **The second-order QMEs are accurate at short times.** Redfield and the
   `:const_int` ansatz follow the exact dynamics through the first
   oscillation at every coupling strength. Their error grows with ``S``, as
   expected for second-order theories, and the time-nonlocal ansatz degrades
   faster than Redfield at ``S = 0.1``.

3. **No dephasing from a single mode.** Each molecule couples to one undamped
   vibrational mode. With the bath frozen at equilibrium, the second-order
   memory kernel never decays, so the QMEs oscillate indefinitely, while the
   exact dynamics dephase through the vibronic structure of the finite bath.
   A continuous spectral density (an infinite bath) would damp the memory.

4. **The corrected kernel** is currently evaluated as a rate model at a fixed
   reference time `t_ref` (here 0.1). The zeroth-order rates vanish for this
   system and the first-order correction overshoots into negative
   populations. The time-nonlocal iterative QME with the corrected kernel
   (thesis, Ch. 3) is not yet implemented.

!!! note "Runtime"
    The memory-kernel solvers re-integrate the whole history at every step,
    so their cost grows roughly as ``t_\text{max}^2``. Generating these plots
    (`docs/generate_comparison_plots.jl`, `nvib = 5`, ``t_\text{max} = 0.5``)
    takes about 30–80 minutes per coupling strength.

## Full Comparison Script

This script runs all methods on a single dimer and collects the site-1
population for comparison:

```julia
using OpenQuantumSystems
using LinearAlgebra

S = 0.05
J = 100.0
T = 300.0
nvib = 5

mode = Mode(; omega=200.0, hr_factor=S)
mol1 = Molecule([mode], nvib, [0.0, 12500.0])
mol2 = Molecule([mode], nvib, [0.0, 12700.0])
aggCore = AggregateCore([mol1, mol2])
aggCore.coupling[2, 3] = J
aggCore.coupling[3, 2] = J

agg = setup_aggregate(aggCore)
Ham = agg.operators.Ham
elLen = aggCore.molCount + 1
t_max = 0.5
tspan = collect(0.0:0.001:t_max)

# --- Exact ---
W0 = thermal_state_composite(T, [0.0, 1.0, 0.0], agg)
_, W_t = evolution_approximate(W0, tspan, Ham)
rho11_exact = [real(trace_bath(W_t[i].data, agg)[2, 2]) for i in 1:length(tspan)]

# --- Redfield (memory-kernel solvers on a shorter window) ---
tspan_qme = collect(0.0:0.001:0.05)
W0 = thermal_state_composite(T, [0.0, 1.0, 0.0], agg)
t_out_r, rho_int_r = QME_sI_Redfield(W0, tspan_qme, agg; reltol=1e-6, abstol=1e-8)
rho_sch_r = interaction_pic_to_schroedinger_pic(rho_int_r, t_out_r, agg)
rho11_redfield = [real(rho_sch_r[i, 2, 2]) for i in 1:length(t_out_r)]

# --- Ansatz ---
W0 = thermal_state_composite(T, [0.0, 1.0, 0.0], agg)
t_out_a, rho_int_a = QME_sI_ansatz(W0, tspan_qme, agg; ansatz=:const_int,
    reltol=1e-6, abstol=1e-8)
rho_sch_a = interaction_pic_to_schroedinger_pic(rho_int_a, t_out_a, agg)
rho11_ansatz = [real(rho_sch_a[i, 2, 2]) for i in 1:length(t_out_a)]

# --- Modified Redfield (exciton → site) ---
_, coeffs = let
    at = AggregateTools(aggCore)
    exciton_basis(aggCore, at)
end
_, pop_mr = modified_redfield_dynamics(aggCore, [1.0, 0.0], tspan; T=T)
rho11_modred = [sum(coeffs[m, 1]^2 * pop_mr[i, m] for m in 1:2) for i in 1:length(tspan)]

# --- Corrected memory kernel (zeroth + first order) ---
t_corr, pop_corr, all_iters = corrected_qme_rdm(aggCore, [1.0, 0.0], tspan;
    T=T, t_ref=0.1, max_iter=2, rtol=1e-3, atol=1e-5)
rho11_corr0 = [sum(coeffs[m, 1]^2 * all_iters[1][i, m] for m in 1:2) for i in 1:length(tspan)]
rho11_corr1 = [sum(coeffs[m, 1]^2 * all_iters[2][i, m] for m in 1:2) for i in 1:length(tspan)]

# --- Förster ---
K = forster_rate_matrix(aggCore; T=T, sigma=30.0)
p = [1.0, 0.0]
rho11_forster = zeros(Float64, length(tspan))
rho11_forster[1] = 1.0
for i in 2:length(tspan)
    dt = tspan[i] - tspan[i-1]
    dp = K * p
    p .= p .+ dt .* dp
    clamp!(p, 0.0, 1.0)
    s = sum(p); s > 0 && (p ./= s)
    rho11_forster[i] = p[1]
end
```

To visualize:

```julia
using Plots
plot(tspan, rho11_exact, label="Exact", lw=2, color=:black)
plot!(t_out_r, rho11_redfield, label="Redfield", ls=:dash, color=:blue)
plot!(t_out_a, rho11_ansatz, label="Ansatz (const)", ls=:dash, color=:cyan)
plot!(tspan, rho11_modred, label="Modified Redfield", ls=:dot, color=:red)
plot!(tspan, rho11_corr0, label="Corrected (0th)", ls=:dashdot, color=:orange)
plot!(tspan, rho11_corr1, label="Corrected (1st)", ls=:dashdot, color=:green)
plot!(tspan, rho11_forster, label="Förster", ls=:dot, color=:purple)
xlabel!("t (a.u.)")
ylabel!("ρ₁₁(t)")
title!("Dimer population dynamics — S = $S")
ylims!(0.0, 1.0)
```

## Method Summary

| Method | Coupling regime | Basis | Coherences | Cost | Key function |
|--------|----------------|-------|------------|------|-------------|
| Exact evolution | Any | Site | Yes | Very high | `evolution_approximate` |
| QME Redfield | Weak | Site | Yes | Medium | `QME_sI_Redfield` |
| QME Ansatz | Moderate | Site | Yes | Medium | `QME_sI_ansatz` |
| Modified Redfield | Intermediate | Exciton | No (pops only) | Medium | `modified_redfield_dynamics` |
| Förster | Weak coupling | Site | No (rates only) | Low | `forster_rate_matrix` |
| Corrected kernel (0th) | Weak-moderate | Exciton | No (pops only) | Medium | `corrected_rates_cf` |
| Corrected kernel (1st) | Weak-moderate | Exciton | No (pops only) | High | `corrected_qme_rdm` |

## References

- [Her23] D. Herman, *Open Quantum Systems in Molecular Aggregates*, Master's thesis (2023).
  Sections 3.4 (Numerical results for QME approaches) and 3.6-3.7 (Iterative QME and corrected memory kernel).
