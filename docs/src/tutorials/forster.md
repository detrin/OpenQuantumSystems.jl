```@meta
CurrentModule=OpenQuantumSystems
```

# Förster Theory and Dipole Couplings

This tutorial shows how to compute absorption and emission spectra, Förster
excitation energy transfer rates, and how to set up aggregates with dipole-dipole
couplings.

## Absorption and Emission Spectra

Each molecule in an aggregate has vibronic structure. The absorption spectrum
shows transitions from thermally-weighted ground vibronic states to excited
vibronic states; the emission spectrum shows the reverse.

```julia
using OpenQuantumSystems

# Define a molecule with one vibrational mode
mode = Mode(; omega=200.0, hr_factor=0.05)
mol = Molecule([mode], 5, [0.0, 12500.0])

# Compute spectra at T=300K
abs_spec = absorption_spectrum(mol; T=300.0, sigma=20.0)
em_spec  = emission_spectrum(mol; T=300.0, sigma=20.0)

# abs_spec.freqs and abs_spec.intensities are vectors
# The spectra are normalized: ∫ intensities dν ≈ 1
```

The absorption spectrum peaks near the 0-0 transition energy with a vibronic
progression to higher frequencies. The emission spectrum is the mirror image,
red-shifted by the Stokes shift.

## Spectral Overlap and Förster Rate

The Förster rate between a donor and acceptor depends on the overlap of the
donor's emission with the acceptor's absorption:

```julia
# Two molecules with different excitation energies
mol_donor    = Molecule([mode], 5, [0.0, 12500.0])
mol_acceptor = Molecule([mode], 5, [0.0, 12700.0])

# Electronic coupling
J = 100.0  # cm⁻¹

# Single-pair Förster rate
k = forster_rate(J, mol_donor, mol_acceptor; T=300.0, sigma=30.0)
```

## Förster Rate Matrix for an Aggregate

For an aggregate with N molecules, `forster_rate_matrix` computes all pairwise
rates from the coupling matrix:

```julia
mol1 = Molecule([mode], 5, [0.0, 12500.0])
mol2 = Molecule([mode], 5, [0.0, 12700.0])
mol3 = Molecule([mode], 5, [0.0, 12600.0])

aggCore = AggregateCore([mol1, mol2, mol3])
aggCore.coupling[2, 3] = 100.0  # J₁₂
aggCore.coupling[3, 2] = 100.0
aggCore.coupling[2, 4] = 50.0   # J₁₃
aggCore.coupling[4, 2] = 50.0
aggCore.coupling[3, 4] = 80.0   # J₂₃
aggCore.coupling[4, 3] = 80.0

K = forster_rate_matrix(aggCore; T=300.0, sigma=30.0)
# K[i,j] is the rate from molecule j to molecule i
```

## Dipole-Dipole Couplings

Instead of manually setting couplings, you can compute them from molecular
positions and transition dipole vectors:

```julia
using LinearAlgebra

dipoles = [
    TransitionDipole([0.0,  0.0, 0.0], [0.0, 1.0, 0.0]),  # mol 1
    TransitionDipole([10.0, 0.0, 0.0], [0.0, 1.0, 0.0]),  # mol 2
    TransitionDipole([20.0, 0.0, 0.0], [0.0, 1.0, 0.0]),  # mol 3
]

# Compute the coupling matrix
J = coupling_from_dipoles(dipoles)

# Or construct the AggregateCore directly
aggCore = AggregateCore([mol1, mol2, mol3], dipoles)
```

The coupling follows the point-dipole formula:
```math
J_{DA} = \frac{\text{DIPOLE\_COUPLING\_FACTOR}}{R^3}
\left(\boldsymbol{\mu}_D \cdot \boldsymbol{\mu}_A
- 3(\boldsymbol{\mu}_D \cdot \hat{\mathbf{R}})(\boldsymbol{\mu}_A \cdot \hat{\mathbf{R}})\right)
```

where positions are in Å, dipole magnitudes in Debye, and the result is in cm⁻¹.

## Modified Redfield Theory

For intermediate coupling regimes, modified Redfield theory provides population
transfer rates in the exciton basis:

```julia
mode = Mode(; omega=200.0, hr_factor=0.3)
mol1 = Molecule([mode], 3, [0.0, 12500.0])
mol2 = Molecule([mode], 3, [0.0, 12700.0])
aggCore = AggregateCore([mol1, mol2])
aggCore.coupling[2, 3] = 100.0
aggCore.coupling[3, 2] = 100.0

# Compute rates and exciton basis
rates, energies, coefficients = modified_redfield_rates(aggCore; T=300.0)

# Propagate exciton populations starting from exciton 2
tspan = collect(0.0:0.01:5.0)
ts, pop = modified_redfield_dynamics(aggCore, [0.0, 1.0], tspan; T=300.0)
# pop[:,1] = population of lower exciton over time
# pop[:,2] = population of upper exciton over time
```

## Spectral Density and Lineshape Functions

The spectral density and lineshape functions used by modified Redfield theory
are also available directly:

```julia
# Extract spectral density from a molecule
sd = spectral_density(mol1)

# Reorganization energy
λ = reorganization_energy(sd)

# Lineshape function g(t) and its derivatives
g  = lineshape_function(sd, 0.01, 300.0)
gd = lineshape_derivative(sd, 0.01, 300.0)
```
