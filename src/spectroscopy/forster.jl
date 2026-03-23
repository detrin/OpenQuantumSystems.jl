
# Boltzmann constant in cm⁻¹/K, for use when molecular energies are in wavenumbers.
const BOLTZMANN_CM = 0.6950387

"""
    LineshapeResult

Container for a normalized molecular lineshape (absorption or emission spectrum).

# Fields
* `freqs`: frequency grid (same units as molecule energies, typically cm⁻¹).
* `intensities`: normalized spectral intensities (∫ intensities dν ≈ 1).
"""
struct LineshapeResult
    freqs::Vector{Float64}
    intensities::Vector{Float64}
end

function _build_sticks(
    mol::Molecule,
    init_el::Int,
    final_el::Int,
    init_weights::Vector{Float64},
    init_vib::Vector{Vector{Int}},
    final_vib::Vector{Vector{Int}},
    init_energies::Vector{Float64},
    final_energies::Vector{Float64},
)
    stick_freqs = Float64[]
    stick_ints = Float64[]
    for (n_idx, n_vib) in enumerate(init_vib)
        w = init_weights[n_idx]
        w < 1e-15 && continue
        for (m_idx, m_vib) in enumerate(final_vib)
            fc = get_mol_state_fc(mol, final_el, m_vib, init_el, n_vib)
            intensity = w * abs2(fc)
            intensity < 1e-15 && continue
            push!(stick_freqs, abs(final_energies[m_idx] - init_energies[n_idx]))
            push!(stick_ints, intensity)
        end
    end
    return stick_freqs, stick_ints
end

function _boltzmann_weights(energies::Vector{Float64}, T::Real)::Vector{Float64}
    T <= 0.0 && (weights = zeros(length(energies)); weights[1] = 1.0; return weights)
    kT = BOLTZMANN_CM * T
    w = exp.(-(energies .- energies[1]) ./ kT)
    return w ./ sum(w)
end

function _broaden_and_normalize(
    stick_freqs::Vector{Float64},
    stick_ints::Vector{Float64},
    freqs::Vector{Float64},
    sigma::Float64,
)::Vector{Float64}
    isempty(stick_freqs) && return zeros(length(freqs))
    ints = zeros(Float64, length(freqs))
    inv2s2 = 1.0 / (2 * sigma^2)
    norm_factor = 1.0 / (sigma * sqrt(2π))
    for (ω0, I0) in zip(stick_freqs, stick_ints)
        @. ints += I0 * norm_factor * exp(-((freqs - ω0)^2) * inv2s2)
    end
    dω = freqs[2] - freqs[1]
    total = sum(ints) * dω
    total > 0 && (ints ./= total)
    return ints
end

function _auto_freq_grid(
    stick_freqs::Vector{Float64},
    stick_ints::Vector{Float64},
    sigma::Float64,
    n_points::Int,
)::Vector{Float64}
    ω_avg = sum(stick_freqs .* stick_ints) / sum(stick_ints)
    half_range = max(maximum(abs.(stick_freqs .- ω_avg)), sigma) + 5 * sigma
    return collect(LinRange(ω_avg - half_range, ω_avg + half_range, n_points))
end

"""
    absorption_spectrum(mol; T=0.0, sigma=50.0, n_points=2000)

Compute the normalized absorption spectrum of a [`Molecule`](@ref).

Transitions from thermally-weighted ground-electronic vibronic states to all
excited-electronic vibronic states, broadened by a Gaussian of width `sigma`.

# Arguments
* `mol`: [`Molecule`](@ref) instance.
* `T`: temperature in Kelvin (energies assumed in cm⁻¹; uses `BOLTZMANN_CM`).
* `sigma`: Gaussian broadening standard deviation (same units as energies; default 50 cm⁻¹).
* `n_points`: frequency grid resolution (default 2000).
"""
function absorption_spectrum(
    mol::Molecule;
    T::Real = 0.0,
    sigma::Real = 50.0,
    n_points::Int = 2000,
)::LineshapeResult
    vib_inds = vibrational_indices(fill(mol.Nvib, length(mol.modes)))
    gnd_E = [get_mol_state_energy(mol, ELECTRONIC_GROUND, v) for v in vib_inds]
    exc_E = [get_mol_state_energy(mol, ELECTRONIC_EXCITED, v) for v in vib_inds]
    gnd_w = _boltzmann_weights(gnd_E, T)

    stick_freqs, stick_ints = _build_sticks(
        mol, ELECTRONIC_GROUND, ELECTRONIC_EXCITED, gnd_w, vib_inds, vib_inds, gnd_E, exc_E
    )

    freqs = _auto_freq_grid(stick_freqs, stick_ints, Float64(sigma), n_points)
    ints = _broaden_and_normalize(stick_freqs, stick_ints, freqs, Float64(sigma))
    return LineshapeResult(freqs, ints)
end

"""
    emission_spectrum(mol; T=0.0, sigma=50.0, n_points=2000)

Compute the normalized emission spectrum of a [`Molecule`](@ref).

Transitions from thermally-weighted excited-electronic vibronic states (Kasha's rule)
to all ground-electronic vibronic states, broadened by a Gaussian of width `sigma`.

# Arguments
* `mol`: [`Molecule`](@ref) instance.
* `T`: temperature in Kelvin.
* `sigma`: Gaussian broadening standard deviation (default 50 cm⁻¹).
* `n_points`: frequency grid resolution (default 2000).
"""
function emission_spectrum(
    mol::Molecule;
    T::Real = 0.0,
    sigma::Real = 50.0,
    n_points::Int = 2000,
)::LineshapeResult
    vib_inds = vibrational_indices(fill(mol.Nvib, length(mol.modes)))
    gnd_E = [get_mol_state_energy(mol, ELECTRONIC_GROUND, v) for v in vib_inds]
    exc_E = [get_mol_state_energy(mol, ELECTRONIC_EXCITED, v) for v in vib_inds]
    exc_w = _boltzmann_weights(exc_E, T)

    stick_freqs, stick_ints = _build_sticks(
        mol, ELECTRONIC_EXCITED, ELECTRONIC_GROUND, exc_w, vib_inds, vib_inds, exc_E, gnd_E
    )

    freqs = _auto_freq_grid(stick_freqs, stick_ints, Float64(sigma), n_points)
    ints = _broaden_and_normalize(stick_freqs, stick_ints, freqs, Float64(sigma))
    return LineshapeResult(freqs, ints)
end

function _interp_spectrum(ls::LineshapeResult, freqs::Vector{Float64})::Vector{Float64}
    src = ls.freqs
    isempty(src) && return zeros(length(freqs))
    dω = src[2] - src[1]
    result = zeros(Float64, length(freqs))
    for (i, ω) in enumerate(freqs)
        (ω < src[1] || ω > src[end]) && continue
        t = (ω - src[1]) / dω
        lo = max(1, floor(Int, t) + 1)
        hi = min(length(src), lo + 1)
        frac = t - (lo - 1)
        result[i] = (1 - frac) * ls.intensities[lo] + frac * ls.intensities[hi]
    end
    return result
end

"""
    spectral_overlap(donor_emission, acceptor_absorption)

Compute the spectral overlap integral:

``J = \\int I_D(\\omega) A_A(\\omega)\\, d\\omega``

Both spectra must be normalized (∫ dν = 1). Returns the overlap in units of
1/[frequency], where frequency is in the same units as the molecule energies.

Uses linear interpolation onto the donor emission frequency grid.
"""
function spectral_overlap(donor_em::LineshapeResult, acceptor_abs::LineshapeResult)::Float64
    freqs = donor_em.freqs
    I_D = donor_em.intensities
    A_A = _interp_spectrum(acceptor_abs, freqs)
    dω = freqs[2] - freqs[1]
    return sum(I_D .* A_A) * dω
end

"""
    forster_rate(J, mol_donor, mol_acceptor; T=0.0, sigma=50.0, n_points=2000)

Compute the Förster excitation energy transfer rate:

``k_{DA} = 2\\pi |J_{DA}|^2 \\int I_D(\\omega) A_A(\\omega)\\, d\\omega``

# Arguments
* `J`: electronic coupling (cm⁻¹).
* `mol_donor`: donor [`Molecule`](@ref).
* `mol_acceptor`: acceptor [`Molecule`](@ref).
* `T`: temperature in Kelvin (default 0).
* `sigma`: Gaussian broadening width in cm⁻¹ (default 50).
* `n_points`: frequency grid resolution (default 2000).

# Returns
Rate in cm⁻¹. To convert to 1/fs: multiply by `2π × c` where `c = 3×10⁻⁵` cm/fs.
"""
function forster_rate(
    J::Real,
    mol_donor::Molecule,
    mol_acceptor::Molecule;
    T::Real = 0.0,
    sigma::Real = 50.0,
    n_points::Int = 2000,
)::Float64
    I_D = emission_spectrum(mol_donor; T = T, sigma = sigma, n_points = n_points)
    A_A = absorption_spectrum(mol_acceptor; T = T, sigma = sigma, n_points = n_points)
    S = spectral_overlap(I_D, A_A)
    return 2π * J^2 * S
end

"""
    forster_rate_matrix(aggCore; T=0.0, sigma=50.0, n_points=2000)

Compute the matrix of pairwise Förster excitation energy transfer rates for
all molecules in an [`AggregateCore`](@ref).

Returns an N×N matrix where entry (i,j) is the Förster rate from molecule i
to molecule j. Diagonal entries are zero.

The electronic coupling ``J_{ij}`` is read from `aggCore.coupling[i+1, j+1]`
(the +1 offset accounts for the ground state row/column in the padded coupling matrix).

# Arguments
* `aggCore`: [`AggregateCore`](@ref) instance.
* `T`: temperature in Kelvin (default 0).
* `sigma`: Gaussian broadening width in cm⁻¹ (default 50).
* `n_points`: frequency grid resolution (default 2000).
"""
function forster_rate_matrix(
    aggCore::AggregateCore;
    T::Real = 0.0,
    sigma::Real = 50.0,
    n_points::Int = 2000,
)::Matrix{Float64}
    N = aggCore.molCount
    em_spectra = [emission_spectrum(aggCore.molecules[i]; T = T, sigma = sigma, n_points = n_points) for i in 1:N]
    abs_spectra = [absorption_spectrum(aggCore.molecules[i]; T = T, sigma = sigma, n_points = n_points) for i in 1:N]

    K = zeros(Float64, N, N)
    for i in 1:N
        for j in 1:N
            i == j && continue
            J_ij = aggCore.coupling[i+1, j+1]
            S = spectral_overlap(em_spectra[i], abs_spectra[j])
            K[i, j] = 2π * J_ij^2 * S
        end
    end
    return K
end
