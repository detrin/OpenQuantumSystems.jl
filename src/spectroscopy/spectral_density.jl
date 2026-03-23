
"""
    SpectralDensity(frequencies, huang_rhys)

Discrete spectral density defined by a set of vibrational mode frequencies
and their Huang-Rhys factors. Used to compute lineshape functions for
Förster and modified Redfield theories.

# Fields
* `frequencies`: mode frequencies (cm⁻¹).
* `huang_rhys`: Huang-Rhys factors ``S_\\alpha`` for each mode.
"""
struct SpectralDensity
    frequencies::Vector{Float64}
    huang_rhys::Vector{Float64}
    function SpectralDensity(frequencies::AbstractVector{<:Real}, huang_rhys::AbstractVector{<:Real})
        length(frequencies) == length(huang_rhys) ||
            throw(DimensionMismatch("frequencies and huang_rhys must have the same length"))
        new(collect(Float64, frequencies), collect(Float64, huang_rhys))
    end
end

"""
    spectral_density(mol::Molecule)

Extract a [`SpectralDensity`](@ref) from a [`Molecule`](@ref) by reading the
frequency and Huang-Rhys factor (shift/2) from each mode.
"""
function spectral_density(mol::Molecule)::SpectralDensity
    freqs = [mode.omega for mode in mol.modes]
    hrs = [mode.shift / 2.0 for mode in mol.modes]
    return SpectralDensity(freqs, hrs)
end

"""
    reorganization_energy(sd::SpectralDensity)

Compute the reorganization energy ``\\lambda = \\sum_\\alpha S_\\alpha \\omega_\\alpha``.
"""
function reorganization_energy(sd::SpectralDensity)::Float64
    return sum(sd.huang_rhys .* sd.frequencies)
end

function _bose_einstein(omega::Float64, T::Real)::Float64
    T <= 0.0 && return 0.0
    kT = BOLTZMANN_CM * T
    x = omega / kT
    x > 500.0 && return 0.0
    return 1.0 / (exp(x) - 1.0)
end

"""
    lineshape_function(sd, t, T)

Compute the lineshape function ``g(t)`` for a discrete spectral density at
temperature `T` (Kelvin):

``g(t) = \\sum_\\alpha S_\\alpha \\left[(2\\bar{n}_\\alpha + 1)(1 - \\cos\\omega_\\alpha t) + i(\\sin\\omega_\\alpha t - \\omega_\\alpha t)\\right]``

Returns a complex number.
"""
function lineshape_function(sd::SpectralDensity, t::Real, T::Real)::ComplexF64
    g = zero(ComplexF64)
    for (ω, S) in zip(sd.frequencies, sd.huang_rhys)
        n = _bose_einstein(ω, T)
        coth_part = 2n + 1.0
        ωt = ω * t
        g += S * (coth_part * (1.0 - cos(ωt)) + im * (sin(ωt) - ωt))
    end
    return g
end

"""
    lineshape_derivative(sd, t, T)

Compute the first derivative ``\\dot{g}(t)`` of the lineshape function.

``\\dot{g}(t) = \\sum_\\alpha S_\\alpha \\omega_\\alpha \\left[(2\\bar{n}_\\alpha + 1)\\sin\\omega_\\alpha t + i(\\cos\\omega_\\alpha t - 1)\\right]``
"""
function lineshape_derivative(sd::SpectralDensity, t::Real, T::Real)::ComplexF64
    gd = zero(ComplexF64)
    for (ω, S) in zip(sd.frequencies, sd.huang_rhys)
        n = _bose_einstein(ω, T)
        coth_part = 2n + 1.0
        ωt = ω * t
        gd += S * ω * (coth_part * sin(ωt) + im * (cos(ωt) - 1.0))
    end
    return gd
end

"""
    lineshape_second_derivative(sd, t, T)

Compute the second derivative ``\\ddot{g}(t)`` of the lineshape function.

``\\ddot{g}(t) = \\sum_\\alpha S_\\alpha \\omega_\\alpha^2 \\left[(2\\bar{n}_\\alpha + 1)\\cos\\omega_\\alpha t - i\\sin\\omega_\\alpha t\\right]``
"""
function lineshape_second_derivative(sd::SpectralDensity, t::Real, T::Real)::ComplexF64
    gdd = zero(ComplexF64)
    for (ω, S) in zip(sd.frequencies, sd.huang_rhys)
        n = _bose_einstein(ω, T)
        coth_part = 2n + 1.0
        ωt = ω * t
        gdd += S * ω^2 * (coth_part * cos(ωt) - im * sin(ωt))
    end
    return gdd
end
