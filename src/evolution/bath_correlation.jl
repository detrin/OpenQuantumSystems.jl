"""
    ExponentialCorrelation(a, z)

Bath correlation function of a single site in exponential-sum form,

``C(\\tau) = \\sum_k a_k e^{-z_k \\tau}, \\quad \\tau \\ge 0, \\qquad C(-\\tau) = C(\\tau)^*``,

with ``\\operatorname{Re} z_k \\ge 0``. Discrete vibrational modes are represented exactly
(``z = \\pm i\\omega``), a Drude-Lorentz bath up to the chosen number of Matsubara terms.
The convention is ``C(\\tau) = \\frac{1}{\\pi}\\int_0^\\infty d\\omega\\, J(\\omega)
[\\coth(\\beta\\omega/2)\\cos\\omega\\tau - i \\sin\\omega\\tau]``, i.e. the same as
[`analytic_correlation_nn`](@ref) for discrete modes.
"""
struct ExponentialCorrelation
    a::Vector{ComplexF64}
    z::Vector{ComplexF64}
    function ExponentialCorrelation(a::AbstractVector, z::AbstractVector)
        length(a) == length(z) || throw(DimensionMismatch("a and z must have the same length"))
        all(real.(z) .>= 0) || throw(ArgumentError("exponents must satisfy Re z ≥ 0"))
        new(collect(ComplexF64, a), collect(ComplexF64, z))
    end
end

function (C::ExponentialCorrelation)(τ::Real)
    τ < 0 && return conj(C(-τ))
    val = zero(ComplexF64)
    for k in eachindex(C.a)
        val += C.a[k] * exp(-C.z[k] * τ)
    end
    return val
end

"""
    bath_correlation(sd::SpectralDensity, T)

Exponential-sum correlation function of discrete modes at temperature `T` (K):
``C(\\tau) = \\sum_\\alpha \\omega_\\alpha^2 S_\\alpha [(\\bar n_\\alpha + 1) e^{-i\\omega_\\alpha\\tau}
+ \\bar n_\\alpha e^{i\\omega_\\alpha\\tau}]``.
"""
function bath_correlation(sd::SpectralDensity, T::Real)
    a = ComplexF64[]
    z = ComplexF64[]
    for (ω, S) in zip(sd.frequencies, sd.huang_rhys)
        n = _bose_einstein(ω, T)
        push!(a, ω^2 * S * (n + 1)); push!(z, im * ω)
        n > 0 && (push!(a, ω^2 * S * n); push!(z, -im * ω))
    end
    return ExponentialCorrelation(a, z)
end

"""
    DrudeLorentz(λ, γ)

Overdamped Brownian (Drude-Lorentz) spectral density
``J(\\omega) = 2\\lambda\\gamma\\omega / (\\omega^2 + \\gamma^2)`` with reorganization energy
`λ` and cut-off frequency `γ` (both in cm⁻¹). It describes an infinite bath whose
correlation function decays on the time scale ``1/\\gamma``.
"""
struct DrudeLorentz
    λ::Float64
    γ::Float64
end

spectral_density(dl::DrudeLorentz, ω::Real) = 2 * dl.λ * dl.γ * ω / (ω^2 + dl.γ^2)
reorganization_energy(dl::DrudeLorentz) = dl.λ

"""
    bath_correlation(dl::DrudeLorentz, T; n_matsubara=4)

Exponential-sum correlation function of a Drude-Lorentz bath at temperature `T` (K):

``C(\\tau) = \\lambda\\gamma[\\cot(\\beta\\gamma/2) - i] e^{-\\gamma\\tau}
+ \\frac{4\\lambda\\gamma}{\\beta}\\sum_{k=1}^{K}\\frac{\\nu_k}{\\nu_k^2-\\gamma^2} e^{-\\nu_k\\tau}``,

with Matsubara frequencies ``\\nu_k = 2\\pi k/\\beta``.
"""
function bath_correlation(dl::DrudeLorentz, T::Real; n_matsubara::Int = 4)
    T > 0 || throw(ArgumentError("Drude-Lorentz correlation requires T > 0"))
    λ, γ = dl.λ, dl.γ
    β = 1 / (BOLTZMANN_CM * T)
    a = ComplexF64[λ * γ * (cot(β * γ / 2) - im)]
    z = ComplexF64[γ]
    for k in 1:n_matsubara
        ν = 2π * k / β
        push!(a, 4 * λ * γ / β * ν / (ν^2 - γ^2))
        push!(z, ν)
    end
    return ExponentialCorrelation(a, z)
end

"""
    ExcitonModel(H_site, baths)
    ExcitonModel(aggCore; T)

Single-excitation manifold of an aggregate in the exciton basis: exciton energies
`energies`, coefficients `coefficients[α, n]` (site `n` in exciton `α`) and one
[`ExponentialCorrelation`](@ref) per site. `H_site` is the site Hamiltonian of the
excited manifold including the reorganization energies, i.e. the system part of
``H_0`` in the ground-state vibrational basis. The `aggCore` form uses the
discrete modes of the molecules.
"""
struct ExcitonModel
    energies::Vector{Float64}
    coefficients::Matrix{Float64}
    baths::Vector{ExponentialCorrelation}
end

function ExcitonModel(H_site::AbstractMatrix, baths::AbstractVector)
    size(H_site, 1) == length(baths) ||
        throw(DimensionMismatch("one bath per site is required"))
    vals, vecs = LinearAlgebra.eigen(LinearAlgebra.Hermitian(Matrix{Float64}(H_site)))
    return ExcitonModel(vals, Matrix(vecs'), collect(ExponentialCorrelation, baths))
end

function ExcitonModel(aggCore::AggregateCore; T::Real)
    energies, coefficients = exciton_basis(aggCore, AggregateTools(aggCore))
    baths = [bath_correlation(spectral_density(mol), T) for mol in aggCore.molecules]
    return ExcitonModel(energies, Matrix(coefficients), baths)
end

site_to_exciton(model::ExcitonModel, ρ_site::AbstractMatrix) =
    model.coefficients * ρ_site * model.coefficients'
exciton_to_site(model::ExcitonModel, ρ_exc::AbstractMatrix) =
    model.coefficients' * ρ_exc * model.coefficients
