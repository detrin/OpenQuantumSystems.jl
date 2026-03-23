
import QuadGK

"""
    exciton_basis(aggCore, aggTools)

Diagonalize the system Hamiltonian (excited-state block) to obtain exciton
energies and transformation coefficients.

Returns `(energies, coefficients)` where:
* `energies`: N-element vector of exciton energies (cm⁻¹).
* `coefficients`: N×N matrix where `coefficients[m, k]` is the contribution
  of site `k` to exciton state `m`.
"""
function exciton_basis(aggCore::AggregateCore, aggTools::AggregateTools;
                       vib_basis::VibBasisLike = GroundGround())
    Ham_sys = get_agg_ham_system_small(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = true)
    N = aggCore.molCount
    H_exc = Ham_sys.data[2:N+1, 2:N+1]
    vals, vecs = LinearAlgebra.eigen(LinearAlgebra.Hermitian(H_exc))
    return vals, vecs'
end

"""
    exciton_spectral_density(site_sds, coefficients, m, n)

Compute the spectral density in the exciton basis for exciton pair (m, n):

``\\tilde{S}_{mn,\\alpha,k} = c_{mk}^2 c_{nk}^2 S_{k\\alpha}``

Returns a combined [`SpectralDensity`](@ref) with all site contributions.
"""
function exciton_spectral_density(
    site_sds::Vector{SpectralDensity},
    coefficients::AbstractMatrix,
    m::Int,
    n::Int,
)::SpectralDensity
    all_freqs = Float64[]
    all_hrs = Float64[]
    N = length(site_sds)
    for k in 1:N
        c2_m = coefficients[m, k]^2
        c2_n = coefficients[n, k]^2
        weight = c2_m * c2_n
        abs(weight) < 1e-15 && continue
        for (ω, S) in zip(site_sds[k].frequencies, site_sds[k].huang_rhys)
            push!(all_freqs, ω)
            push!(all_hrs, weight * S)
        end
    end
    isempty(all_freqs) && return SpectralDensity([1.0], [0.0])
    return SpectralDensity(all_freqs, all_hrs)
end

"""
    modified_redfield_rates(aggCore; T=300.0, t_max=1000.0, gamma=0.0, rtol=1e-6, atol=1e-8)

Compute the modified Redfield population transfer rate matrix in the exciton basis.

Returns `(rates, energies, coefficients)` where:
* `rates`: N×N rate matrix. `rates[m, n]` is the rate from exciton `n` to exciton `m`.
  Diagonal elements are set so that columns sum to zero (population conservation).
* `energies`: exciton energies.
* `coefficients`: exciton transformation coefficients.

The rate from exciton state n to exciton state m is:

``R_{mn} = 2 \\operatorname{Re} \\int_0^\\infty dt\\, \\left[\\ddot{g}_{mn}(t) + \\left(\\dot{g}_{mm}(t) - \\dot{g}_{nn}(t)\\right) \\dot{g}_{mn}(t)\\right] e^{+i\\varepsilon_{mn}t - f_{mn}(t) - \\gamma t}``

where ``f_{mn}(t) = g_{mm}(t) + g_{nn}(t) - 2g_{mn}(t)`` and
``\\varepsilon_{mn} = \\varepsilon_m - \\varepsilon_n``.

# Arguments
* `aggCore`: [`AggregateCore`](@ref) instance.
* `T`: temperature in Kelvin (default 300).
* `t_max`: upper integration limit (default 1000, increase for slow dynamics).
* `gamma`: phenomenological damping rate for convergence (default 0, set > 0 for weakly damped systems).
* `rtol`, `atol`: tolerances for numerical integration.
"""
function modified_redfield_rates(
    aggCore::AggregateCore;
    T::Real = 300.0,
    t_max::Real = 1000.0,
    gamma::Real = 0.0,
    rtol::Real = 1e-6,
    atol::Real = 1e-8,
)
    aggTools = AggregateTools(aggCore)
    energies, coefficients = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount

    site_sds = [spectral_density(aggCore.molecules[k]) for k in 1:N]

    sd_cache = Dict{Tuple{Int,Int}, SpectralDensity}()
    function get_sd(m, n)
        get!(sd_cache, (m, n)) do
            exciton_spectral_density(site_sds, coefficients, m, n)
        end
    end

    rates = zeros(Float64, N, N)

    for m in 1:N
        for n in 1:N
            m == n && continue
            ε_mn = energies[m] - energies[n]
            sd_mn = get_sd(m, n)
            sd_mm = get_sd(m, m)
            sd_nn = get_sd(n, n)

            function integrand(t)
                g_mn = lineshape_function(sd_mn, t, T)
                g_mm = lineshape_function(sd_mm, t, T)
                g_nn = lineshape_function(sd_nn, t, T)
                gd_mn = lineshape_derivative(sd_mn, t, T)
                gd_mm = lineshape_derivative(sd_mm, t, T)
                gd_nn = lineshape_derivative(sd_nn, t, T)
                gdd_mn = lineshape_second_derivative(sd_mn, t, T)

                f_mn = g_mm + g_nn - 2.0 * g_mn
                bracket = gdd_mn + (gd_mm - gd_nn) * gd_mn
                phase = exp(im * ε_mn * t - f_mn - gamma * t)
                return bracket * phase
            end

            val, _ = QuadGK.quadgk(integrand, 0.0, Float64(t_max); rtol = rtol, atol = atol)
            rates[m, n] = 2.0 * real(val)
        end
    end

    for n in 1:N
        rates[n, n] = -sum(rates[m, n] for m in 1:N if m != n)
    end

    return rates, energies, coefficients
end

"""
    modified_redfield_dynamics(aggCore, rho0_exciton, tspan; T=300.0, t_max=1000.0, gamma=0.0)

Propagate exciton populations using modified Redfield rates.

# Arguments
* `aggCore`: [`AggregateCore`](@ref) instance.
* `rho0_exciton`: initial population vector in the exciton basis (N-element vector).
* `tspan`: time points for output.
* `T`: temperature in Kelvin.
* `t_max`: integration limit for rate calculation.
* `gamma`: phenomenological damping for rate convergence.

Returns `(tspan, populations)` where `populations` is a Matrix{Float64} of
size `(length(tspan), N)`.
"""
function modified_redfield_dynamics(
    aggCore::AggregateCore,
    rho0_exciton::AbstractVector{<:Real},
    tspan::AbstractVector;
    T::Real = 300.0,
    t_max::Real = 1000.0,
    gamma::Real = 0.0,
)
    rates, energies, coefficients = modified_redfield_rates(aggCore; T = T, t_max = t_max, gamma = gamma)
    N = length(energies)
    length(rho0_exciton) == N || throw(
        DimensionMismatch("rho0_exciton must have $N elements, got $(length(rho0_exciton))")
    )

    R = rates
    pop = Matrix{Float64}(undef, length(tspan), N)
    p = collect(Float64, rho0_exciton)
    pop[1, :] .= p

    for i in 2:length(tspan)
        dt = tspan[i] - tspan[i-1]
        dp = R * p
        p .= p .+ dt .* dp
        clamp!(p, 0.0, 1.0)
        s = sum(p)
        s > 0 && (p ./= s)
        pop[i, :] .= p
    end

    return collect(Float64, tspan), pop
end
