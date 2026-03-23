
function analytic_correlation_nn(sd::SpectralDensity, tau::Real, T::Real)::ComplexF64
    val = zero(ComplexF64)
    for (ω, S) in zip(sd.frequencies, sd.huang_rhys)
        n = _bose_einstein(ω, T)
        val += ω^2 * S * (n * exp(im * ω * tau) + (n + 1) * exp(-im * ω * tau))
    end
    return val
end

function analytic_correlation_nn_high_T(sd::SpectralDensity, tau::Real, T::Real)::Float64
    val = 0.0
    for (ω, S) in zip(sd.frequencies, sd.huang_rhys)
        n = _bose_einstein(ω, T)
        val += ω^2 * S * 2n * cos(ω * tau)
    end
    return val
end

function is_high_temperature(sd::SpectralDensity, T::Real; threshold::Real = 5.0)::Bool
    T <= 0.0 && return false
    kT = BOLTZMANN_CM * T
    for ω in sd.frequencies
        ω / kT > 1.0 / threshold && return false
    end
    return true
end

function exciton_correlation(
    C_nn_vals::AbstractVector{ComplexF64},
    coefficients::AbstractMatrix,
    a::Int, b::Int, c::Int, d::Int,
)::ComplexF64
    N = length(C_nn_vals)
    val = zero(ComplexF64)
    for n in 1:N
        val += coefficients[a, n] * coefficients[b, n] *
               coefficients[c, n] * coefficients[d, n] * C_nn_vals[n]
    end
    return val
end
