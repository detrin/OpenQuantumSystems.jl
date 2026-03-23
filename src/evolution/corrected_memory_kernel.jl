
function equilibrium_bath_state(
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    aggOperators::AggregateOperators,
    T::Real,
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    bSize = aggTools.bSize

    a1, a2 = indicesMap[1][1], indicesMap[1][end]
    H_0_g = real.(aggOperators.Ham_0.data[a1:a2, a1:a2])
    rho_B = exp(-H_0_g / (BOLTZMANN_CM * T))
    rho_B ./= LinearAlgebra.tr(rho_B)

    W_eq = zeros(ComplexF64, bSize, bSize)
    for a in 1:elLen, b in 1:elLen
        ai1, ai2 = indicesMap[a][1], indicesMap[a][end]
        bi1, bi2 = indicesMap[b][1], indicesMap[b][end]
        W_eq[ai1:ai2, bi1:bi2] .= rho_B
    end
    return W_eq
end

function zeroth_order_memory_kernel_cf(
    t1::Real, t2::Real,
    aggCore::AggregateCore;
    T::Real = 300.0,
)
    aggTools = AggregateTools(aggCore)
    energies, coefficients = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount
    elLen = N + 1

    site_sds = [spectral_density(aggCore.molecules[k]) for k in 1:N]

    tau = t1 - t2
    C_fwd = ComplexF64[analytic_correlation_nn(sd, tau, T) for sd in site_sds]
    C_bwd = ComplexF64[analytic_correlation_nn(sd, -tau, T) for sd in site_sds]

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools;
        vib_basis = GroundGround(), groundEnergy = true)
    E = zeros(Float64, elLen)
    E[1] = real(Ham_sys.data[1, 1])
    E[2:end] .= energies

    c = coefficients
    M = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    for a in 2:elLen, b in 2:elLen, ci in 2:elLen, d in 2:elLen
        ai, bi, cci, di = a - 1, b - 1, ci - 1, d - 1
        val = zero(ComplexF64)

        if d == b
            for e in 2:elLen
                ei = e - 1
                phase = exp(im * ((E[a] - E[e]) * t1 + (E[e] - E[ci]) * t2))
                val += phase * exciton_correlation(C_fwd, c, ai, ei, ei, cci)
            end
        end

        if a == ci
            for e in 2:elLen
                ei = e - 1
                phase = exp(im * ((E[d] - E[e]) * t2 + (E[e] - E[b]) * t1))
                val += phase * exciton_correlation(C_bwd, c, di, ei, ei, bi)
            end
        end

        phase3 = exp(im * ((E[a] - E[ci]) * t1 + (E[d] - E[b]) * t2))
        val -= phase3 * exciton_correlation(C_bwd, c, di, bi, ai, cci)

        phase4 = exp(im * ((E[a] - E[ci]) * t2 + (E[d] - E[b]) * t1))
        val -= phase4 * exciton_correlation(C_fwd, c, di, bi, ai, cci)

        M[a, b, ci, d] = val
    end

    return M
end

function site_to_exciton_kernel(
    M_site::Array{ComplexF64, 4},
    coefficients::AbstractMatrix,
)
    elLen = size(M_site, 1)
    N = elLen - 1

    T_mat = zeros(Float64, elLen, elLen)
    T_mat[1, 1] = 1.0
    T_mat[2:end, 2:end] .= coefficients

    M_exc = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    for α in 1:elLen, β in 1:elLen, γ in 1:elLen, δ in 1:elLen
        s = zero(ComplexF64)
        for a in 1:elLen, b in 1:elLen, cc in 1:elLen, d in 1:elLen
            s += T_mat[α, a] * T_mat[β, b] * T_mat[γ, cc] * T_mat[δ, d] * M_site[a, b, cc, d]
        end
        M_exc[α, β, γ, δ] = s
    end
    return M_exc
end
