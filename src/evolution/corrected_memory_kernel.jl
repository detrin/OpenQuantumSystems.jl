import QuadGK

function _wick_second_order_cf(
    c::AbstractMatrix, site_sds::Vector{SpectralDensity}, T::Real, N::Int,
    i1::Int, i2::Int, i3::Int, i4::Int, i5::Int, i6::Int, i7::Int, i8::Int,
    ta::Real, tb::Real, tc::Real, td::Real,
)::ComplexF64
    val = zero(ComplexF64)
    for n in 1:N, m in 1:N
        K_nnmm = c[i1,n]*c[i2,n]*c[i3,n]*c[i4,n] * c[i5,m]*c[i6,m]*c[i7,m]*c[i8,m]
        K_nmnm = c[i1,n]*c[i2,n]*c[i5,n]*c[i6,n] * c[i3,m]*c[i4,m]*c[i7,m]*c[i8,m]
        K_nmmn = c[i1,n]*c[i2,n]*c[i7,n]*c[i8,n] * c[i3,m]*c[i4,m]*c[i5,m]*c[i6,m]
        if abs(K_nnmm) + abs(K_nmnm) + abs(K_nmmn) < 1e-30
            continue
        end
        Cn_ab = analytic_correlation_nn(site_sds[n], ta - tb, T)
        Cm_cd = analytic_correlation_nn(site_sds[m], tc - td, T)
        Cn_ac = analytic_correlation_nn(site_sds[n], ta - tc, T)
        Cm_bd = analytic_correlation_nn(site_sds[m], tb - td, T)
        Cn_ad = analytic_correlation_nn(site_sds[n], ta - td, T)
        Cm_bc = analytic_correlation_nn(site_sds[m], tb - tc, T)
        val += K_nnmm * Cn_ab * Cm_cd + K_nmnm * Cn_ac * Cm_bd + K_nmmn * Cn_ad * Cm_bc
    end
    return val
end

# [H_x, [H_y, O]] = H_x H_y O - H_x O H_y - H_y O H_x + O H_y H_x,
# each term as (sign, time slots left of O, time slots right of O)
_double_commutator_terms(x::Int, y::Int) = (
    (+1, (x, y), ()),
    (-1, (x,), (y,)),
    (-1, (y,), (x,)),
    (+1, (), (y, x)),
)

# Electronic paths α → β through `nops` single-exciton operators.
function _exciton_paths(α::Int, β::Int, nops::Int, elLen::Int)
    nops == 0 && return α == β ? [(α,)] : NTuple{1,Int}[]
    nops == 1 && return [(α, β)]
    return [(α, g, β) for g in 2:elLen]
end

# Integrand of the first-order correction of the memory kernel,
#   tr_B <a| [H(t1), [H(t2), (1/ρ_cd) <c|[H(t3), [H(t4), ρ ⊕ w_eq]]|d> |c><d|]] |b>,
# for all exciton indices (a, b, c, d) of the excited manifold. The sixteen operator
# orderings are generated from both double commutators; the four bath operators enter
# `wk` in cyclic trace order tr{L_o L_i w R_i R_o} = <R_i R_o L_o L_i>.
function _first_order_integrand!(
    result::Array{ComplexF64, 4},
    t1::Real, t2::Real, t3::Real, t4::Real,
    E::Vector{Float64},
    rho_I::AbstractMatrix,
    elLen::Int,
    wk::Function,
)
    times = (t1, t2, t3, t4)
    phase(tslot, α, β) = exp(im * (E[α] - E[β]) * times[tslot])
    outer = _double_commutator_terms(1, 2)
    inner = _double_commutator_terms(3, 4)
    ops = Vector{NTuple{3,Int}}(undef, 4)

    for a in 2:elLen, b in 2:elLen, c in 2:elLen, d in 2:elLen
        rho_cd = rho_I[c, d]
        if abs(rho_cd) < _SAFE_DIV_TOL
            result[a, b, c, d] = zero(ComplexF64)
            continue
        end
        val = zero(ComplexF64)
        for e in 2:elLen, f in 2:elLen
            rho_ef = rho_I[e, f]
            abs(rho_ef) < _SAFE_DIV_TOL && continue
            val_ef = zero(ComplexF64)
            for (s_o, L_o, R_o) in outer, (s_i, L_i, R_i) in inner
                for p_Lo in _exciton_paths(a, c, length(L_o), elLen),
                    p_Li in _exciton_paths(c, e, length(L_i), elLen),
                    p_Ri in _exciton_paths(f, d, length(R_i), elLen),
                    p_Ro in _exciton_paths(d, b, length(R_o), elLen)

                    k = 0
                    for (slots, path) in ((R_i, p_Ri), (R_o, p_Ro), (L_o, p_Lo), (L_i, p_Li))
                        for (j, tslot) in enumerate(slots)
                            k += 1
                            ops[k] = (tslot, path[j], path[j+1])
                        end
                    end
                    ph = prod(phase(op...) for op in ops)
                    (o1, o2, o3, o4) = ops
                    val_ef += s_o * s_i * ph * wk(
                        o1[2]-1, o1[3]-1, o2[2]-1, o2[3]-1, o3[2]-1, o3[3]-1, o4[2]-1, o4[3]-1,
                        times[o1[1]], times[o2[1]], times[o3[1]], times[o4[1]],
                    )
                end
            end
            val += (rho_ef / rho_cd) * val_ef
        end
        result[a, b, c, d] = val
    end
end

function _wick_second_order_cf_high_T(
    c::AbstractMatrix, site_sds::Vector{SpectralDensity}, T::Real, N::Int,
    i1::Int, i2::Int, i3::Int, i4::Int, i5::Int, i6::Int, i7::Int, i8::Int,
    ta::Real, tb::Real, tc::Real, td::Real,
)::Float64
    val = 0.0
    for n in 1:N, m in 1:N
        K_nnmm = c[i1,n]*c[i2,n]*c[i3,n]*c[i4,n] * c[i5,m]*c[i6,m]*c[i7,m]*c[i8,m]
        K_nmnm = c[i1,n]*c[i2,n]*c[i5,n]*c[i6,n] * c[i3,m]*c[i4,m]*c[i7,m]*c[i8,m]
        K_nmmn = c[i1,n]*c[i2,n]*c[i7,n]*c[i8,n] * c[i3,m]*c[i4,m]*c[i5,m]*c[i6,m]
        if abs(K_nnmm) + abs(K_nmnm) + abs(K_nmmn) < 1e-30
            continue
        end
        Cn_ab = analytic_correlation_nn_high_T(site_sds[n], ta - tb, T)
        Cm_cd = analytic_correlation_nn_high_T(site_sds[m], tc - td, T)
        Cn_ac = analytic_correlation_nn_high_T(site_sds[n], ta - tc, T)
        Cm_bd = analytic_correlation_nn_high_T(site_sds[m], tb - td, T)
        Cn_ad = analytic_correlation_nn_high_T(site_sds[n], ta - td, T)
        Cm_bc = analytic_correlation_nn_high_T(site_sds[m], tb - tc, T)
        val += K_nnmm * Cn_ab * Cm_cd + K_nmnm * Cn_ac * Cm_bd + K_nmmn * Cn_ad * Cm_bc
    end
    return val
end

function first_order_memory_kernel_high_T(
    t1::Real, t2::Real,
    aggCore::AggregateCore,
    rho_I::AbstractMatrix;
    T::Real = 300.0,
    rtol::Real = 1e-6,
    atol::Real = 1e-8,
)
    M0 = zeroth_order_memory_kernel_high_T(t1, t2, aggCore; T = T)
    if t2 <= 0.0
        return M0
    end

    aggTools = AggregateTools(aggCore)
    energies, coefficients = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount
    elLen = N + 1

    site_sds = [spectral_density(aggCore.molecules[k]) for k in 1:N]

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools;
        vib_basis = GroundGround(), groundEnergy = true)
    E = zeros(Float64, elLen)
    E[1] = real(Ham_sys.data[1, 1])
    E[2:end] .= energies

    cc = coefficients

    wk(i1,i2,i3,i4,i5,i6,i7,i8, ta,tb,tc,td) =
        _wick_second_order_cf_high_T(cc, site_sds, T, N, i1,i2,i3,i4,i5,i6,i7,i8, ta,tb,tc,td)

    correction, _ = QuadGK.quadgk(0.0, t2; rtol=rtol, atol=atol) do t3
        inner, _ = QuadGK.quadgk(0.0, t3; rtol=rtol, atol=atol) do t4
            result = zeros(ComplexF64, elLen, elLen, elLen, elLen)
            _first_order_integrand!(result, t1, t2, t3, t4, E, rho_I, elLen, wk)
            return result
        end
        return inner
    end

    return M0 .- correction
end

function first_order_memory_kernel_cf(
    t1::Real, t2::Real,
    aggCore::AggregateCore,
    rho_I::AbstractMatrix;
    T::Real = 300.0,
    rtol::Real = 1e-6,
    atol::Real = 1e-8,
)
    M0 = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = T)
    if t2 <= 0.0
        return M0
    end

    aggTools = AggregateTools(aggCore)
    energies, coefficients = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount
    elLen = N + 1

    site_sds = [spectral_density(aggCore.molecules[k]) for k in 1:N]

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools;
        vib_basis = GroundGround(), groundEnergy = true)
    E = zeros(Float64, elLen)
    E[1] = real(Ham_sys.data[1, 1])
    E[2:end] .= energies

    cc = coefficients

    wk(i1,i2,i3,i4,i5,i6,i7,i8, ta,tb,tc,td) =
        _wick_second_order_cf(cc, site_sds, T, N, i1,i2,i3,i4,i5,i6,i7,i8, ta,tb,tc,td)

    correction, _ = QuadGK.quadgk(0.0, t2; rtol=rtol, atol=atol) do t3
        inner, _ = QuadGK.quadgk(0.0, t3; rtol=rtol, atol=atol) do t4
            result = zeros(ComplexF64, elLen, elLen, elLen, elLen)
            _first_order_integrand!(result, t1, t2, t3, t4, E, rho_I, elLen, wk)
            return result
        end
        return inner
    end

    return M0 .- correction
end

function _thermal_rho_exciton(aggCore::AggregateCore, T::Real)
    aggTools = AggregateTools(aggCore)
    energies, _ = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount
    elLen = N + 1
    rho = zeros(ComplexF64, elLen, elLen)
    rho[1, 1] = 1.0
    boltz = [exp(-E / (BOLTZMANN_CM * T)) for E in energies]
    Z = sum(boltz)
    for i in 1:N
        rho[i+1, i+1] = boltz[i] / Z
    end
    return rho
end

function corrected_rates_cf(
    aggCore::AggregateCore;
    T::Real = 300.0,
    t_ref::Real = 1.0,
    rho_I::Union{AbstractMatrix, Nothing} = nothing,
    order::Int = 1,
    rtol::Real = 1e-4,
    atol::Real = 1e-6,
)
    N = aggCore.molCount
    _rho = rho_I === nothing ? _thermal_rho_exciton(aggCore, T) : rho_I

    val, _ = QuadGK.quadgk(0.0, t_ref; rtol = rtol, atol = atol) do t2
        M = if order == 0
            zeroth_order_memory_kernel_cf(t_ref, t2, aggCore; T = T)
        else
            first_order_memory_kernel_cf(t_ref, t2, aggCore, _rho;
                T = T, rtol = rtol, atol = atol)
        end
        result = zeros(Float64, N, N)
        for m in 1:N, n in 1:N
            result[m, n] = real(M[m+1, m+1, n+1, n+1])
        end
        return result
    end

    rates = val
    for n in 1:N
        rates[n, n] = -sum(rates[m, n] for m in 1:N if m != n; init = 0.0)
    end
    return rates
end

function _propagate_populations(
    rates::AbstractMatrix, rho0::AbstractVector, tspan::AbstractVector,
)
    N = length(rho0)
    pop = Matrix{Float64}(undef, length(tspan), N)
    p = collect(Float64, rho0)
    pop[1, :] .= p
    for i in 2:length(tspan)
        dt = tspan[i] - tspan[i-1]
        dp = rates * p
        p .= p .+ dt .* dp
        clamp!(p, 0.0, 1.0)
        s = sum(p)
        s > 0 && (p ./= s)
        pop[i, :] .= p
    end
    return pop
end

function corrected_qme_rdm(
    aggCore::AggregateCore,
    rho0_exciton::AbstractVector{<:Real},
    tspan::AbstractVector;
    T::Real = 300.0,
    t_ref::Real = 1.0,
    max_iter::Int = 2,
    rtol::Real = 1e-4,
    atol::Real = 1e-6,
)
    N = aggCore.molCount
    elLen = N + 1
    length(rho0_exciton) == N || throw(
        DimensionMismatch("rho0_exciton must have $N elements, got $(length(rho0_exciton))")
    )

    rho_I = _thermal_rho_exciton(aggCore, T)
    all_pops = Matrix{Float64}[]

    for iter in 0:max_iter-1
        order = iter == 0 ? 0 : 1
        rates = corrected_rates_cf(aggCore;
            T = T, t_ref = t_ref, rho_I = rho_I, order = order,
            rtol = rtol, atol = atol)

        pop = _propagate_populations(rates, rho0_exciton, tspan)
        push!(all_pops, pop)

        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        for i in 1:N
            rho_I[i+1, i+1] = pop[end, i]
        end
    end

    return collect(Float64, tspan), all_pops[end], all_pops
end

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

function zeroth_order_memory_kernel_high_T(
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
    C_real = Float64[analytic_correlation_nn_high_T(sd, tau, T) for sd in site_sds]

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools;
        vib_basis = GroundGround(), groundEnergy = true)
    E = zeros(Float64, elLen)
    E[1] = real(Ham_sys.data[1, 1])
    E[2:end] .= energies

    c = coefficients
    M = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    C_fwd_cx = ComplexF64.(C_real)
    C_bwd_cx = ComplexF64.(C_real)

    for a in 2:elLen, b in 2:elLen, ci in 2:elLen, d in 2:elLen
        ai, bi, cci, di = a - 1, b - 1, ci - 1, d - 1
        val = zero(ComplexF64)

        if d == b
            for e in 2:elLen
                ei = e - 1
                phase = exp(im * ((E[a] - E[e]) * t1 + (E[e] - E[ci]) * t2))
                val += phase * exciton_correlation(C_fwd_cx, c, ai, ei, ei, cci)
            end
        end

        if a == ci
            for e in 2:elLen
                ei = e - 1
                phase = exp(im * ((E[d] - E[e]) * t2 + (E[e] - E[b]) * t1))
                val += phase * exciton_correlation(C_bwd_cx, c, di, ei, ei, bi)
            end
        end

        phase3 = exp(im * ((E[a] - E[ci]) * t1 + (E[d] - E[b]) * t2))
        val -= phase3 * exciton_correlation(C_bwd_cx, c, di, bi, ai, cci)

        phase4 = exp(im * ((E[a] - E[ci]) * t2 + (E[d] - E[b]) * t1))
        val -= phase4 * exciton_correlation(C_fwd_cx, c, di, bi, ai, cci)

        M[a, b, ci, d] = val
    end

    return M
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
