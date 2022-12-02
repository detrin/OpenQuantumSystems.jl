using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq, Interpolations

function interpolate_with_tspan(itp, tspan, t)
    t0 = tspan[1]
    t1 = tspan[end]
    N = length(tspan)
    if t < t0
        t_i = 1
    elseif t > t1
        t_i = N
    else
        t_i = (t-t0)/t1 * (N-1) + 1
    end
    return itp[t_i]
end

function normalize_bath(W_bath, aggCore, aggTools, aggOperators)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    W_bath_tr = trace_bath(W_bath, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
    for b=1:elLen
        for a=1:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            W_bath[a1:a2, b1:b2] /= W_bath_tr[a, b]
        end
    end
    return W_bath
end


function W_aabb_1_bath_core_old(t, s, p, tmp1, tmp2; bath_evolution=:none, K_rtol=1e-12, K_atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, _ = p
    
    # W_1_bath = deepcopy(W0_bath.data)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham = aggOperators.Ham
    
    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    K_aabb_s = K_aabb_W_bath_intp(s, p, tmp1, tmp2; rtol=K_rtol, atol=K_atol)
    rho_s = interpolate_with_tspan(rho_0_int_t_itp, tspan, s)
    W_bath_s = interpolate_with_tspan(W_0_bath_t_itp, tspan, s)
    W_1_bath[:, :] = W_bath_s[:, :]
    for b=2:elLen
        if bath_evolution == :interaction_picture
            U_bb_ = evolution_el_part(Ham.data, s, b, b, indicesMap)
            U_0_bb = evolution_el_part(Ham_0.data, s, b, b, indicesMap)
            U_bb = U_0_bb * U_bb_
        elseif bath_evolution == :shroedinger_picture
            U_bb = evolution_el_part(Ham.data, s, b, b, indicesMap)
        end
        for a=2:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            if bath_evolution == :interaction_picture
                U_aa_ = evolution_el_part(Ham.data, t-s, a, a, indicesMap)
                U_0_aa = evolution_el_part(Ham_0.data, t-s, a, a, indicesMap)
                U_aa = U_0_aa * U_aa_
            elseif bath_evolution == :shroedinger_picture
                U_aa = evolution_el_part(Ham.data, t-s, a, a, indicesMap)
            end
            # 1 / rho_s[a, a] is already in K rate constant
            if bath_evolution == :none 
                W_1_bath[a1:a2, a1:a2] = W_1_bath[a1:a2, a1:a2] + 
                    K_aabb_s[a, b] * rho_s[b, b] * W_bath_s[b1:b2, b1:b2] 
            else
                U_ = U_aa * U_bb
                W_1_bath[a1:a2, a1:a2] = W_1_bath[a1:a2, a1:a2] + 
                    K_aabb_s[a, b] * rho_s[b, b]  * U_ * W_bath_s[b1:b2, b1:b2] * adjoint(U_)
            end
        end
    end
    return W_1_bath
end

function W_abcd_1_bath_core_old(t, s, p, tmp1, tmp2; bath_evolution=:none, K_rtol=1e-12, K_atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, _ = p
    
    # W_1_bath = deepcopy(W0_bath.data)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham = aggOperators.Ham
    
    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    K_abcd_s = K_abcd_W_bath_intp(s, p, tmp1, tmp2; rtol=K_rtol, atol=K_atol)
    rho_s = interpolate_with_tspan(rho_0_int_t_itp, tspan, s)
    W_bath_s = interpolate_with_tspan(W_0_bath_t_itp, tspan, s)
    W_1_bath[:, :] = W_bath_s[:, :]
    for c=2:elLen, d=2:elLen
        c1 = indicesMap[c][1]
        c2 = indicesMap[c][end]
        d1 = indicesMap[d][1]
        d2 = indicesMap[d][end]
        if bath_evolution == :interaction_picture
            U_cd_ = evolution_el_part(Ham.data, s, c, d, indicesMap)
            U_0_cd = evolution_el_part(Ham_0.data, s, c, d, indicesMap)
            U_cd = U_0_cd * U_cd_
        elseif bath_evolution == :shroedinger_picture
            U_cd = evolution_el_part(Ham.data, s, c, d, indicesMap)
        end
        for a=2:elLen, b=1:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            if bath_evolution == :interaction_picture
                U_ab_ = evolution_el_part(Ham.data, t-s, a, b, indicesMap)
                U_0_ab = evolution_el_part(Ham_0.data, t-s, a, b, indicesMap)
                U_ab = U_0_ab * U_ab_
            elseif bath_evolution == :shroedinger_picture
                U_aa = evolution_el_part(Ham.data, t-s, a, b, indicesMap)
            end
            # 1 / rho_s[a, b] is already in K rate constant
            if bath_evolution == :none 
                W_1_bath[a1:a2, b1:b2] = W_1_bath[a1:a2, b1:b2] + 
                    K_abcd_s[a, b, c, d] * rho_s[c, d] * W_bath_s[c1:c2, d1:d2] 
            else
                U_ = U_aa * U_bb
                W_1_bath[a1:a2, a1:a2] = W_1_bath[a1:a2, a1:a2] + 
                K_abcd_s[a, b, c, d] * rho_s[c, d] * U_ * W_bath_s[c1:c2, d1:d2]  * adjoint(U_)
            end
        end
    end
    return W_1_bath
end

function W_1_bath_old(
        t, p, tmp1, tmp2; 
        bath_evolution=:none, bath_ansatz=:population, normalize=false,
        W_1_rtol=1e-12, W_1_atol=1e-12, K_rtol=1e-12, K_atol=1e-12
    )
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    if bath_ansatz == :population
        W_1_diff, err = QuadGK.quadgk(
            s -> W_aabb_1_bath_core_old(t, s, p, tmp1, tmp2; bath_evolution=bath_evolution, K_rtol=K_rtol, K_atol=K_atol),
            0,
            t,
            rtol = W_1_rtol,
            atol = W_1_atol,
        )
    elseif bath_ansatz == :population_coherences
        W_1_diff, err = QuadGK.quadgk(
            s -> W_abcd_1_bath_core_old(t, s, p, tmp1, tmp2; bath_evolution=bath_evolution, K_rtol=K_rtol, K_atol=K_atol),
            0,
            t,
            rtol = W_1_rtol,
            atol = W_1_atol,
        )
    end
    
    W_1_bath = deepcopy(W0_bath.data) - W_1_diff
    if normalize
        W_1_bath = normalize_bath(W_1_bath, aggCore, aggTools, aggOperators)
    end
    return W_1_bath
end

function W_1_bath_old(t, W0, W0_bath, agg::Aggregate; W_1_rtol=1e-12, W_1_atol=1e-12, K_rtol=1e-12, K_atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))
    return W_1_bath(t, p, tmp1, tmp2; W_1_rtol=W_1_rtol, W_1_atol=W_1_atol, K_rtol=K_rtol, K_atol=K_atol)
end

function QME_sI_iterative_old(
    W0::T,
    rho_0_int_t,
    W_0_bath_t,
    tspan::Array,
    agg::Aggregate;
    bath_evolution=:none, 
    bath_ansatz=:population, 
    normalize=false,
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    W_1_rtol::AbstractFloat = 1e-12, 
    W_1_atol::AbstractFloat = 1e-12, 
    K_rtol::AbstractFloat = 1e-12, 
    K_atol::AbstractFloat = 1e-12,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    rho0 = trace_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    
    # Calculate and interpolate rho_0_int_t, W_0_bath_t, W_1_bath_t    
    if ndims(rho_0_int_t) == 1 && typeof(rho_0_int_t[1]) <: Operator
        rho_0_int_t = operator_recast(rho_0_int_t)
    end
    if ndims(rho_0_int_t) == 3
        rho_0_int_t = [rho_0_int_t[i, :, :] for i=1:size(rho_0_int_t, 1)]
    end
    rho_0_int_t_itp = Interpolations.interpolate(
        rho_0_int_t,
        Interpolations.BSpline(Interpolations.Linear())
    )
    
    if ndims(rho_0_int_t) == 1 && typeof(rho_0_int_t[1]) <: Operator
        W_0_bath_t = operator_recast(W_0_bath_t)
    end
    if ndims(W_0_bath_t) == 3
        W_0_bath_t = [W_0_bath_t[i, :, :] for i=1:size(W_0_bath_t, 1)]
    end
    W_0_bath_t_itp = Interpolations.interpolate(
        W_0_bath_t, 
        Interpolations.BSpline(Interpolations.Linear())
    )
    
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, eltype(W0))
    elLen = agg.core.molCount+1
    W_1_bath_t = []
    for t_i=1:length(tspan)
        t = tspan[t_i]
        W_1_bath_ = W_1_bath_old(
            t, p, tmp1, tmp2; 
            bath_evolution=bath_evolution, bath_ansatz=bath_ansatz, normalize=normalize,
            W_1_rtol=W_1_rtol, W_1_atol=W_1_atol, K_rtol=K_rtol, K_atol=K_atol
        )
        push!(W_1_bath_t, W_1_bath_)
    end
    W_1_bath_itp = Interpolations.interpolate(W_1_bath_t, Interpolations.BSpline(Interpolations.Linear()))
    
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, W_1_bath_itp, tspan, eltype(W0))
    
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_iterative_old(
        t,
        rho,
        drho,
        history_fun,
        tmp1,
        tmp2,
        p,
        int_reltol,
        int_abstol,
    )
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    tspan, rho_int_1_t = OpenQuantumSystems.integrate_delayed(
        tspan_,
        dmaster_,
        history_fun,
        x0,
        state,
        dstate,
        fout;
        p = p,
        reltol = reltol,
        abstol = abstol,
        alg = alg,
        kwargs...,
    )
    return tspan, rho_int_1_t, W_1_bath_t
end

function dQME_sI_iterative_old(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp1::Array,
    tmp2::Array,
    p,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
) where {B<:Basis,T<:Operator{B,B}}
    aggCore, aggTools, aggOperators, W0, _, _, _, elementtype = p
        
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_iterative_old(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_iterative_old(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, W_1_bath_itp, tspan, _ = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)

    W0_int_s = interpolate_with_tspan(W_1_bath_itp, tspan, s)
    tmp1[:, :] = ad(rho, W0_int_s, aggCore, aggTools)
    
    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
end