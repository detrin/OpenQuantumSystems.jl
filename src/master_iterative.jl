using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq


function W_1_bath_core(t, s, p, tmp1, tmp2; K_rtol=1e-12, K_atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    # W_1_bath = deepcopy(W0_bath.data)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham = aggOperators.Ham
    
    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    K_ab_s = K_ab(s, p, tmp1, tmp2; rtol=K_rtol, atol=K_atol)
    for b=1:elLen
        U_bb_ = evolution_el_part(Ham.data, s, b, b, indicesMap)
        U_0_bb = evolution_el_part(Ham_0.data, s, b, b, indicesMap)
        U_bb = U_0_bb * U_bb_
        for a=1:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            U_aa_ = evolution_el_part(Ham.data, t-s, a, a, indicesMap)
            U_0_aa = evolution_el_part(Ham_0.data, t-s, a, a, indicesMap)
            U_aa = U_0_aa * U_aa_
            U_ = U_aa * U_bb
            W_1_bath[a1:a2, a1:a2] = W_1_bath[a1:a2, a1:a2] + 
                K_ab_s[a, b] * U_ * W0_bath.data[b1:b2, b1:b2] * adjoint(U_)
        end
    end
    return W_1_bath
end

function W_1_bath(t, p, tmp1, tmp2; W_1_rtol=1e-12, W_1_atol=1e-12, K_rtol=1e-12, K_atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    W_1_diff, err = QuadGK.quadgk(
        s -> W_1_bath_core(t, s, p, tmp1, tmp2; K_rtol=K_rtol, K_atol=K_atol),
        0,
        t,
        rtol = W_1_rtol,
        atol = W_1_atol,
    )
    
    W_1_bath = deepcopy(W0_bath.data) + W_1_diff
    return W_1_bath
end

function W_1_bath(t, W0, W0_bath, agg::Aggregate; W_1_rtol=1e-12, W_1_atol=1e-12, K_rtol=1e-12, K_atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))
    return W_1_bath(t, p, tmp1, tmp2; W_1_rtol=W_1_rtol, W_1_atol=W_1_atol, K_rtol=K_rtol, K_atol=K_atol)
end

function normalize_bath(W_bath, aggCore, aggTools)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    W_bath_tr = trace_bath(W_bath, aggCore, aggTools)
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

function QME_sI_ansatz_const(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    rho0 = trace_bath(W0, agg.core, agg.tools)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_const(
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
    OpenQuantumSystems.integrate_delayed(
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
end

function dQME_sI_ansatz_const(
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
    aggCore, aggTools, aggOperators, W0, _, elementtype = p
        
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggTools)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_const(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_const(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)

    
    U_0_op = evolutionOperator(Ham_0, s)
    W0_int_s = U_0_op' * W0_bath * U_0_op
    tmp1[:, :] = ad(rho, W0_int_s.data, aggCore, aggTools)
    

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools)
end