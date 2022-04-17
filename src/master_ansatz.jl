
using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq


function QME_sI_ansatz_test(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_test(
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

function dQME_sI_ansatz_test(
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
    aggCore, aggTools, aggOperators, W0, elementtype = p
        
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_test(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_test(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, _ = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)

    U_0_s = evolutionOperator(Ham, s)
    W_s = U_0_s * W0 * U_0_s'
    U_0_op = evolutionOperator(Ham_0, s)
    W_int_s = U_0_op' * W_s * U_0_op

    W_bath_s = get_rho_bath(W_int_s, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
    rho_s = trace_bath(W_int_s, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
    tmp1[:, :] = ad(rho_s.data, W_bath_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end

function QME_sI_ansatz_const_int(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_const_int(
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

function dQME_sI_ansatz_const_int(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_const_int(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_const_int(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)

    tmp1[:, :] = ad(rho, W0_bath.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end

function QME_sI_ansatz_const_sch(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_const_sch(
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

function dQME_sI_ansatz_const_sch(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_const_sch(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_const_sch(t, s, h, p, tmp1, tmp2, Ham_II_t)
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

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end


function QME_sI_ansatz_linear_sch(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_linear_sch(
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

function dQME_sI_ansatz_linear_sch(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_linear_sch(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_linear_sch(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, _, W0_bath, elementtype = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    
    U_0_op = evolutionOperator(Ham_0, s)
    # linear element
    W_bath_s = W0_bath - elementtype(im) * (Ham * W0_bath - W0_bath * Ham) * s
    W0_int_s = U_0_op' * W_bath_s * U_0_op
    tmp1[:, :] = ad(rho, W0_int_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end


function QME_sI_ansatz_linear2_sch(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    t_mk_bath_count::Integer = 100,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    t_mk_bath_step = (tspan[end] - tspan[1]) / t_mk_bath_count
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, t_mk_bath_step, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_linear2_sch(
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

function dQME_sI_ansatz_linear2_sch(
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
    aggCore, aggTools, aggOperators, W0, _, _, elementtype = p
        
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_linear2_sch(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_linear2_sch(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, _, W0_bath, t_mk_bath_step, elementtype = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    
    U_0_op = evolutionOperator(Ham_0, s)
    # linear element on intervals
    W_bath_s = deepcopy(W0_bath)
    s_rest = s
    while s_rest > t_mk_bath_step
        W_bath_s = W_bath_s - elementtype(im) * (Ham * W_bath_s - W_bath_s * Ham) * t_mk_bath_step
        s_rest -= t_mk_bath_step
    end
    if s_rest > 0
        W_bath_s = W_bath_s - elementtype(im) * (Ham * W_bath_s - W_bath_s * Ham) * s_rest
    end

    W0_int_s = U_0_op' * W_bath_s * U_0_op
    tmp1[:, :] = ad(rho, W0_int_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end

function QME_sI_ansatz_upart1_sch(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_upart1_sch(
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

function dQME_sI_ansatz_upart1_sch(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_upart1_sch(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_upart1_sch(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, _, W0_bath, elementtype = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    
    U_0_op = evolutionOperator(Ham_0, s)
    # linear element on intervals
    W_bath_s = deepcopy(W0_bath)
    for el_st in 1:aggCore.molCount+1
        a1 = indicesMap[el_st][1]
        a2 = indicesMap[el_st][end]
        U_el_part = evolution_el_part(Ham.data, s, el_st, el_st, indicesMap)
        W_bath_s.data[a1:a2, a1:a2] = U_el_part * W_bath_s.data[a1:a2, a1:a2] * adjoint(U_el_part)
    end

    W0_int_s = U_0_op' * W_bath_s * U_0_op
    tmp1[:, :] = ad(rho, W0_int_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end

function QME_sI_ansatz_upart1_int(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_upart1_int(
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

function dQME_sI_ansatz_upart1_int(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_upart1_int(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_upart1_int(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, _, W0_bath, elementtype = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    
    # linear element on intervals
    W_bath_s = deepcopy(W0_bath)
    for el_st in 1:aggCore.molCount+1
        a1 = indicesMap[el_st][1]
        a2 = indicesMap[el_st][end]
        U_el_part = evolution_el_part(Ham.data, s, el_st, el_st, indicesMap)
        W_bath_s.data[a1:a2, a1:a2] = U_el_part * W_bath_s.data[a1:a2, a1:a2] * adjoint(U_el_part)
    end

    tmp1[:, :] = ad(rho, W_bath_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end

function QME_sI_ansatz_upart2_sch(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_upart2_sch(
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

function dQME_sI_ansatz_upart2_sch(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_upart2_sch(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_upart2_sch(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, _, W0_bath, elementtype = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    
    U_0_op = evolutionOperator(Ham_0, s)
    # linear element on intervals
    W_bath_s = deepcopy(W0_bath)
    for el_st1 in 1:aggCore.molCount+1
        a1 = indicesMap[el_st1][1]
        a2 = indicesMap[el_st1][end]
        for el_st2 in 1:aggCore.molCount+1
            b1 = indicesMap[el_st2][1]
            b2 = indicesMap[el_st2][end]
            U_el_part = evolution_el_part(Ham.data, s, el_st1, el_st2, indicesMap)
            W_bath_s.data[a1:a2, b1:b2] = U_el_part * W_bath_s.data[a1:a2, b1:b2] * adjoint(U_el_part)
        end
    end

    W0_int_s = U_0_op' * W_bath_s * U_0_op
    tmp1[:, :] = ad(rho, W0_int_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end

function QME_sI_ansatz_upart2_int(
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
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz_upart2_int(
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

function dQME_sI_ansatz_upart2_int(
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
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz_upart2_int(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_ansatz_upart2_int(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, _, W0_bath, elementtype = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham = aggOperators.Ham
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    
    # linear element on intervals
    W_bath_s = deepcopy(W0_bath)
    for el_st1 in 1:aggCore.molCount+1
        a1 = indicesMap[el_st1][1]
        a2 = indicesMap[el_st1][end]
        for el_st2 in 1:aggCore.molCount+1
            b1 = indicesMap[el_st2][1]
            b2 = indicesMap[el_st2][end]
            U_el_part = evolution_el_part(Ham.data, s, el_st1, el_st2, indicesMap)
            W_bath_s.data[a1:a2, b1:b2] = U_el_part * W_bath_s.data[a1:a2, b1:b2] * adjoint(U_el_part)
        end
    end

    tmp1[:, :] = ad(rho, W_bath_s.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end