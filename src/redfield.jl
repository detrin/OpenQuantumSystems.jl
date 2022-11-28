using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq

function QME_sI_Redfield(
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
    rho0 = trace_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_Redfield(
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

function dQME_sI_Redfield(
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
    K_traced = trace_bath(K, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_Redfield(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_Redfield(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p

    rho = h(p, t)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, t-s)

    tmp1[:, :] = ad(rho, W0_bath.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
end