using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq

function QME_sI_Redfield(
    W0::T,
    tspan::AbstractVector,
    agg::Aggregate;
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    setup = _setup_delayed_integration(W0, tspan, agg)
    (; history_fun, tmp1, tmp2, tspan_, x0, state, dstate) = setup
    W0_bath = get_rho_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)
    p = (aggCore=agg.core, aggTools=agg.tools, aggOperators=agg.operators, W0=W0, W0_bath=W0_bath, elementtype=eltype(W0))

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
    integrate_delayed(
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
    tmp1::AbstractMatrix,
    tmp2::AbstractMatrix,
    p,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
) where {B<:Basis,T<:Operator{B,B}}
    (; aggCore, aggTools, aggOperators, W0, elementtype) = p
        
    Ham_II_t = get_interaction_ham_i_picture(aggOperators.Ham_0, aggOperators.Ham_I, t)
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
    (; aggCore, aggTools, aggOperators, W0, W0_bath) = p

    rho = h(p, t)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = get_interaction_ham_i_picture(Ham_0, Ham_I, t-s)

    tmp1[:, :] = ad(rho, W0_bath.data, aggCore, aggTools)

    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
end