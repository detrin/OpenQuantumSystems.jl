
using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq


# ---------------------------------------------------------------------------
# Bath-state helper functions — one per ansatz variant
# Each returns (rho_for_ad, W_bath_data) given (s, p, rho)
# where rho is already a plain matrix (output of _to_matrix)
# ---------------------------------------------------------------------------

function _bath_state_test(s, p, rho)
    (; aggCore, aggTools, aggOperators, W0) = p
    Ham_0 = aggOperators.Ham_0
    Ham = aggOperators.Ham
    U_0_s = evolutionOperator(Ham, s)
    W_s = U_0_s * W0 * U_0_s'
    U_0_op = evolutionOperator(Ham_0, s)
    W_int_s = U_0_op' * W_s * U_0_op
    W_bath_s = get_rho_bath(W_int_s, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
    rho_s = trace_bath(W_int_s, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
    return rho_s.data, W_bath_s.data
end

function _bath_state_const_int(s, p, rho)
    (; W0_bath) = p
    return rho, W0_bath.data
end

function _bath_state_const_sch(s, p, rho)
    (; aggOperators, W0_bath) = p
    U_0_op = evolutionOperator(aggOperators.Ham_0, s)
    W0_int_s = U_0_op' * W0_bath * U_0_op
    return rho, W0_int_s.data
end

function _bath_state_linear_sch(s, p, rho)
    (; aggOperators, W0_bath, elementtype) = p
    Ham = aggOperators.Ham
    U_0_op = evolutionOperator(aggOperators.Ham_0, s)
    W_bath_s = W0_bath - elementtype(im) * (Ham * W0_bath - W0_bath * Ham) * s
    W0_int_s = U_0_op' * W_bath_s * U_0_op
    return rho, W0_int_s.data
end

function _bath_state_linear2_sch(s, p, rho)
    (; aggOperators, W0_bath, t_mk_bath_step, elementtype) = p
    Ham = aggOperators.Ham
    U_0_op = evolutionOperator(aggOperators.Ham_0, s)
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
    return rho, W0_int_s.data
end

function _bath_state_upart1_sch(s, p, rho)
    (; aggCore, aggTools, aggOperators, W0_bath) = p
    indicesMap = aggTools.indicesMap
    Ham = aggOperators.Ham
    Ham_0 = aggOperators.Ham_0
    U_0_op = evolutionOperator(Ham_0, s)
    W_bath_s = deepcopy(W0_bath)
    for el_st in 1:aggCore.molCount+1
        a1 = indicesMap[el_st][1]
        a2 = indicesMap[el_st][end]
        U_el_part = evolution_el_part(Ham.data, s, el_st, el_st, indicesMap)
        W_bath_s.data[a1:a2, a1:a2] = U_el_part * W_bath_s.data[a1:a2, a1:a2] * adjoint(U_el_part)
    end
    W0_int_s = U_0_op' * W_bath_s * U_0_op
    return rho, W0_int_s.data
end

function _bath_state_upart1_int(s, p, rho)
    (; aggCore, aggTools, aggOperators, W0_bath) = p
    indicesMap = aggTools.indicesMap
    Ham = aggOperators.Ham
    W_bath_s = deepcopy(W0_bath)
    for el_st in 1:aggCore.molCount+1
        a1 = indicesMap[el_st][1]
        a2 = indicesMap[el_st][end]
        U_el_part = evolution_el_part(Ham.data, s, el_st, el_st, indicesMap)
        W_bath_s.data[a1:a2, a1:a2] = U_el_part * W_bath_s.data[a1:a2, a1:a2] * adjoint(U_el_part)
    end
    return rho, W_bath_s.data
end

function _bath_state_upart2_sch(s, p, rho)
    (; aggCore, aggTools, aggOperators, W0_bath) = p
    indicesMap = aggTools.indicesMap
    Ham = aggOperators.Ham
    Ham_0 = aggOperators.Ham_0
    U_0_op = evolutionOperator(Ham_0, s)
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
    return rho, W0_int_s.data
end

function _bath_state_upart2_int(s, p, rho)
    (; aggCore, aggTools, aggOperators, W0_bath) = p
    indicesMap = aggTools.indicesMap
    Ham = aggOperators.Ham
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
    return rho, W_bath_s.data
end

# Dispatch table: Symbol -> bath-state function
function _bath_state_fn(ansatz::Symbol)
    ansatz == :test        && return _bath_state_test
    ansatz == :const_int   && return _bath_state_const_int
    ansatz == :const_sch   && return _bath_state_const_sch
    ansatz == :linear_sch  && return _bath_state_linear_sch
    ansatz == :linear2_sch && return _bath_state_linear2_sch
    ansatz == :upart1_sch  && return _bath_state_upart1_sch
    ansatz == :upart1_int  && return _bath_state_upart1_int
    ansatz == :upart2_sch  && return _bath_state_upart2_sch
    ansatz == :upart2_int  && return _bath_state_upart2_int
    throw(ArgumentError(
        "Unknown ansatz: :$ansatz. Valid options: :test, :const_int, :const_sch, " *
        ":linear_sch, :linear2_sch, :upart1_sch, :upart1_int, :upart2_sch, :upart2_int"
    ))
end


# ---------------------------------------------------------------------------
# Generic kernel
# ---------------------------------------------------------------------------

function kernel_sI_ansatz(t, s, h, p, tmp1, tmp2, Ham_II_t, ansatz::Symbol)
    (; aggCore, aggTools, aggOperators) = p
    rho = _to_matrix(h(p, s))
    Ham_II_s = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, s)
    bath_state_fn = _bath_state_fn(ansatz)
    rho_for_ad, W_bath_data = bath_state_fn(s, p, rho)
    tmp1[:, :] = ad(rho_for_ad, W_bath_data, aggCore, aggTools)
    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data
    return trace_bath(tmp1, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
end


# ---------------------------------------------------------------------------
# Generic dQME RHS
# ---------------------------------------------------------------------------

function dQME_sI_ansatz(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp1::AbstractMatrix,
    tmp2::AbstractMatrix,
    p,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
    ansatz::Symbol,
) where {B<:Basis,T<:Operator{B,B}}
    (; aggCore, aggTools, aggOperators, W0, elementtype) = p

    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)

    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_ansatz(t, s, history_fun, p, tmp1, tmp2, Ham_II_t, ansatz),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end


# ---------------------------------------------------------------------------
# Generic top-level solver
# ---------------------------------------------------------------------------

function QME_sI_ansatz(
    W0::T,
    tspan::AbstractVector,
    agg::Aggregate;
    ansatz::Symbol = :test,
    t_mk_bath_count::Int = 10,
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    _bath_state_fn(ansatz)  # validate ansatz early
    rho0 = trace_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    W0_bath = get_rho_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)
    t_mk_bath_step = (tspan[end] - tspan[1]) / t_mk_bath_count
    p = (
        aggCore=agg.core, aggTools=agg.tools, aggOperators=agg.operators,
        W0=W0, W0_bath=W0_bath, t_mk_bath_step=t_mk_bath_step, elementtype=eltype(W0),
    )

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_ansatz(
        t, rho, drho, history_fun, tmp1, tmp2, p, int_reltol, int_abstol, ansatz,
    )
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    tspan, rho_t = integrate_delayed(
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
    return tspan, rho_t
end


# ---------------------------------------------------------------------------
# Backward-compatible thin wrappers (one per original variant)
# ---------------------------------------------------------------------------

QME_sI_ansatz_test(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:test, kwargs...)

QME_sI_ansatz_const_int(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:const_int, kwargs...)

QME_sI_ansatz_const_sch(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:const_sch, kwargs...)

QME_sI_ansatz_linear_sch(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:linear_sch, kwargs...)

QME_sI_ansatz_linear2_sch(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:linear2_sch, kwargs...)

QME_sI_ansatz_upart1_sch(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:upart1_sch, kwargs...)

QME_sI_ansatz_upart1_int(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:upart1_int, kwargs...)

QME_sI_ansatz_upart2_sch(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:upart2_sch, kwargs...)

QME_sI_ansatz_upart2_int(W0, tspan, agg; kwargs...) =
    QME_sI_ansatz(W0, tspan, agg; ansatz=:upart2_int, kwargs...)
