
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq

"""
    master_int(W0, tspan, Ham_0, Ham_I; 
    \treltol=1.0e-12, abstol=1.0e-12, int_reltol=1.0e-8, int_abstol=0.0, 
    \tfout=nothing, alg=DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()))

Integrate Quantum Master equation 

``\\frac{d}{d t} \\rho^{(I)}(t) = - \\frac{i}{\\hbar} [ \\hat{H}_I^{(I)}(t), \\rho^{(I)}(t_0) ] 
-\\frac{1}{\\hbar^2} \\int_{t_0}^{t_1} \\text{d} \\tau \\: [ \\hat{H}_I^{(I)}(t), [ \\hat{H}_I^{(I)}(\\tau), \\rho^{(I)}(\\tau) ]] ``

``H = H_S + H_B + H_I = H_0 + H_I, \\quad \\hbar = 1. ``

# Arguments
* `W0`: Initial state vector (can be a bra or a ket) or initial propagator.
* `tspan`: Vector specifying the points of time for which output should be displayed.s
* `Ham_0`: System and bath Hamiltonian as Operator.
* `Ham_I`: Interaction Hamiltonian as Operator.
* `reltol`: Relative tolerance for DiffEqCallbacks solver and its inner states.
* `abstol`: Absolute tolerance for DiffEqCallbacks solver and its inner states.
* `int_reltol`: Relative tolerance for QuadGK solver and its inner states.
* `int_abstol`: Absolute tolerance for QuadGK solver and its inner states.
* `fout=nothing`: If given, this function `fout(t, rho)` is called every time
        an output should be displayed. ATTENTION: The state `rho` is neither
        normalized nor permanent! It is still in use by the ode solver and
        therefore must not be changed.
* `alg`: Algorithm with which DiffEqCallbacks will solve QME equation.
"""
function master_int(
    rho0::T,
    tspan::Array,
    p;
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B},V<:Operator{B,B}}
    history_fun(p, t) = T(W0.basis_l, W0.basis_r, zeros(ComplexF64, size(W0.data)))
    Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, Ham_B, W0, W0_bath, agg, FCProd, indices, indicesMap, elementtype, aggCore, aggTools = p
    # (du,u,h,p,t)
    tmp = copy(W0.data)
    dmaster_(t, W::T, dW::T, history_fun, p) = dmaster_int(
        t,
        W,
        dW,
        history_fun,
        tmp,
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

function dmaster_int(
    t::AbstractFloat,
    W::T,
    dW::T,
    history_fun,
    tmp::Array,
    p,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B},V<:Operator{B,B}}
    Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, Ham_B, W0, W0_bath, agg, FCProd, indices, indicesMap, elementtype, aggCore, aggTools = p
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    elementtype = eltype(W.data)
    kernel_integrated, err = QuadGK.quadgk(
        s -> kernel_int(t, s, tmp, history_fun, p, Ham_II_t.data),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )
    # LinearAlgebra.mul!(tmp1, Ham_II_t.data, W0.data, -elementtype(im), zero(elementtype))
    # LinearAlgebra.mul!(tmp1, W0.data, Ham_II_t.data, elementtype(im), one(elementtype))
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggTools)
    kernel_integrated_traced = trace_bath(kernel_integrated, aggCore, aggTools)
    dW.data[:, :] = -elementtype(im) * K_traced[:, :] - kernel_integrated_traced[:, :] 
    # println(D(kernel_integrated , zero(kernel_integrated)))
    return dW
end

function kernel_int(t, s, tmp, h, p, Ham_II_t)
    Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, Ham_B, W0, W0_bath, agg, FCProd, indices, indicesMap, elementtype, aggCore, aggTools = p
    # Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, W0, W0_bath, agg, FCProd, indices, indicesMap, elementtype = p
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    rho_s = h(p, s)

    if (typeof(rho_s) <: Operator)
        rho_s = rho_s.data
    end
    # commutator(Ham.data, commutator(Ham.data, tmp))
    # QuantumOpticsBase.mul!(tmp, Ham_II_s, W_s, elementtype(1), zero(elementtype))
    # QuantumOpticsBase.mul!(tmp, W_s, Ham_II_s, -elementtype(1), one(elementtype))

    # QuantumOpticsBase.mul!(dW, Ham_II_t, tmp, elementtype(1), zero(elementtype))
    # QuantumOpticsBase.mul!(dW, tmp, Ham_II_t, -elementtype(1), one(elementtype))
    U_0_s = evolutionOperator(Ham_0+Ham_I, s)
    W_bath_s = U_0_s * W0 * W0_bath'
    U_0_op = evolutionOperator(Ham_0, s)
    W0_int_s = U_0_op' * W_bath_s * U_0_op
    # tmp2 .= ad(rho_s, W0_int_s.data, aggCore, aggTools)
    C1 = Ham_II_s.data * W0_int_s.data - W0_int_s.data * Ham_II_s.data
    C2 = Ham_II_t.data * C1 - C1 * Ham_II_t.data
    return C2
end

"""
    master(W0, tspan, Ham; 
    \treltol=1.0e-12, abstol=1.0e-12, int_reltol=1.0e-8, int_abstol=0.0, 
    \tfout=nothing, alg=DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()))

Integrate Quantum Master equation 

``\\frac{d}{d t} \\rho(t) = - \\frac{i}{\\hbar} [ \\hat{H}, \\rho(t_0) ] 
-\\frac{1}{\\hbar^2} \\int_{t_0}^{t_1} \\text{d} \\tau \\: [ \\hat{H}, [ \\hat{H}, \\rho(\\tau) ]] 
,\\quad \\hbar = 1. ``

# Arguments
* `W0`: Initial state vector (can be a bra or a ket) or initial propagator.
* `tspan`: Vector specifying the points of time for which output should be displayed.s
* `Ham`: Arbitrary operator specifying the Hamiltonian.
* `reltol`: Relative tolerance for DiffEqCallbacks solver and its inner states.
* `abstol`: Absolute tolerance for DiffEqCallbacks solver and its inner states.
* `int_reltol`: Relative tolerance for QuadGK solver and its inner states.
* `int_abstol`: Absolute tolerance for QuadGK solver and its inner states.
* `fout=nothing`: If given, this function `fout(t, rho)` is called every time
        an output should be displayed. ATTENTION: The state `rho` is neither
        normalized nor permanent! It is still in use by the ode solver and
        therefore must not be changed.
* `alg`: Algorithm with which DiffEqCallbacks will solve QME equation.
"""
function master(
    W0::T,
    tspan,
    agg::Aggregate;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-8,
    int_abstol::AbstractFloat = 0.0,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    p = (agg.operators, W0, eltype(W0))
    history_fun(p, t) = T(W0.basis_l, W0.basis_r, zeros(ComplexF64, size(W0.data)))
    # (du,u,h,p,t)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) =
        dmaster(t, rho, drho, history_fun, tmp1, tmp2, p, int_reltol, int_abstol)
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = W0.data
    state = T(W0.basis_l, W0.basis_r, W0.data)
    dstate = T(W0.basis_l, W0.basis_r, W0.data)
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

function dmaster(
    t::AbstractFloat,
    W::T,
    dW::T,
    history_fun,
    tmp1::Array,
    tmp2::Array,
    p,
    int_reltol::AbstractFloat = 1.0e-5,
    int_abstol::AbstractFloat = 1.0e-5,
) where {B<:Basis,T<:Operator{B,B}}
    aggOperators, W0, elementtype = p
    Ham = aggOperators.Ham.data

    kernel_integrated, err = QuadGK.quadgk(
        s -> kernel(t, s, tmp1, tmp2, history_fun, p),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )

    tmp1[:, :] = -elementtype(im) * (Ham * W0.data - W0.data * Ham)
    dW.data[:, :] = tmp1 - kernel_integrated
    return dW
end

function kernel(t, s, tmp1, tmp2, h, p)
    aggOperators, _, _ = p
    Ham = aggOperators.Ham.data
    W = h(p, s)
    if (typeof(W) <: Operator)
        W = W.data
    end

    tmp1[:, :] = Ham * W - W * Ham
    tmp2[:, :] = Ham * tmp1 - tmp1 * Ham

    return tmp2
end
