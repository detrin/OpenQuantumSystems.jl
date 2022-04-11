
import OrdinaryDiffEq
import QuantumOpticsBase

"""
    liouvilleVonNeumann(rho0, tspan, H;
    \treltol=1.0e-12, abstol=1.0e-12, fout=nothing, alg=OrdinaryDiffEq.Tsit5())

Integrate Liouville-von Neumann equation to evolve states or compute propagators

``\\frac{d}{d t} \\rho(t) = - \\frac{i}{\\hbar} [ \\hat H, \\rho(t) ], \\quad \\hbar = 1``.

# Arguments
* `rho0`: Initial state vector (can be a bra or a ket) or initial propagator.
* `tspan`: Vector specifying the points of time for which output should be displayed.
* `H`: Arbitrary operator specifying the Hamiltonian.
* `reltol`: Relative tolerance for OrdinaryDiffEq solver and its inner states.
* `abstol`: Absolute tolerance for OrdinaryDiffEq solver and its inner states.
* `fout=nothing`: If given, this function `fout(t, rho)` is called every time
        an output should be displayed. ATTENTION: The state `rho` is neither
        normalized nor permanent! It is still in use by the ode solver and
        therefore must not be changed.
* `alg`: Algorithm with which OrdinaryDiffEq will solve LvN equation.
"""
function LvN_sI(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.Tsit5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    p = (agg.core, agg.tools, agg.operators, W0, eltype(W0))
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    dLvN_(t, rho, drho, p) = dLvN_sI(t, rho, drho, p)

    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate(
        tspan_,
        dLvN_,
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

function dLvN_sI(
    t::AbstractFloat,
    rho::T,
    drho::T,
    p
) where {B<:Basis,T<:Operator{B,B}}
    aggCore, aggTools, aggOperators, W0, elementtype = p
    Ham = aggOperators.Ham

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)

    U_t = evolutionOperator(Ham, t)
    W_t = U_t * W0 * U_t'
    U_0_op = evolutionOperator(Ham_0, t)
    W_int_t = U_0_op' * W_t * U_0_op

    K = -elementtype(im) * (Ham_II_t.data * W_int_t.data - W_int_t.data * Ham_II_t.data)
    drho.data[:, :] = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
    return drho
end

function LvN_sS(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.Tsit5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    p = (agg.core, agg.tools, agg.operators, W0, eltype(W0))
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    dLvN_(t, rho, drho, p) = dLvN_sS(t, rho, drho, p)

    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate(
        tspan_,
        dLvN_,
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

function dLvN_sS(
    t::AbstractFloat,
    rho::T,
    drho::T,
    p
) where {B<:Basis,T<:Operator{B,B}}
    aggCore, aggTools, aggOperators, W0, elementtype = p
    Ham = aggOperators.Ham

    U_t = evolutionOperator(Ham, t)
    W_t = U_t * W0 * U_t'

    K = -elementtype(im) * (Ham.data * W_t.data - W_t.data * Ham.data)
    drho.data[:, :] = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
    return drho
end


function LvN_SI(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.Tsit5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    p = (agg.operators, eltype(W0))
    dLvN_(t, W::T, dW::T, p) = dLvN_SI(t, W, dW, p)

    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = W0.data
    state = T(W0.basis_l, W0.basis_r, W0.data)
    dstate = T(W0.basis_l, W0.basis_r, W0.data)

    OpenQuantumSystems.integrate(
        tspan_,
        dLvN_,
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

function dLvN_SI(
    t::AbstractFloat,
    W::T,
    dW::T,
    p
) where {B<:Basis,T<:Operator{B,B}}
    aggOperators, elementtype = p

    Ham_0 = aggOperators.Ham_0.data
    Ham_I = aggOperators.Ham_I.data
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)

    dW.data[:, :] = -elementtype(im) * (Ham_II_t * W.data - W.data * Ham_II_t)
    return dW
end


function LvN_SS(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.Tsit5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    p = (agg.operators, eltype(W0))
    dLvN_(t, W::T, dW::T, p) = dLvN_SS(t, W, dW, p)

    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = W0.data
    state = T(W0.basis_l, W0.basis_r, W0.data)
    dstate = T(W0.basis_l, W0.basis_r, W0.data)

    OpenQuantumSystems.integrate(
        tspan_,
        dLvN_,
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

function dLvN_SS(
    t::AbstractFloat,
    W::T,
    dW::T,
    p
) where {B<:Basis,T<:Operator{B,B}}
    aggOperators, elementtype = p
    Ham = aggOperators.Ham.data
    dW.data[:, :] = -elementtype(im) * (Ham * W.data - W.data * Ham)
    return dW
end
