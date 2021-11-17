
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
    W0::T,
    tspan::Array,
    Ham_0::U,
    Ham_I::V;
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-8,
    int_abstol::AbstractFloat = 0.0,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B},V<:Operator{B,B}}
    history_fun(p, t) = T(W0.basis_l, W0.basis_r, zeros(ComplexF64, size(W0.data)))
    # (du,u,h,p,t)
    tmp = copy(W0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) = dmaster_int(
        t,
        rho,
        drho,
        history_fun,
        tmp,
        p,
        W0,
        Ham_0,
        Ham_I,
        int_reltol,
        int_abstol,
    )
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
        reltol = reltol,
        abstol = abstol,
        alg = alg,
        kwargs...,
    )
end

function dmaster_int(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp::Array,
    p,
    W0,
    Ham_0::U,
    Ham_I::V,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B},V<:Operator{B,B}}
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    QuantumOpticsBase.mul!(drho, Ham_II_t, W0, -eltype(rho)(im), zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho, W0, Ham_II_t, eltype(rho)(im), one(eltype(rho)))

    kernel_integrated, err = QuadGK.quadgk(
        s -> kernel_int(t, s, tmp, history_fun, p, Ham_II_t.data, Ham_0, Ham_I),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )
    LinearAlgebra.mul!(
        drho.data,
        -one(eltype(rho)),
        kernel_integrated,
        one(eltype(rho)),
        one(eltype(rho)),
    )
    return drho
end

function kernel_int(t, s, tmp, h, p, Ham_II_t, Ham_0, Ham_I)
    Ham_II_s = getInteractionHamIPictureA(Ham_0, Ham_I, s)
    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end
    drho = deepcopy(rho)
    # commutator(Ham.data, commutator(Ham.data, tmp))
    QuantumOpticsBase.mul!(tmp, Ham_II_s, rho, eltype(rho)(1), zero(eltype(rho)))
    QuantumOpticsBase.mul!(tmp, rho, Ham_II_s, -eltype(rho)(1), one(eltype(rho)))

    QuantumOpticsBase.mul!(drho, Ham_II_t, tmp, eltype(rho)(1), zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho, tmp, Ham_II_t, -eltype(rho)(1), one(eltype(rho)))
    return drho
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
    Ham::U;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-8,
    int_abstol::AbstractFloat = 0.0,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B}}
    history_fun(p, t) = T(W0.basis_l, W0.basis_r, zeros(ComplexF64, size(W0.data)))
    # (du,u,h,p,t)
    tmp = copy(W0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) =
        dmaster(t, rho, drho, history_fun, tmp, p, W0, Ham)
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
        reltol = reltol,
        abstol = abstol,
        alg = alg,
        kwargs...,
    )
end

function dmaster(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp::Array,
    p,
    W0,
    Ham::U,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B}}
    QuantumOpticsBase.mul!(drho, Ham, W0, -eltype(rho)(im), zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho, W0, Ham, eltype(rho)(im), one(eltype(rho)))

    kernel_integrated, err =
        QuadGK.quadgk(s -> kernel(t, s, tmp, history_fun, p, Ham.data), 0, t, rtol = 1e-8)
    LinearAlgebra.mul!(
        drho.data,
        -one(eltype(rho)),
        kernel_integrated,
        one(eltype(rho)),
        one(eltype(rho)),
    )
    return drho
end

function kernel(t, s, tmp, h, p, Ham)
    rho = h(p, s)
    if (typeof(rho) <: Operator)
        rho = rho.data
    end
    drho = deepcopy(rho)
    # commutator(Ham.data, commutator(Ham.data, tmp))
    QuantumOpticsBase.mul!(tmp, Ham, rho, eltype(rho)(1), zero(eltype(rho)))
    QuantumOpticsBase.mul!(tmp, rho, Ham, -eltype(rho)(1), one(eltype(rho)))

    QuantumOpticsBase.mul!(drho, Ham, tmp, eltype(rho)(1), zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho, tmp, Ham, -eltype(rho)(1), one(eltype(rho)))
    return drho
end
