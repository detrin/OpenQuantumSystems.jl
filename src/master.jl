
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq

function master_int(rho0::T, tspan::Array, Ham_S::U, Ham_int::V;
        reltol::Float64=1.0e-12, abstol::Float64=1.0e-12,
        alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
        fout::Union{Function,Nothing}=nothing,
        kwargs...) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B},V<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    # (du,u,h,p,t)
    tmp = copy(rho0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) = dmaster_int(t, rho, drho, history_fun, tmp, p, rho0, Ham_S, Ham_int)
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate_delayed(tspan_, dmaster_, history_fun, x0, state, dstate, fout; reltol=reltol, abstol=abstol, alg=alg, kwargs...)
end

function dmaster_int(t::AbstractFloat, rho::T,
        drho::T, history_fun::Function, tmp::Array, p, rho0, Ham_S::U, Ham_int::V) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B},V<:Operator{B,B}}
    Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    QuantumOpticsBase.mul!(drho,Ham_II_t,rho0,-eltype(rho)(im),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,rho0,Ham_II_t,eltype(rho)(im),one(eltype(rho)))

    kernel_integrated, err = QuadGK.quadgk(s -> kernel_int(t, s, tmp, history_fun, p, Ham_II_t.data, Ham_S, Ham_int), 0, t, rtol=1e-8)
    LinearAlgebra.mul!(drho.data,-one(eltype(rho)),kernel_integrated,one(eltype(rho)),one(eltype(rho)))
    return drho
end

function kernel_int(t, s, tmp, h, p, Ham_II_t, Ham_S, Ham_int)
    Ham_II_s = getInteractionHamIPictureA(Ham_S, Ham_int, s)
    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end
    drho = deepcopy(rho)
    # commutator(Ham.data, commutator(Ham.data, tmp))
    QuantumOpticsBase.mul!(tmp,Ham_II_s,rho,eltype(rho)(1),zero(eltype(rho)))
    QuantumOpticsBase.mul!(tmp,rho,Ham_II_s,-eltype(rho)(1),one(eltype(rho)))

    QuantumOpticsBase.mul!(drho,Ham_II_t,tmp,eltype(rho)(1),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,tmp,Ham_II_t,-eltype(rho)(1),one(eltype(rho)))
    return drho
end


function master(rho0::T, tspan, Ham::U;
        reltol::Float64=1.0e-12, abstol::Float64=1.0e-12,
        alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
        fout::Union{Function,Nothing}=nothing,
        kwargs...) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    # (du,u,h,p,t)
    tmp = copy(rho0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) = dmaster(t, rho, drho, history_fun, tmp, p, rho0, Ham)
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate_delayed(tspan_, dmaster_, history_fun, x0, state, dstate, fout; reltol=reltol, abstol=abstol, alg=alg, kwargs...)
    end

function dmaster(t::AbstractFloat, rho::T,
        drho::T, history_fun::Function, tmp::Array, p, rho0, Ham::U) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B}}
    QuantumOpticsBase.mul!(drho,Ham,rho0,-eltype(rho)(im),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,rho0,Ham,eltype(rho)(im),one(eltype(rho)))

    kernel_integrated, err = QuadGK.quadgk(s -> kernel(t, s, tmp, history_fun, p, Ham.data), 0, t, rtol=1e-8)
    LinearAlgebra.mul!(drho.data,-one(eltype(rho)),kernel_integrated,one(eltype(rho)),one(eltype(rho)))
    return drho
    end

function kernel(t, s, tmp, h, p, Ham)
    rho = h(p, s)
    if (typeof(rho) <: Operator)
        rho = rho.data
    end
    drho = deepcopy(rho)
    # commutator(Ham.data, commutator(Ham.data, tmp))
    QuantumOpticsBase.mul!(tmp,Ham,rho,eltype(rho)(1),zero(eltype(rho)))
    QuantumOpticsBase.mul!(tmp,rho,Ham,-eltype(rho)(1),one(eltype(rho)))

    QuantumOpticsBase.mul!(drho,Ham,tmp,eltype(rho)(1),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,tmp,Ham,-eltype(rho)(1),one(eltype(rho)))
    return drho
end
