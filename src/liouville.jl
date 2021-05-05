
import OrdinaryDiffEq
import QuantumOpticsBase

function liouvilleVonNeumann(rho0::T, H::AbstractOperator{B,B}, tspan;
        reltol::Float64=1.0e-12, abstol::Float64=1.0e-12,
        alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.DP5(),
        fout::Union{Function,Nothing}=nothing,
        kwargs...) where {B<:Basis,T<:Operator{B,B}}
    tmp = copy(rho0)
    dliouvilleVonNeumann_(t, rho::T, drho::T) = dliouvilleVonNeumann(rho, H, drho, tmp)
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate(tspan_, dliouvilleVonNeumann_, x0, state, dstate, fout; reltol=reltol, abstol=abstol, alg=alg, kwargs...)
end

function dliouvilleVonNeumann(rho::T, H::AbstractOperator{B,B},
        drho::T, tmp::T) where {B<:Basis,T<:Operator{B,B}}
    QuantumOpticsBase.mul!(drho,H,rho,-eltype(rho)(im),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,rho,H,eltype(rho)(im),one(eltype(rho)))
    return drho
end
