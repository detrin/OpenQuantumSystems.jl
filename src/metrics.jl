
import QuantumOpticsBase:
    tracenorm, tracenorm_h, tracenorm_nh, tracedistance, tracedistance_h, tracedistance_nh
import LinearAlgebra: ishermitian

"""
    tracenorm(rho)
Trace norm of `rho`.
It is defined as
```math
T(ρ) = Tr\\{\\sqrt{ρ^† ρ}\\}.
```
Depending if `rho` is hermitian either [`tracenorm_h`](@ref) or
[`tracenorm_nh`](@ref) is called.
"""
function tracenorm(rho::DenseSuperOpType)
    ishermitian(rho) ? tracenorm_h(rho) : tracenorm_nh(rho)
end

function tracenorm!(rho::DenseSuperOpType)
    if ishermitian(rho)
        rho.data ./= tracenorm_h(rho)
    else
        rho.data ./= tracenorm_nh(rho)
    end
end

"""
    tracenorm_h(rho)
Trace norm of `rho`.
It uses the identity
```math
T(ρ) = Tr\\{\\sqrt{ρ^† ρ}\\} = \\sum_i |λ_i|
```
where ``λ_i`` are the eigenvalues of `rho`.
"""
function tracenorm_h(rho::DenseSuperOpType{B,B,T}) where {B<:Tuple{Basis,Basis},T}
    s = eigvals(Hermitian(rho.data))
    sum(abs.(s))
end

#=
function tracenorm_h(rho::T) where {T<:AbstractOperator}
    throw(ArgumentError("tracenorm_h not implemented for $(typeof(rho)). Use dense operators instead."))
end
=#

"""
    tracenorm_nh(rho)
Trace norm of `rho`.
Note that in this case `rho` doesn't have to be represented by a square
matrix (i.e. it can have different left-hand and right-hand bases).
It uses the identity
```math
    T(ρ) = Tr\\{\\sqrt{ρ^† ρ}\\} = \\sum_i σ_i
```
where ``σ_i`` are the singular values of `rho`.
"""
tracenorm_nh(rho::SuperOperator{B,B,T}) where {B<:Tuple{Basis,Basis},T} =
    sum(svdvals(rho.data))

"""
    tracedistance(rho, sigma)
Trace distance between `rho` and `sigma`.
It is defined as
```math
T(ρ,σ) = \\frac{1}{2} Tr\\{\\sqrt{(ρ - σ)^† (ρ - σ)}\\}.
```
It calls [`tracenorm`](@ref) which in turn either uses [`tracenorm_h`](@ref)
or [`tracenorm_nh`](@ref) depending if ``ρ-σ`` is hermitian or not.
"""
tracedistance(
    rho::SuperOperator{B,B,T},
    sigma::SuperOperator{B,B,T},
) where {B<:Tuple{Basis,Basis},T} = 0.5 * tracenorm(rho - sigma)

"""
    tracedistance_h(rho, sigma)
Trace distance between `rho` and `sigma`.
It uses the identity
```math
T(ρ,σ) = \\frac{1}{2} Tr\\{\\sqrt{(ρ - σ)^† (ρ - σ)}\\} = \\frac{1}{2} \\sum_i |λ_i|
```
where ``λ_i`` are the eigenvalues of `rho` - `sigma`.
"""
tracedistance_h(
    rho::SuperOperator{B,B,T},
    sigma::SuperOperator{B,B,T},
) where {B<:Tuple{Basis,Basis},T} = 0.5 * tracenorm_h(rho - sigma)


"""
    tracedistance_nh(rho, sigma)
Trace distance between `rho` and `sigma`.
Note that in this case `rho` and `sigma` don't have to be represented by square
matrices (i.e. they can have different left-hand and right-hand bases).
It uses the identity
```math
    T(ρ,σ) = \\frac{1}{2} Tr\\{\\sqrt{(ρ - σ)^† (ρ - σ)}\\}
         = \\frac{1}{2} \\sum_i σ_i
```
where ``σ_i`` are the singular values of `rho` - `sigma`.
"""
tracedistance_nh(
    rho::SuperOperator{B1,B2,T},
    sigma::SuperOperator{B1,B2,T},
) where {B1<:Tuple{Basis,Basis},B2<:Tuple{Basis,Basis},T} = 0.5 * tracenorm_nh(rho - sigma)
