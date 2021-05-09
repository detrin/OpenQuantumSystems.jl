
import QuantumOpticsBase:
    tracenorm, tracenorm_h, tracenorm_nh, tracedistance, tracedistance_h, tracedistance_nh
import LinearAlgebra: ishermitian


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

function tracenorm_h(rho::DenseSuperOpType{B,B,T}) where {B<:Tuple{Basis,Basis},T}
    s = eigvals(Hermitian(rho.data))
    sum(abs.(s))
end

#=
function tracenorm_h(rho::T) where {T<:AbstractOperator}
    throw(ArgumentError("tracenorm_h not implemented for $(typeof(rho)). Use dense operators instead."))
end
=#

tracenorm_nh(rho::SuperOperator{B,B,T}) where {B<:Tuple{Basis,Basis},T} =
    sum(svdvals(rho.data))

tracedistance(
    rho::SuperOperator{B,B,T},
    sigma::SuperOperator{B,B,T},
) where {B<:Tuple{Basis,Basis},T} = 0.5 * tracenorm(rho - sigma)

tracedistance_h(
    rho::SuperOperator{B,B,T},
    sigma::SuperOperator{B,B,T},
) where {B<:Tuple{Basis,Basis},T} = 0.5 * tracenorm_h(rho - sigma)

tracedistance_nh(
    rho::SuperOperator{B1,B2,T},
    sigma::SuperOperator{B1,B2,T},
) where {B1<:Tuple{Basis,Basis},B2<:Tuple{Basis,Basis},T} = 0.5 * tracenorm_nh(rho - sigma)
