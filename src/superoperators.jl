
import QuantumOpticsBase: Operator, SuperOperator, Basis
import LinearAlgebra: ishermitian

ishermitian(A::SuperOperator{B1,B2}) where {B1<:Tuple{Basis,Basis},B2<:Tuple{Basis,Basis}} = ishermitian(A.data)

function Commutator(A::Operator)::SuperOperator
    return spre(A) - spost(A)
end