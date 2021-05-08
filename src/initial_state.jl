
#=
function exp_series(A::Array, N::Integer)
    ret = A^0
    M = deepcopy(A)
    ret = zero(A)
    fac = 1
    for k in 1:N
        # M = M * A / k
        fac *= k
        ret += A^k / fac
    end
    return ret
end

function exp_series(A::T, N::Integer) where {B<:Basis,T<:Operator{B,B}}
    ret = exp_series(A.data, N)
    return DenseOperator(A.basis_l, A.basis_r, ret)
end
=#

function thermal_state(T, mu_array, Ham, aggIndices, boltzmann_const = 0.69503476)
    aggIndsLen = length(aggIndices)
    W0 = exp(-Ham/(T*boltzmann_const))

    for I = 1:aggIndsLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = OpenQuantumSystems.elIndOrder(elind1)
        if elind1 in mu_array
            continue
        end
        for J = 1:aggIndsLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = OpenQuantumSystems.elIndOrder(elind2)
            if elind2 in mu_array
                continue
            end
            W0.data[I, J] = 0.
        end
    end
    normalize!(W0)
    return W0
end