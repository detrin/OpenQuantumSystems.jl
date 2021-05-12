
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

function thermal_state(T, mu_array, Ham, aggIndices; boltzmann_const::Float64 = 0.69503476, diagonalize::Bool=false, diagonal=false)
    aggIndsLen = length(aggIndices)
    data = -Ham.data/(T*boltzmann_const)
    if diagonal
        data = Diagonal(data)
    end
    if !diagonalize
        data = exp(data)
    else
        H_lambda, H_S = eigen(data)
        H_lambda = diagm(H_lambda)
        data = H_S * exp(-H_lambda/(T*boltzmann_const)) * inv(H_S)
    end
    W0 = DenseOperator(Ham.basis_r, Ham.basis_l, data)

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

function thermal_state_composite(T, mu_weighted, Ham, aggIndices; boltzmann_const::Float64 = 0.69503476, diagonalize::Bool=false, diagonal=false)
    aggIndsLen = length(aggIndices)
    data = -Ham.data/(T*boltzmann_const)
    if diagonal
        data = Diagonal(data)
    end
    if !diagonalize
        data = exp(data)
    else
        H_lambda, H_S = eigen(data)
        H_lambda = diagm(H_lambda)
        data = H_S * exp(-H_lambda/(T*boltzmann_const)) * inv(H_S)
    end
    W0 = DenseOperator(Ham.basis_r, Ham.basis_l, zero(data))

    for I = 1:aggIndsLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = OpenQuantumSystems.elIndOrder(elind1)
        mu_w = mu_weighted[elOrder1]
        for J = 1:aggIndsLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = OpenQuantumSystems.elIndOrder(elind2)
            if elOrder1 != elOrder2
                continue
            end
            W0.data[I, J] = data[I, J] * mu_w
        end
    end
    normalize!(W0)
    return W0
end