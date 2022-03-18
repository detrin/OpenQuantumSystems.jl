
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


"""
    thermal_state(T, mu_array, Ham, aggIndices;
    \tboltzmann_const = 0.69503476, diagonalize = false, diagonal = false)

Get initial state as thermal state excited with ultra-fast laser pulse. In this version
we suppose that after the thermal state is excited with laser pulse, the whole population
of ground state is distributed over electric states in `mu_array`. We assume
Condon approximation.

``\\rho_\\text{thermal} = \\exp( -\\frac{i}{\\hbar} H ), \\quad \\hbar = 1```.

# Arguments
* `T`: Temperature of the initial thermal state.
* `mu_array`: Vector of electric states in local basis, see [`electronicIndices`](@ref). The first
        index is for ground state of the aggregate the rest are first excited states in local basis.
* `Ham`: Arbitrary operator specifying the Hamiltonian.
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `boltzmann_const`: Boltzmann const in ``\\mathrm{cm^{-1}}``.
* `diagonalize`: Decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
        and calculate the exponential using eigenvalue decomposition.
* `diagonal`: Return only the diagonal part of the excited density matrix.
"""
function thermal_state(
    T,
    mu_array,
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    aggOperators::AggregateOperators;
    boltzmann_const::Float64 = 0.69503476,
    diagonalize::Bool = false,
    diagonal = false
)
    vibindices = aggTools.indicesMap
    Ham = aggOperators.Ham
    a1 = vibindices[1][1]
    a2 = vibindices[1][end]
    Ham_g = Ham.data[a1:a2, a1:a2]
    data = -Ham_g / (T * boltzmann_const)
    if diagonal
        data = Diagonal(data)
    end
    if !diagonalize
        data = exp(data)
    else
        H_lambda, H_S = eigen(data)
        H_lambda = diagm(H_lambda)
        data = H_S * exp(H_lambda / (T * boltzmann_const)) * inv(H_S)
    end
    W0 = DenseOperator(Ham.basis_r, Ham.basis_l, zero(Ham.data))
    excitedElInd = 0
    for elInd in 1:length(vibindices)
        excitedElInd += elInd * (mu_array[1][elInd] - 1)
    end
    a1 = vibindices[excitedElInd][1]
    a2 = vibindices[excitedElInd][end]
    W0.data[a1:a2, a1:a2] = data

    normalize!(W0)
    return W0
end

function thermal_state_old(
    T,
    mu_array,
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    aggOperators::AggregateOperators;
    boltzmann_const::Float64 = 0.69503476,
    diagonalize::Bool = false,
    diagonal = false
)
    vibindices = aggTools.indicesMap
    aggIndices = aggTools.indices
    Ham = aggOperators.Ham
    aggIndsLen = length(aggIndices)
    data = -Ham.data / (T * boltzmann_const)
    if diagonal
        data = Diagonal(data)
    end
    if !diagonalize
        data = exp(data)
    else
        H_lambda, H_S = eigen(data)
        H_lambda = diagm(H_lambda)
        data = H_S * exp(-H_lambda / (T * boltzmann_const)) * inv(H_S)
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
            W0.data[I, J] = 0.0
        end
    end
    normalize!(W0)
    return W0
end

"""
    thermal_state_composite(T, mu_weighted, Ham, aggIndices;
    \tboltzmann_const::Float64 = 0.69503476, diagonalize::Bool=false, diagonal=false)

Functionality of this method is similar to [`thermal_state`](@ref), but the final
state is constructed from partial `thermal_states` with weight specified in `mu_weighted`.
For example

```julia
thermal_state_composite(T, [0.0, 0.8, 0.2], ...) = 0.8 * thermal_state(T, [1, 2, 1], ...) + 0.2 * thermal_state(T, [1, 1, 2], ...)
````

``\\rho_\\text{thermal} = \\exp( -\\frac{i}{\\hbar} H ), \\quad \\hbar = 1``.

# Arguments
* `T`: Temperature of the initial thermal state.
* `mu_weighted`: Vector of weights of electric states in local basis, see [`electronicIndices`](@ref).
* `Ham`: Arbitrary operator specifying the Hamiltonian.
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `boltzmann_const`: Boltzmann const in ``\\mathrm{cm^{-1}}``.
* `diagonalize`: Decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
        and calculate the exponential using eigenvalue decomposition.
* `diagonal`: Return only the diagonal part of the excited density matrix.
"""
function thermal_state_composite(
    T,
    mu_weighted,
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    aggOperators::AggregateOperators;
    boltzmann_const::Float64 = 0.69503476,
    diagonalize::Bool = false,
    diagonal = false,
)
    vibindices = aggTools.indicesMap
    aggIndices = aggTools.indices
    Ham = aggOperators.Ham
    a1 = vibindices[1][1]
    a2 = vibindices[1][end]
    Ham_g = Ham.data[a1:a2, a1:a2]
    data = -Ham_g / (T * boltzmann_const)
    if diagonal
        data = Diagonal(data)
    end
    if !diagonalize
        data = exp(data)
    else
        H_lambda, H_S = eigen(data)
        H_lambda = diagm(H_lambda)
        data = H_S * exp(H_lambda / (T * boltzmann_const)) * inv(H_S)
    end
    W0 = DenseOperator(Ham.basis_r, Ham.basis_l, zero(Ham.data))

    for excitedElInd = 1:length(mu_weighted)
        a1 = vibindices[excitedElInd][1]
        a2 = vibindices[excitedElInd][end]
        W0.data[a1:a2, a1:a2] = data * mu_weighted[excitedElInd]
    end

    normalize!(W0)
    return W0
end

function thermal_state_composite_old(
    T,
    mu_weighted,
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    aggOperators::AggregateOperators;
    boltzmann_const::Float64 = 0.69503476,
    diagonalize::Bool = false,
    diagonal = false,
)
    vibindices = aggTools.indicesMap
    aggIndices = aggTools.indices
    Ham = aggOperators.Ham
    aggIndsLen = length(aggIndices)
    data = -Ham.data / (T * boltzmann_const)
    if diagonal
        data = Diagonal(data)
    end
    if !diagonalize
        data = exp(data)
    else
        H_lambda, H_S = eigen(data)
        H_lambda = diagm(H_lambda)
        data = H_S * exp(-H_lambda / (T * boltzmann_const)) * inv(H_S)
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

function ultrafast_laser_excitation(T::AbstractFloat, weights::Array, agg::Aggregate; diagonalize = true)
    molCount = agg.core.molCount
    mu_array = [ones(Int64, molCount+1)]

    mu_array_tmp = deepcopy(mu_array)
    mu_array_tmp[1][1] = 2
    W0 = weights[1] * thermal_state(T, mu_array_tmp, agg.core, agg.tools, agg.operators; diagonalize = diagonalize)

    for i in 2:molCount+1
        mu_array_tmp = deepcopy(mu_array)
        mu_array_tmp[1][i] = 2
        W0.data[:, :] += weights[i] * thermal_state(T, mu_array_tmp, agg.core, agg.tools, agg.operators; diagonalize = diagonalize).data
    end
    normalize!(W0)
    W0 = DenseOperator(W0.basis_l, W0.basis_r, complex(W0.data))

    rho0 = trace_bath(W0, agg.core, agg.tools)
    rho0 = DenseOperator(rho0.basis_l, rho0.basis_r, complex(rho0.data))

    W0_bath = get_rho_bath(W0, agg.core, agg.tools)
    W0_bath = DenseOperator(W0_bath.basis_l, W0_bath.basis_r, complex(W0_bath.data))
    return W0, rho0, W0_bath
end
