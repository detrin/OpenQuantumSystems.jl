
function getFCProd(agg, FCFact, aggIndices, vibindices; groundState = false)
    elLen = length(agg.molecules)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    if groundState
        elLen += 1
    end
    FCProd = zeros(eltype(FCFact), aggIndLen, aggIndLen)

    if !groundState
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)
            if elOrder1 == 1
                continue
            end

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)
                if elOrder2 == 1
                    continue
                end

                for m = 1:vibLen
                    K = vibindices[elOrder1][m]
                    L = vibindices[elOrder2][m]
                    FCProd[I, J] += FCFact[K, I] * FCFact[J, L]
                end
            end
        end
    else
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)
                #=
                K1 = vibindices[elOrder1][1]
                K2 = vibindices[elOrder1][end]
                L1 = vibindices[elOrder2][1]
                L2 = vibindices[elOrder2][end]
                FCProd[I, J] = sum(FCFact[K1:K2, I] .* FCFact[J, L1:L2])
                =#
                
                for m = 1:vibLen
                    K = vibindices[elOrder1][m]
                    L = vibindices[elOrder2][m]
                    FCProd[I, J] += FCFact[K, I] * FCFact[J, L]
                end
            end
        end
    end
    return FCProd
end

function trace_bath(rho::Array, agg, FCProd, aggIndices, vibindices; groundState = false)
    elLen = length(agg.molecules)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    if groundState
        elLen += 1
    end
    rho_traced = zeros(eltype(rho), elLen, elLen)

    if !groundState
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)
            if elOrder1 == 1
                continue
            end

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)
                if elOrder2 == 1
                    continue
                end

                rho_traced[elOrder1-1, elOrder2-1] += rho[I, J] * FCProd[I, J]
            end
        end
    else
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)

                rho_traced[elOrder1, elOrder2] += rho[I, J] * FCProd[I, J]
            end
        end
    end
    return rho_traced
end


function trace_bath_slow(rho::Array, agg, FCFact, aggIndices, vibindices; groundState = false)
    elLen = length(agg.molecules)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    if groundState
        elLen += 1
    end
    rho_traced = zeros(eltype(rho), elLen, elLen)

    if !groundState
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)
            if elOrder1 == 1
                continue
            end

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)
                if elOrder2 == 1
                    continue
                end

                for m = 1:vibLen
                    K = vibindices[elOrder1][m]
                    L = vibindices[elOrder2][m]
                    rho_traced[elOrder1-1, elOrder2-1] +=
                        FCFact[K, I] * rho[I, J] * FCFact[J, L]
                end
            end
        end
    else
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)

                for m = 1:vibLen
                    K = vibindices[elOrder1][m]
                    L = vibindices[elOrder2][m]
                    rho_traced[elOrder1, elOrder2] +=
                        FCFact[K, I] * rho[I, J] * FCFact[J, L]
                end
            end
        end
    end
    return rho_traced
end

function trace_bath(
    rho::T,
    agg,
    FCFact,
    aggIndices,
    vibindices;
    groundState = false,
) where {B<:Basis,T<:Operator{B,B}}
    rho_traced =
        trace_bath(rho.data, agg, FCFact, aggIndices, vibindices; groundState = groundState)
    basisLen = size(rho_traced, 1)
    basis = GenericBasis([basisLen])
    return DenseOperator(basis, basis, rho_traced)
end

function trace_bath_slow(
    rho::T,
    agg,
    FCFact,
    aggIndices,
    vibindices;
    groundState = false,
) where {B<:Basis,T<:Operator{B,B}}
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggIndices, vibindices; groundState = groundState)
    basisLen = size(rho_traced, 1)
    basis = GenericBasis([basisLen])
    return DenseOperator(basis, basis, rho_traced)
end