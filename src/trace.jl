
function trace_bath(rho::Array, agg, FCFact, aggIndices, vibindices; groundState = false)
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
