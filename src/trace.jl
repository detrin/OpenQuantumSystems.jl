
"""
    getFCProd(agg, FCFact, aggIndices, vibindices; groundState = false)

Get product of Franck-Condon factors. This way the trace over bath will be faster

"""
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


"""
    trace_bath(rho, agg, FCProd, aggIndices, vibindices; groundState = false)

Trace out bath degrees of freedom from `rho`

``\\rho_\\text{tr} = \\operatorname{tr}_B \\{\\rho\\} = 
\\sum_{k} \\langle k \\vert \\left( \\sum_{ab} \\rho_{am, bn} \\vert am \\rangle \\langle bn \\vert \\right)\\vert k \\rangle``

"""
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


function trace_bath(
    rho::T,
    agg,
    FCProd,
    aggIndices,
    vibindices;
    groundState = false,
) where {B<:Basis,T<:Operator{B,B}}
    rho_traced =
        trace_bath(rho.data, agg, FCProd, aggIndices, vibindices; groundState = groundState)
    basisLen = size(rho_traced, 1)
    basis = GenericBasis([basisLen])
    return DenseOperator(basis, basis, rho_traced)
end

"""
    trace_bath_slow(rho, agg, FCFact, aggIndices, vibindices; groundState = false)

Trace out bath degrees of freedom from `rho` without the product of Franck-Condon factors.

"""
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

"""
    trace_bath(rho, a, b, agg, FCProd, aggIndices, vibindices; groundState = false)

Trace out bath degrees of freedom from `rho` without the product of Franck-Condon factors.
The trace will be done only on the Hilber space for electric bra part `a` and ket part `b`.
Input density matrix `rho` is for the whole Hilber space. This method returns number.

"""
function trace_bath(rho::Array, a, b, agg, FCProd, aggIndices, vibindices)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    rho_traced = eltype(rho)(0)

    a_vibindices = vibindices[a] 
    b_vibindices = vibindices[b]

    for I in a_vibindices
        for J in b_vibindices
            rho_traced += rho[I, J] * FCProd[I, J]
        end
    end
    return rho_traced
end

function trace_bath(
    rho::T,
    a,
    b,
    agg,
    FCProd,
    aggIndices,
    vibindices
) where {B<:Basis,T<:Operator{B,B}}
    rho_traced =
        trace_bath(rho.data, a, b, agg, FCProd, aggIndices, vibindices)
    return rho_traced
end

"""
    trace_bath_part(rho, a, b, agg, FCProd, aggIndices, vibindices; groundState = false)

Trace out bath degrees of freedom from `rho` without the product of Franck-Condon factors.
The trace will be done only on the Hilber space for electric bra part `a` and ket part `b`.
Input density matrix `rho` is only for the subspace. This method returns number.

"""
function trace_bath_part(rho::Array, a, b, agg, FCProd, aggIndices, vibindices; groundState = false)
    vibLen = length(vibindices[2])
    rho_traced = eltype(rho)(0)

    a_vibindices = vibindices[a] 
    b_vibindices = vibindices[b]

    for a_vib in 1:vibLen
        I = a_vibindices[a_vib]
        for b_vib in 1:vibLen
            J = b_vibindices[b_vib]
            rho_traced += rho[a_vib, b_vib] * FCProd[I, J]
        end
    end
    return rho_traced
end

function trace_bath_part(
    rho::T,
    a,
    b,
    agg,
    FCProd,
    aggIndices,
    vibindices; 
    groundState = false
) where {B<:Basis,T<:Operator{B,B}}
    rho_traced =
        trace_bath(rho.data, a, b, agg, FCProd, aggIndices, vibindices; groundState = groundState)
    return rho_traced
end


"""
    get_rho_bath(rho, agg, FCProd, aggIndices, vibindices; groundState=false, justCopy=false)

This method will return the bath part of `rho` knowing the result of [`trace_bath`](@ref) defined as follows

`` \\rho_\\text{bath} = \\operatorname{tr}_S \\{\\rho\\} ``

`` \\rho_{\\text{bath}, ab} = \\rho_{ab} / \\langle a \\vert \\operatorname{tr}_B \\{ \\rho \\}\\vert b \\rangle``

"""
function get_rho_bath(rho::Array, agg, FCProd, aggIndices, vibindices; groundState=false, justCopy=false)
    rho_traced = trace_bath(rho, agg, FCProd, aggIndices, vibindices; groundState=groundState)
    vibLen = length(vibindices[end])
    aggIndLen = length(aggIndices)
    elLen = length(agg.molecules)
    if groundState
        elLen += 1
    end
    rho_bath = zeros(eltype(rho), aggIndLen, aggIndLen)
    rho_bath_ref = zeros(eltype(rho), vibLen, vibLen)
    if groundState
        el1_p = 0; el2_p = 0
        for el1=1:elLen, el2=1:elLen
            if abs(rho_traced[el1, el2]) != 0
                el1_p = el1; el2_p = el2
                break
            end
        end
        vib11 = vibindices[el1_p][1]; vib12 = vibindices[el1_p][end]
        vib21 = vibindices[el2_p][1]; vib22 = vibindices[el2_p][end]
        rho_bath_ref[:, :] = rho[vib11:vib12, vib21:vib22] / rho_traced[el1_p, el2_p]
        for el1=1:elLen, el2=1:elLen
            vib11 = vibindices[el1][1]; vib12 = vibindices[el1][end]
            vib21 = vibindices[el2][1]; vib22 = vibindices[el2][end]
            if abs(rho_traced[el1, el2]) != 0
                rho_bath[vib11:vib12, vib21:vib22] = rho[vib11:vib12, vib21:vib22] / rho_traced[el1, el2]
            else
                if justCopy
                    rho_bath[vib11:vib12, vib21:vib22] = rho[vib11:vib12, vib21:vib22]
                else
                    rho_bath[vib11:vib12, vib21:vib22] = rho_bath_ref[:, :]
                end
            end
        end
    else
        el1_p = 0; el2_p = 0
        for el1=1:elLen, el2=1:elLen
            if abs(rho_traced[el1, el2]) != 0
                el1_p = el1; el2_p = el2
                break
            end
        end
        el1_p += 1; el2_p += 1
        vib11 = vibindices[el1_p][1]; vib12 = vibindices[el1_p][end]
        vib21 = vibindices[el2_p][1]; vib22 = vibindices[el2_p][end]
        rho_bath_ref[:, :] = rho[vib11:vib12, vib21:vib22] / rho_traced[el1_p-1, el2_p-1]
        for el1=2:elLen+1, el2=2:elLen+1
            vib11 = vibindices[el1][1]; vib12 = vibindices[el1][end]
            vib21 = vibindices[el2][1]; vib22 = vibindices[el2][end]
            if abs(rho_traced[el1-1, el2-1]) != 0
                rho_bath[vib11:vib12, vib21:vib22] = rho[vib11:vib12, vib21:vib22] / rho_traced[el1-1, el2-1]
            else
                if justCopy
                    rho_bath[vib11:vib12, vib21:vib22] = rho[vib11:vib12, vib21:vib22]
                else
                    rho_bath[vib11:vib12, vib21:vib22] = rho_bath_ref[:, :]
                end
            end
        end
    end
    return rho_bath
end

function get_rho_bath(
        rho::T, agg, FCProd, aggIndices, vibindices; groundState=false, justCopy=false
        ) where {B<:Basis,T<:Operator{B,B}}
    rho_data = get_rho_bath(rho.data, agg, FCProd, aggIndices, vibindices; groundState=groundState, justCopy=justCopy)
    return DenseOperator(rho.basis_l, rho.basis_r, rho_data)
end


"""
    ad(rho_traced, rho_bath, agg, FCProd, aggIndices, vibindices; groundState=false)

This is the inverse operation to the trace over bath [`trace_bath`](@ref) and [`get_rho_bath`](@ref) 
defined as follows

`` \\rho = \\operatorname{ad}\\{\\rho_\\text{tr}, \\rho_\\text{bath} \\} ``

"""
function ad(rho_traced::Array, rho_bath::Array, agg, FCProd, aggIndices, vibindices; groundState=false)
    elLen = length(agg.molecules)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    if groundState
        elLen += 1
    end
    W = zero(rho_bath)

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

                W[I, J] = rho_traced[elOrder1-1, elOrder2-1] * rho_bath[I, J]
            end
        end
    else
        for I = 1:aggIndLen
            elind1, vibind1 = aggIndices[I]
            elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

            for J = 1:aggIndLen
                elind2, vibind2 = aggIndices[J]
                elOrder2 = OpenQuantumSystems.elIndOrder(elind2)

                W[I, J] = rho_traced[elOrder1, elOrder2] * rho_bath[I, J]
            end
        end
    end
    return W
end

function ad(
    rho_traced::T, 
    W_bath::Array,
    agg,
    FCFact,
    aggIndices,
    vibindices;
    groundState = false,
) where {B<:Basis,T<:Operator{B,B}}
    W = ad(rho_traced.data, W_bath, agg, FCFact, aggIndices, vibindices; groundState = groundState)
    basisLen = size(W, 1)
    basis = GenericBasis([basisLen])
    return DenseOperator(basis, basis, W)
end

function ad(
    rho_traced::Array, 
    W_bath::T,
    agg,
    FCFact,
    aggIndices,
    vibindices;
    groundState = false,
) where {B<:Basis,T<:Operator{B,B}}
    W = ad(rho_traced, W_bath.data, agg, FCFact, aggIndices, vibindices; groundState = groundState)
    basisLen = size(W, 1)
    basis = GenericBasis([basisLen])
    return DenseOperator(basis, basis, W)
end

function ad(
    rho_traced::T, 
    W_bath::T,
    agg,
    FCFact,
    aggIndices,
    vibindices;
    groundState = false,
) where {B<:Basis,T<:Operator{B,B}}
    W = ad(rho_traced.data, W_bath.data, agg, FCFact, aggIndices, vibindices; groundState = groundState)
    basisLen = size(W, 1)
    basis = GenericBasis([basisLen])
    return DenseOperator(basis, basis, W)
end


"""
    correlation_function(t, rho0_bath, Ham_0, Ham_I, agg, FCProd, aggInds, vibindices; groundState = true)

Get time dependent correlation function for a specified time `t` using following definition

`` C(t) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(t) \\hat{H}_I \\rho_\\text{bath}(0) \\} ``

"""
function correlation_function(t, rho0_bath, Ham_0, Ham_I, agg, FCProd, aggInds, vibindices; groundState = true)
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    prod = Ham_II_t * Ham_I * rho0_bath
    return trace_bath(prod, agg, FCProd, aggInds, vibindices; groundState = groundState)
end