



"""
    trace_bath(rho, agg, FCProd, aggIndices, vibindices)

Trace out bath degrees of freedom from `rho`

``\\rho_\\text{tr} = \\operatorname{tr}_B \\{\\rho\\} = 
\\sum_{k} \\langle k \\vert \\left( \\sum_{ab} \\rho_{am, bn} \\vert am \\rangle \\langle bn \\vert \\right)\\vert k \\rangle``

"""
function trace_bath_ground_excited(
    W::Array,
    aggCore::AggregateCore,
    aggTools::AggregateTools
)
    elLen = aggCore.molCount + 1
    rho = zeros(eltype(W), elLen, elLen)

    for I = 1:aggTools.bSize
        elind1, vibind1 = aggTools.indices[I]
        elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

        for J = 1:aggTools.bSize
            elind2, vibind2 = aggTools.indices[J]
            elOrder2 = OpenQuantumSystems.elIndOrder(elind2)

            rho[elOrder1, elOrder2] += W[I, J] * aggTools.FCproduct[I, J]
        end
    end
    return rho
end

function trace_bath_ground_ground(
    W::Array,
    aggCore::AggregateCore,
    aggTools::AggregateTools
)
    elLen = aggCore.molCount + 1
    rho = zeros(eltype(W), elLen, elLen)
    indicesMap = aggTools.indicesMap

    for a=1:length(indicesMap)
        a1 = indicesMap[a][1]
        a2 = indicesMap[a][end]
        for b=1:length(indicesMap)
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            rho[a, b] = tr(W[a1:a2, b1:b2])
        end
    end
    return rho
end

function trace_bath(
    W::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools;
    vib_basis::Symbol=:none
) 
    if vib_basis == :none
        vib_basis = aggOperators.vib_basis
    end

    if vib_basis ∉ (:ground_ground, :ground_excited)
        throw(ArgumentError("Optional argument vib_basis has to be selected from (:ground_ground, :ground_excited)"))
    end
    rho = nothing
    
    if vib_basis == :ground_excited
        rho = trace_bath_ground_excited(W, aggCore, aggTools)
    end
    if vib_basis == :ground_ground
        rho = trace_bath_ground_ground(W, aggCore, aggTools)
    end
    return rho
end

function trace_bath(
    W::Array,
    agg::Aggregate;
    vib_basis::Symbol=:none
) 
    aggCore = agg.core
    aggOperators = agg.operators
    aggTools = agg.tools
    return trace_bath(W, aggCore, aggOperators, aggTools; vib_basis=vib_basis)
end

function trace_bath(
    W::Operator,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools;
    vib_basis::Symbol=:none
) 
    rho = trace_bath(W.data, aggCore, aggOperators, aggTools; vib_basis=vib_basis)
    return DenseOperator(aggTools.basisSystem, aggTools.basisSystem, rho)
end

function trace_bath(
    W::Operator,
    agg::Aggregate;
    vib_basis::Symbol=:none
) 
    aggCore = agg.core
    aggOperators = agg.operators
    aggTools = agg.tools
    return trace_bath(W, aggCore, aggOperators, aggTools; vib_basis=vib_basis)
end

"""
    trace_bath_slow(rho, agg, FCFact, aggIndices, vibindices)

Trace out bath degrees of freedom from `rho` without the product of Franck-Condon factors.

"""
function trace_bath_slow(
    W::Array,
    aggCore::AggregateCore,
    aggTools::AggregateTools
)
    elLen = aggCore.molCount + 1
    rho = zeros(eltype(W), elLen, elLen)

    for I = 1:aggTools.bSize
        elind1, vibind1 = aggTools.indices[I]
        elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

        for J = 1:aggTools.bSize
            elind2, vibind2 = aggTools.indices[J]
            elOrder2 = OpenQuantumSystems.elIndOrder(elind2)

            for m = 1:aggTools.bBathSize
                # according to quantarhei, trace_over_vibrations()
                K = aggTools.indicesMap[1][m]
                L = aggTools.indicesMap[1][m]
                rho[elOrder1, elOrder2] += aggTools.FCfactors[K, I] * W[I, J] * aggTools.FCfactors[J, L]
            end
        end
    end
    return rho
end

function trace_bath_slow(
    W::Operator,
    aggCore::AggregateCore,
    aggTools::AggregateTools
) 
    rho = trace_bath_slow(
        W.data,
        aggCore,
        aggTools
    )
    return DenseOperator(aggTools.basisSystem, aggTools.basisSystem, rho)
end

"""
    trace_bath(rho, a, b, agg, FCProd, aggIndices, vibindices)
Trace out bath degrees of freedom from `rho` without the product of Franck-Condon factors.
The trace will be done only on the Hilber space for electric bra part `a` and ket part `b`.
Input density matrix `rho` is for the whole Hilber space. This method returns number.
"""
function trace_bath(    
    W::Array, 
    a::N, 
    b::N, 
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_excited
) where {N <: Integer}
    if vib_basis ∉ (:ground_ground, :ground_excited)
        throw(ArgumentError("Optional argument vib_basis has to be selected from (:ground_ground, :ground_excited)"))
    end
    rho = eltype(W)(0)

    if vib_basis == :ground_excited
        for I in aggTools.indicesMap[a]
            for J in aggTools.indicesMap[b]
                rho += W[I, J] * aggTools.FCproduct[I, J]
            end
        end
    end
    if vib_basis == :ground_ground
        a1 = aggTools.indicesMap[a][1]
        a2 = aggTools.indicesMap[a][end]
        b1 = aggTools.indicesMap[b][1]
        b2 = aggTools.indicesMap[b][end]
        rho = tr(W[a1:a2, b1:b2])
    end
    
    rho
end

function trace_bath(
    W::Operator,
    a::N, 
    b::N, 
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_excited
) where {N <: Integer}
    rho = trace_bath(W.data, a, b, aggTools; vib_basis=vib_basis)
    return rho
end

"""
    trace_bath_part(rho, a, b, agg, FCProd, aggIndices, vibindices)
Trace out bath degrees of freedom from `rho` without the product of Franck-Condon factors.
The trace will be done only on the Hilber space for electric bra part `a` and ket part `b`.
Input density matrix `rho` is only for the subspace. This method returns number.
"""
function trace_bath_part(
    W::Array,
    a::N, 
    b::N, 
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_ground
) where {N <: Integer}
    rho_traced = eltype(W)(0)
    if vib_basis ∉ (:ground_ground, :ground_excited)
        throw(ArgumentError("Optional argument vib_basis has to be selected from (:ground_ground, :ground_excited)"))
    end

    if vib_basis == :ground_excited
        for a_vib = 1:aggTools.bBathSize
            I = aggTools.indicesMap[a][a_vib]
            for b_vib = 1:aggTools.bBathSize
                J = aggTools.indicesMap[b][b_vib]
                rho_traced += W[a_vib, b_vib] * aggTools.FCproduct[I, J]
            end
        end
    end

    if vib_basis == :ground_ground && a == b
        rho_traced = tr(W)
    end
    
    return rho_traced
end

function trace_bath_part(
    W::Operator,
    a::N, 
    b::N, 
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_ground
) where {N <: Integer}
    rho_traced = trace_bath_part(
        W.data,
        a,
        b,
        aggTools;
        vib_basis=vib_basis
    )
    return rho_traced
end


"""
    get_rho_bath(rho, agg, FCProd, aggIndices, vibindices; justCopy=false)

This method will return the bath part of `rho` knowing the result of [`trace_bath`](@ref) defined as follows

`` \\rho_\\text{bath} = \\operatorname{tr}_S \\{\\rho\\} ``

`` \\rho_{\\text{bath}, ab} = \\rho_{ab} / \\langle a \\vert \\operatorname{tr}_B \\{ \\rho \\}\\vert b \\rangle``

"""
function get_rho_bath(
    W::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools;
    justCopy::Bool = false,
    vib_basis::Symbol=:none
)
    rho = trace_bath(W, aggCore, aggOperators, aggTools; vib_basis=vib_basis)
    vibindices = aggTools.indicesMap
    vibLen = length(vibindices[end])
    elLen = aggCore.molCount + 1
    rho_bath = zeros(eltype(rho), aggTools.bSize, aggTools.bSize)
    rho_bath_ref = zeros(eltype(rho), vibLen, vibLen)

    el1_p = 0
    el2_p = 0
    for el1 = 1:elLen, el2 = 1:elLen
        if abs(rho[el1, el2]) != 0
            el1_p = el1
            el2_p = el2
            break
        end
    end
    vib11 = vibindices[el1_p][1]
    vib12 = vibindices[el1_p][end]
    vib21 = vibindices[el2_p][1]
    vib22 = vibindices[el2_p][end]
    rho_bath_ref[:, :] = W[vib11:vib12, vib21:vib22] / rho[el1_p, el2_p]
    for el1 = 1:elLen, el2 = 1:elLen
        vib11 = vibindices[el1][1]
        vib12 = vibindices[el1][end]
        vib21 = vibindices[el2][1]
        vib22 = vibindices[el2][end]
        if abs(rho[el1, el2]) != 0
            rho_bath[vib11:vib12, vib21:vib22] =
                W[vib11:vib12, vib21:vib22] / rho[el1, el2]
        else
            if justCopy
                rho_bath[vib11:vib12, vib21:vib22] = W[vib11:vib12, vib21:vib22]
            else
                rho_bath[vib11:vib12, vib21:vib22] = rho_bath_ref[:, :]
            end
        end
    end

    return rho_bath
end

function get_rho_bath(
    W::Operator,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools;
    justCopy = false,
    vib_basis::Symbol=:none
)
    rho_data = get_rho_bath(
        W.data,
        aggCore,
        aggOperators,
        aggTools;
        justCopy = justCopy,
        vib_basis = vib_basis
    )
    return DenseOperator(W.basis_l, W.basis_r, rho_data)
end


"""
    ad(rho, rho_bath, agg, FCProd, aggIndices, vibindices)

This is the inverse operation to the trace over bath [`trace_bath`](@ref) and [`get_rho_bath`](@ref) 
defined as follows

`` \\rho = \\operatorname{ad}\\{\\rho_\\text{tr}, \\rho_\\text{bath} \\} ``

"""
function ad(
    rho::Array,
    rho_bath::Array,
    aggCore::AggregateCore,
    aggTools::AggregateTools
)
    W = zero(rho_bath)

    for I = 1:aggTools.bSize
        elind1, vibind1 = aggTools.indices[I]
        elOrder1 = OpenQuantumSystems.elIndOrder(elind1)

        for J = 1:aggTools.bSize
            elind2, vibind2 = aggTools.indices[J]
            elOrder2 = OpenQuantumSystems.elIndOrder(elind2)

            W[I, J] = rho[elOrder1, elOrder2] * rho_bath[I, J]
        end
    end
    return W
end

function ad(
    rho::Operator,
    W_bath::Array,
    aggCore::AggregateCore,
    aggTools::AggregateTools;
) 
    W = ad(
        rho.data,
        W_bath,
        aggCore,
        aggTools
    )
    return DenseOperator(aggTools.basis, aggTools.basis, W)
end

function ad(
    rho::Array,
    W_bath::Operator,
    aggCore::AggregateCore,
    aggTools::AggregateTools
) 
    W = ad(
        rho,
        W_bath.data,
        aggCore,
        aggTools
    )
    return DenseOperator(aggTools.basis, aggTools.basis, W)
end

function ad(
    rho::Operator,
    W_bath::Operator,
    aggCore::AggregateCore,
    aggTools::AggregateTools;
) 
    W = ad(
        rho.data,
        W_bath.data,
        aggCore,
        aggTools
    )
    return DenseOperator(aggTools.basis, aggTools.basis, W)
end


"""
    correlation_function(t, rho0_bath, Ham_0, Ham_I, agg, FCProd, aggInds, vibindices)

Get time dependent correlation function for a specified time `t` using following definition

`` C(t) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(t) \\hat{H}_I \\rho_\\text{bath}(0) \\} ``

"""
function correlation_function(
    t,
    rho0_bath,
    Ham_0,
    Ham_I,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools;
    vib_basis::Symbol=:none
)
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    prod = Ham_II_t * Ham_I * rho0_bath
    if vib_basis == :none
        vib_basis = aggOperators.vib_basis
    end
    return trace_bath(prod, aggCore, aggOperators, aggTools; vib_basis=vib_basis)
end
