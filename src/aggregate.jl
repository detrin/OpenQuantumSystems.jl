include("core.jl")

mutable struct Aggregate{T<:Integer, C1<:ComputableType, C2<:ComputableType}
    molecules::Vector{Molecule{T, C1, C2}}
    coupling::Matrix{C1}
    function Aggregate{T,C1,C2}(
        molecules::Vector{Molecule{T, C1, C2}},
        coupling::Matrix{C1}
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
        new(molecules, coupling)
    end
end

Aggregate(molecules::Vector{Molecule{T, C1, C2}}, coupling::Matrix{C1}) where {T, C1, C2} = Aggregate{T, C1, C2}(molecules, coupling)
Aggregate(molecules::Vector{Molecule{T, C1, C2}}) where {T, C1, C2} = Aggregate{T, C1, C2}(molecules, zeros(C1, (length(molecules)+1, length(molecules)+1)))

function getNvib(agg::Aggregate{T, C1, C2}) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    NvibMols = Array{Array{T, 1}, 1}(undef, 0)
    for mol in agg.molecules
        push!(NvibMols, fill(mol.Nvib, length(mol.modes)))
    end
    return NvibMols
end

function vibrationalIndices(agg::Aggregate{T, C1, C2}) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    NvibMols = getNvib(agg)
    molLen = length(agg.molecules)
    NvibIndMols = Array{Array{Array{T, 1}, 1}, 1}(undef, 0)
    for mol_i in 1:molLen
        push!(NvibIndMols, vibrationalIndices(NvibMols[mol_i]))
    end
    molNvib = map((vibInds) -> length(vibInds), NvibIndMols)
    molInds = vibrationalIndices(molNvib)
    aggInds = Array{Array{Array{T, 1}, 1}, 1}(undef, 0)
    for molInd in molInds
        aggInd = Array{Array{T, 1}, 1}(undef, 0)
        for mol_i in 1:molLen
            push!(aggInd, NvibIndMols[mol_i][molInd[mol_i]])
        end
        push!(aggInds, aggInd)
    end
    return aggInds
end

function electronicIndices(agg::Aggregate{T, C1, C2}; groundState=true) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    vibInds = Array{Array{T, 1}, 1}(undef, 0)
    molCount = length(agg.molecules)
    currentInds = fill(1, (molCount))
    if groundState
        push!(vibInds, copy(currentInds))
    end
    for indPos in 1:molCount
        currentInds[indPos] += 1
        push!(vibInds, copy(currentInds))
        currentInds[indPos] -= 1
    end
    return vibInds
end

function getIndices(agg::Aggregate{T, C1, C2}; groundState=true) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    vibInds = vibrationalIndices(agg)
    vibIndsLen = length(vibInds)
    elInds = electronicIndices(agg; groundState=groundState)
    elIndsLen = length(elInds)
    indices = Array{Array{Array{T,1} where T,1}, 1}(undef, 0)
    for el_i in 1:elIndsLen
        for vib_i in 1:vibIndsLen
            push!(indices, [elInds[el_i], vibInds[vib_i]])
        end
    end
    return indices
end

function getFranckCondonFactors(
        agg::Aggregate{T, C1, C2},
        aggIndices::Any; 
        groundState::Bool=true, 
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg; groundState=groundState)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    FC = zeros(C2, (aggIndLen, aggIndLen))
    for I in 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        for J in 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            if elind1 == elind2 
                if vibind1 == vibind2
                    FC[I, J] = 1.
                end
            else
                fc = 1.::C2
                for mi in 1:molLen
                    mol = agg.molecules[mi]
                    fc *= getMolStateFC(mol, elind1[mi], vibind1[mi], elind2[mi], vibind2[mi])
                end
                FC[I, J] = fc
            end
        end
    end
    return FC
end

getFranckCondonFactors(
    agg::Aggregate{T, C1, C2};
    groundState::Bool=true
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType
    } = getFranckCondonFactors(agg::Aggregate{T, C1, C2}, nothing; groundState=groundState)

function getAggStateEnergy(
        agg::Aggregate{T, C1, C2}, 
        aggElState::Vector{U}, aggVibState::Vector{Vector{U}}
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType, U<:Integer}
    energy = 0.::C1
    for mol_i in 1:length(agg.molecules)
        molElState = aggElState[mol_i]
        molVibState = aggVibState[mol_i]
        energy += getMolStateEnergy(agg.molecules[mol_i], molElState, molVibState)
    end
    return energy
end

function elIndOrder(elInd::Vector{T}) where T<:Integer
    len = length(elInd)
    ind = 1::T
    for i in 1:len
        ind += (elInd[i] - 1) * i
    end
    return ind
end

function getAggHamiltonian(
        agg::Aggregate{T, C1, C2},
        aggIndices::Any,
        franckCondonFactors::Any; 
        groundState::Bool=false,
        groundEnergy::Bool=false
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg; groundState=groundState)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham = zeros(C2, (aggIndLen, aggIndLen))
    for I in 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = elIndOrder(elind1)
        for J in 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = elIndOrder(elind2)
            if I == J
                Ham[I, J] = getAggStateEnergy(agg, elind1, vibind1)
            else
                if elind1 != elind2
                    Ham[I, J] = agg.coupling[elOrder1, elOrder2] * franckCondonFactors[I, J]
                end
            end
        end
    end
    if !groundEnergy
        E0 = Ham[1, 1]
        for I in 1:aggIndLen
            Ham[I, I] -= E0
        end
    end
    b = GenericBasis([aggIndLen])
    return DenseOperator(b, b, Ham)
end

getAggHamiltonian(
    agg::Aggregate{T, C1, C2};
    groundState::Bool=true,
    groundEnergy::Bool=false
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType
    } = getAggHamiltonian(agg::Aggregate{T, C1, C2}, nothing, nothing; groundState=groundState, groundEnergy=groundEnergy)

getAggHamiltonian(
    agg::Aggregate{T, C1, C2},
    aggIndices::Any;
    groundState::Bool=true,
    groundEnergy::Bool=false
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType
    } = getAggHamiltonian(agg::Aggregate{T, C1, C2}, aggIndices, nothing; groundState=groundState, groundEnergy=groundEnergy)

function getAggHamiltonianSystem(
        agg::Aggregate{T, C1, C2};
        groundState::Bool=false
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    molLen = length(agg.molecules)
    if groundState
        Ham_sys = zeros(C2, (molLen+1, molLen+1))
    else
        Ham_sys = zeros(C2, (molLen, molLen))
    end
    elInds = electronicIndices(agg; groundState=groundState)
    E_agg = zeros(C1, (2, molLen))
    E_agg[1, :] = map(mol -> mol.E[1], agg.molecules)
    E_agg[2, :] = map(mol -> mol.E[2], agg.molecules)
    for elInd in elInds
        ind = OpenQuantumSystems.elIndOrder(elInd)
        if !groundState
            ind -= 1
        end
        E_state = 0
        for mol_i in 1:molLen
            E_state += E_agg[elInd[mol_i], mol_i]
        end
        Ham_sys[ind, ind] = E_state
    end
    
    if groundState
        Ham_sys[:, :] += agg.coupling[:, :]
    else
        Ham_sys[:, :] += agg.coupling[2:molLen+1, 2:molLen+1]
    end
    b = GenericBasis([size(Ham_sys, 1)])
    return DenseOperator(b, b, Ham_sys)
end

function getAggHamiltonianBath(
        agg::Aggregate{T, C1, C2}
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    molLen = length(agg.molecules)
    vibInds = vibrationalIndices(agg)
    vibIndsLen = length(vibInds)
    Ham_bath = zeros(C2, (vibIndsLen, vibIndsLen))

    agg2 = deepcopy(agg)
    for mol_i in molLen
        agg2.molecules[mol_i].E[1] = 0.
    end
    elind = fill(1, (molLen+1))
    for I in 1:vibIndsLen
        vibind = vibInds[I]
        Ham_bath[I, I] = getAggStateEnergy(agg, elind, vibind)
    end
    b = GenericBasis([vibIndsLen])
    return DenseOperator(b, b, Ham_bath)
end

function getAggHamiltonianInteraction(
        agg::Aggregate{T, C1, C2},
        aggIndices::Any,
        franckCondonFactors::Any; 
        groundState::Bool=false
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg; groundState=groundState)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham = getAggHamiltonian(agg, aggIndices, franckCondonFactors; groundState=groundState, groundEnergy=true)
    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundState=groundState)
    b = GenericBasis([aggIndLen])
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_S = tensor(OneDenseOperator(b_bath), Ham_sys) + tensor(Ham_bath, OneDenseOperator(b_sys))
    H_int = Ham.data - Ham_S.data
    return DenseOperator(b, b, H_int)
end

### Sparse versions

function getFranckCondonFactorsSparse(
        agg::Aggregate{T, C1, C2},
        aggIndices::Any; 
        groundState::Bool=true, 
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg; groundState=groundState)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    FC = spzeros(C2, aggIndLen, aggIndLen)
    for I in 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        for J in 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            if elind1 == elind2 
                if vibind1 == vibind2
                    FC[I, J] = 1.
                end
            else
                fc = 1.::C2
                for mi in 1:molLen
                    mol = agg.molecules[mi]
                    fc *= getMolStateFC(mol, elind1[mi], vibind1[mi], elind2[mi], vibind2[mi])
                end
                if fc != 0.
                    FC[I, J] = fc
                end
            end
        end
    end
    return FC
end

getFranckCondonFactorsSparse(
    agg::Aggregate{T, C1, C2};
    groundState::Bool=true
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType
    } = getFranckCondonFactorsSparse(agg::Aggregate{T, C1, C2}, nothing; groundState=groundState)

function getAggHamiltonianSparse(
        agg::Aggregate{T, C1, C2},
        aggIndices::Any,
        franckCondonFactors::Any; 
        groundState::Bool=false
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg; groundState=groundState)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham = spzeros(C2, aggIndLen, aggIndLen)
    for I in 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = elIndOrder(elind1)
        for J in 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = elIndOrder(elind2)
            if I == J
                val = getAggStateEnergy(agg, elind1, vibind1)
                if val != 0.
                    Ham[I, J] = val
                end
            else
                if elind1 != elind2
                    val = agg.coupling[elOrder1, elOrder2] * franckCondonFactors[I, J]
                    if val != 0.
                        Ham[I, J] = val
                    end
                end
            end
        end
    end
    E0 = Ham[1, 1]
    for I in 1:aggIndLen
        Ham[I, I] -= E0
    end
    b = GenericBasis([aggIndLen])
    return SparseOperator(b, b, Ham)
end

getAggHamiltonianSparse(
    agg::Aggregate{T, C1, C2};
    groundState::Bool=true
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType
    } = getAggHamiltonianSparse(agg::Aggregate{T, C1, C2}, nothing, nothing; groundState=groundState)

getAggHamiltonianSparse(
    agg::Aggregate{T, C1, C2},
    aggIndices::Any;
    groundState::Bool=true
    ) where {T<:Integer, C1<:ComputableType, C2<:ComputableType
    } = getAggHamiltonianSparse(agg::Aggregate{T, C1, C2}, aggIndices, nothing; groundState=groundState)