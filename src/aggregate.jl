
# include("core.jl")

"""
    Aggregate{T,C1,C2}(molecules, coupling)
    Aggregate{T,C1,C2}(molecules)

Mutable stuct which purpose is to store important information about the whole system.

# Arguments
* `molecules`: Vector of molecules ([`Molecule`](@ref)).
* `coupling`: Matrix of couplings ``J_{nm}`` between molecules.
"""
mutable struct Aggregate{T<:Integer,C1<:ComputableType,C2<:ComputableType}
    molecules::Vector{Molecule{T,C1,C2}}
    coupling::Matrix{C1}
    function Aggregate{T,C1,C2}(
        molecules::Vector{Molecule{T,C1,C2}},
        coupling::Matrix{C1},
    ) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
        new(molecules, coupling)
    end
end

Aggregate(molecules::Vector{Molecule{T,C1,C2}}, coupling::Matrix{C1}) where {T,C1,C2} =
    Aggregate{T,C1,C2}(molecules, coupling)
Aggregate(molecules::Vector{Molecule{T,C1,C2}}) where {T,C1,C2} =
    Aggregate{T,C1,C2}(molecules, zeros(C1, (length(molecules) + 1, length(molecules) + 1)))

"""
    getNvib(agg)

Get maximum number of vibrational states of each molecule in the [`Aggregate`](@ref).

"""
function getNvib(
    agg::Aggregate{T,C1,C2},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    NvibMols = Array{Array{T,1},1}(undef, 0)
    for mol in agg.molecules
        push!(NvibMols, fill(mol.Nvib, length(mol.modes)))
    end
    return NvibMols
end

"""
    getShifts(agg)

Get shifts for every mode on each molecule in the [`Aggregate`](@ref).

"""
function getShifts(
    agg::Aggregate{T,C1,C2},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    shifts = Array{Array{C1,1},1}(undef, 0)
    for mol in agg.molecules
        push!(shifts, getMolShifts(mol))
    end
    return shifts
end

"""
    getFrequencies(agg)

Get frequencies for every mode on each molecule in the [`Aggregate`](@ref).

"""
function getFrequencies(
    agg::Aggregate{T,C1,C2},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    frequencies = Array{Array{C1,1},1}(undef, 0)
    for mol in agg.molecules
        push!(frequencies, getMolFrequencies(mol))
    end
    return frequencies
end

"""
    vibrationalIndices(agg)

Get all vibrational indices of the [`Aggregate`](@ref).

"""
function vibrationalIndices(
    agg::Aggregate{T,C1,C2},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    NvibMols = getNvib(agg)
    molLen = length(agg.molecules)
    NvibIndMols = Array{Array{Array{T,1},1},1}(undef, 0)
    for mol_i = 1:molLen
        push!(NvibIndMols, vibrationalIndices(NvibMols[mol_i]))
    end
    molNvib = map((vibInds) -> length(vibInds), NvibIndMols)
    molInds = vibrationalIndices(molNvib)
    aggInds = Array{Array{Array{T,1},1},1}(undef, 0)
    for molInd in molInds
        aggInd = Array{Array{T,1},1}(undef, 0)
        for mol_i = 1:molLen
            push!(aggInd, NvibIndMols[mol_i][molInd[mol_i]])
        end
        push!(aggInds, aggInd)
    end
    return aggInds
end

"""
    vibrationalIndices(agg)

Get all electronic indices of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function electronicIndices(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    vibInds = Array{Array{T,1},1}(undef, 0)
    molCount = length(agg.molecules)
    currentInds = fill(1, (molCount))
    push!(vibInds, copy(currentInds))
    for indPos = 1:molCount
        currentInds[indPos] += 1
        push!(vibInds, copy(currentInds))
        currentInds[indPos] -= 1
    end
    return vibInds
end

"""
    getIndices(agg)

Get all indices (electronic and vibrational) of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function getIndices(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    vibInds = vibrationalIndices(agg)
    vibIndsLen = length(vibInds)
    elInds = electronicIndices(agg)
    elIndsLen = length(elInds)
    indices = Array{Array{Array{T,1} where T,1},1}(undef, 0)
    for el_i = 1:elIndsLen
        for vib_i = 1:vibIndsLen
            push!(indices, [elInds[el_i], vibInds[vib_i]])
        end
    end
    return indices
end

"""
    getVibIndices(agg)

Get pointers (integers) to the indices of the [`Aggregate`](@ref) separated by 
electronic states (e.g. [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function getVibIndices(agg, aggIndices)
    aggIndLen = length(aggIndices)
    vibindices = Array{Array{Int,1},1}(undef, 0)
    elLen = length(agg.molecules) + 1
    for el_i = 1:elLen
        push!(vibindices, Array{Int,1}(undef, 0))
    end

    for I = 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = OpenQuantumSystems.elIndOrder(elind1)
        push!(vibindices[elOrder1], I)
    end
    return vibindices
end

"""
    getFranckCondonFactors(agg, aggIndices)
    getFranckCondonFactors(agg)

Get Frack-Condon factors of the [`Aggregate`](@ref) in a form of matrix.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices generated by [`getIndices`](@ref).
"""
function getFranckCondonFactors(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    FC = zeros(C2, (aggIndLen, aggIndLen))
    for I = 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        for J = 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            if elind1 == elind2
                if vibind1 == vibind2
                    FC[I, J] = 1.0
                end
            else
                fc = 1.0::C2
                for mi = 1:molLen
                    mol = agg.molecules[mi]
                    fc *=
                        getMolStateFC(mol, elind1[mi], vibind1[mi], elind2[mi], vibind2[mi])
                end
                FC[I, J] = fc
            end
        end
    end
    return FC
end

getFranckCondonFactors(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} =
    getFranckCondonFactors(agg::Aggregate{T,C1,C2}, nothing)

"""
    getAggStateEnergy(agg, aggElState, aggVibState)

Get Hamiltonian of the [`Aggregate`](@ref) in a form of DenseOperator.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggElState`: Aggregate electric state (e.g. [1, 1, 2]).
* `aggVibState`: Aggregate vibrational state (e.g. [[9, 1], [2, 5, 3], [10]]).
"""
function getAggStateEnergy(
    agg::Aggregate{T,C1,C2},
    aggElState::Vector{U},
    aggVibState::Vector{Vector{U}},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType,U<:Integer}
    energy = 0.0::C1
    for mol_i = 1:length(agg.molecules)
        molElState = aggElState[mol_i]
        molVibState = aggVibState[mol_i]
        energy += getMolStateEnergy(agg.molecules[mol_i], molElState, molVibState)
    end
    return energy
end

function elIndOrder(elInd::Vector{T}) where {T<:Integer}
    len = length(elInd)
    ind = 1::T
    for i = 1:len
        ind += (elInd[i] - 1) * i
    end
    return ind
end

"""
    getAggHamSystemSmall(agg)

Get Hamiltonian of the system, ``H_S`` only for ``\\mathcal{H}_S``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function getAggHamSystemSmall(
    agg::Aggregate{T,C1,C2};
    groundEnergy::Bool = true,
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    molLen = length(agg.molecules)
    Ham_sys = zeros(C2, (molLen + 1, molLen + 1))

    agg_shifts = getShifts(agg)
    agg_frequencies = getShifts(agg)
    reorganisation_energies = []
    for mol_i = 1:length(agg.molecules)
        mol = agg.molecules[mol_i]
        reorganisation_energy = 0.
        for mode_i = 1:length(mol.modes)
            reorganisation_energy += agg_frequencies[mol_i][mode_i] * agg_shifts[mol_i][mode_i]^2 / 2.
        end
        push!(reorganisation_energies, reorganisation_energy)
    end

    elInds = electronicIndices(agg)
    E_agg = zeros(C1, (2, molLen))
    E_agg[1, :] = map(mol -> mol.E[1], agg.molecules)
    E_agg[2, :] = map(mol -> mol.E[2], agg.molecules)
    for elInd in elInds
        ind = OpenQuantumSystems.elIndOrder(elInd)
        E_state = 0
        for mol_i = 1:molLen
            E_state += E_agg[elInd[mol_i], mol_i]
        end
        # add reorganisation energy
        if ind > 1
            E_state += reorganisation_energies[ind-1]
        end
        Ham_sys[ind, ind] = E_state 
    end

    Ham_sys[:, :] += agg.coupling[:, :]
    if !groundEnergy
        E0 = Ham_sys[1, 1]
        for i = 1:size(Ham_sys, 1)
            Ham_sys[i, i] -= E0
        end
    end
    b = GenericBasis([size(Ham_sys, 1)])
    return DenseOperator(b, b, Ham_sys)
end

"""
    getAggHamSystemBig(agg)

Get Hamiltonian of the system, ``H_S`` for ``\\mathcal{H}_S \\otimes \\mathcal{H}_B``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function getAggHamSystemBig(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any,
    franckCondonFactors::Any;
    groundEnergy::Bool = true,
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham_sys = getAggHamSystemSmall(agg; groundEnergy = groundEnergy)
    Ham = zeros(C2, (aggIndLen, aggIndLen))
    for I = 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = elIndOrder(elind1)
        for J = 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = elIndOrder(elind2)

            Ham[I, J] = Ham_sys.data[elOrder1, elOrder2] * franckCondonFactors[I, J]
        end
    end
    if !groundEnergy
        E0 = Ham[1, 1]
        for I = 1:aggIndLen
            Ham[I, I] -= E0
        end
    end
    b = GenericBasis([aggIndLen])
    return DenseOperator(b, b, Ham)
end

"""
    getAggHamBathSmall(agg)

Get Hamiltonian of the bath, ``H_B`` only for ``\\mathcal{H}_B``.

"""
function getAggHamBathSmall(
    agg::Aggregate{T,C1,C2};
    groundEnergy::Bool = true
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    molLen = length(agg.molecules)
    vibInds = vibrationalIndices(agg)
    vibIndsLen = length(vibInds)
    Ham_bath = zeros(C2, (vibIndsLen, vibIndsLen))

    agg2 = deepcopy(agg)
    for mol_i in molLen
        agg2.molecules[mol_i].E[1] = 0.0
    end
    elind = fill(1, (molLen + 1))
    Ham_sys = getAggHamSystemSmall(agg)
    for I = 1:vibIndsLen
        vibind = vibInds[I]
        Ham_bath[I, I] = getAggStateEnergy(agg, elind, vibind) - Ham_sys.data[1, 1]
    end
    b = GenericBasis([vibIndsLen])
    if !groundEnergy
        E0 = Ham_bath[1, 1]
        for i = 1:size(Ham_bath, 1)
            Ham_bath[i, i] -= E0
        end
    end
    return DenseOperator(b, b, Ham_bath)
end

"""
    getAggHamBathBig(agg)

Get Hamiltonian of the bath, ``H_B`` only for ``\\mathcal{H}_S \\otimes \\mathcal{H}_B``.

"""
function getAggHamBathBig(
    agg::Aggregate{T,C1,C2};
    groundEnergy::Bool = true
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    aggIndices = getIndices(agg)
    aggIndLen = length(aggIndices)
    b = GenericBasis([aggIndLen])

    Ham_bath = getAggHamBathSmall(agg)
    Ham_sys = getAggHamSystemSmall(agg)
    b_sys = GenericBasis([size(Ham_sys, 1)])
    one_sys = OneDenseOperator(b_sys)
    # TODO: optimise types in future
    one_sys = DenseOperator(b_sys, b_sys, real(one_sys.data))
    Ham_bath = tensor(Ham_bath, one_sys)
    if !groundEnergy
        E0 = Ham_bath.data[1, 1]
        for I = 1:size(Ham_bath.data, 1)
            Ham_bath.data[I, I] -= E0
        end
    end
    return DenseOperator(b, b, Ham_bath.data)
end

"""
    getAggHamSystemBath(agg, aggIndices, 
    \tfranckCondonFactors; groundEnergy = false)
    getAggHamSystemBath(agg, aggIndices; groundEnergy = false)
    getAggHamSystemBath(agg; groundEnergy = false)

Get Hamiltonian of the [`Aggregate`](@ref) in a form of DenseOperator.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices generated by [`getIndices`](@ref).
* `franckCondonFactors`: Franck-Condon factors generated by [`getFranckCondonFactors`](@ref).
"""
function getAggHamSystemBath(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any,
    franckCondonFactors::Any;
    groundEnergy::Bool = false,
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    
    Ham_B = getAggHamBathBig(agg; groundEnergy=true)
    Ham_S = getAggHamSystemBig(agg, aggIndices, franckCondonFactors; groundEnergy=true)
    Ham_0 = Ham_B + Ham_S
    if !groundEnergy
        E0 = Ham_0.data[1, 1]
        for I = 1:aggIndLen
            Ham_0.data[I, I] -= E0
        end
    end
    return Ham_0
end

getAggHamSystemBath(
    agg::Aggregate{T,C1,C2};
    groundEnergy::Bool = false,
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamSystemBath(
    agg::Aggregate{T,C1,C2},
    nothing,
    nothing;
    groundEnergy = groundEnergy,
)

getAggHamSystemBath(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any;
    groundEnergy::Bool = false,
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamSystemBath(
    agg::Aggregate{T,C1,C2},
    aggIndices,
    nothing;
    groundEnergy = groundEnergy,
)

"""
    getAggHamInteraction(agg, aggIndices, franckCondonFactors; 
    \t, groundEnergy = false)
    getAggHamInteraction(agg, aggIndices; groundEnergy = false)
    getAggHamInteraction(agg; groundEnergy = false)

Get interation Hamiltonian of the [`Aggregate`](@ref), ``H_I`` by substracting 
Hamiltonian and Hamiltonian of the systems and of the bath, ``H - H_S - H_B``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices generated by [`getIndices`](@ref).
* `franckCondonFactors`: Franck-Condon factors generated by [`getFranckCondonFactors`](@ref).
"""
function getAggHamInteraction(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any,
    franckCondonFactors::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham_I = zeros(C2, (aggIndLen, aggIndLen))
    agg_shifts = getShifts(agg)
    agg_frequencies = getShifts(agg)
    agg_coeffs = deepcopy(agg_frequencies)
    for mol_i = 1:length(agg.molecules)
        mol = agg.molecules[mol_i]
        for mode_i = 1:length(mol.modes)
            agg_coeffs[mol_i][mode_i] = agg_shifts[mol_i][mode_i] * sqrt(agg_frequencies[mol_i][mode_i] / 2)
        end
    end

    for I = 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = elIndOrder(elind1)
        if elOrder1 == 1
            continue
        end
        for J = 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = elIndOrder(elind2)
            if elOrder2 == 1 || elOrder1 != elOrder2
                continue
            end
            diff_num = 0; mol_j = 0; mode_j = 0
            for mol_i = 1:length(agg.molecules)
                mol = agg.molecules[mol_i]
                for mode_i = 1:length(mol.modes)
                    diff_vib = abs(vibind1[mol_i][mode_i] - vibind2[mol_i][mode_i])
                    diff_num += diff_vib
                    if diff_vib == 1 && diff_num == 1
                        mol_j = mol_i; mode_j = mode_i;
                    end
                end
            end

            if diff_num != 1 || mol_j != elOrder1-1
                continue
            end
            vib_n = vibind1[mol_j][mode_j]
            vib_m = vibind2[mol_j][mode_j]
            coeff = agg_coeffs[mol_j][mode_j]
            Ham_I[I, J] = - coeff * sqrt(min(vib_n, vib_m))
        end
    end
    b = GenericBasis([aggIndLen])
    return DenseOperator(b, b, Ham_I)
end

getAggHamInteraction(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamInteraction(
    agg::Aggregate{T,C1,C2},
    nothing,
    nothing
)

getAggHamInteraction(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamInteraction(
    agg::Aggregate{T,C1,C2},
    aggIndices,
    nothing
)

"""
    getAggHamiltonian(agg, aggIndices, franckCondonFactors; 
    \t, groundEnergy = false)
    getAggHamiltonian(agg, aggIndices; groundEnergy = false)
    getAggHamiltonian(agg; groundEnergy = false)

Get Hamiltonian of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices generated by [`getIndices`](@ref).
* `franckCondonFactors`: Franck-Condon factors generated by [`getFranckCondonFactors`](@ref).
"""
function getAggHamiltonian(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any,
    franckCondonFactors::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham_I = getAggHamInteraction(agg, aggIndices, franckCondonFactors)
    Ham_0 = getAggHamSystemBath(agg, aggIndices, franckCondonFactors)
    return Ham_0 + Ham_I 
end

getAggHamiltonian(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamiltonian(
    agg::Aggregate{T,C1,C2},
    nothing,
    nothing
)

getAggHamiltonian(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamiltonian(
    agg::Aggregate{T,C1,C2},
    aggIndices,
    nothing
)

### Sparse versions

"""
    getFranckCondonFactorsSparse(agg, aggIndices)
    getFranckCondonFactorsSparse(agg)

Sparse version of [`getFranckCondonFactors`](@ref).

"""
function getFranckCondonFactorsSparse(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    FC = spzeros(C2, aggIndLen, aggIndLen)
    for I = 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        for J = 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            if elind1 == elind2
                if vibind1 == vibind2
                    FC[I, J] = 1.0
                end
            else
                fc = 1.0::C2
                for mi = 1:molLen
                    mol = agg.molecules[mi]
                    fc *=
                        getMolStateFC(mol, elind1[mi], vibind1[mi], elind2[mi], vibind2[mi])
                end
                if fc != 0.0
                    FC[I, J] = fc
                end
            end
        end
    end
    return FC
end

getFranckCondonFactorsSparse(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getFranckCondonFactorsSparse(
    agg::Aggregate{T,C1,C2},
    nothing
)

"""
    getAggHamiltonianSparse(agg, aggIndices)
    getAggHamiltonianSparse(agg)

Sparse version of [`getAggHamiltonian`](@ref).

"""
function getAggHamiltonianSparse(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any,
    franckCondonFactors::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    if aggIndices === nothing
        aggIndices = getIndices(agg)
    end
    aggIndLen = length(aggIndices)
    molLen = length(agg.molecules)
    if franckCondonFactors === nothing
        franckCondonFactors = getFranckCondonFactors(agg, aggIndices)
    end
    Ham = spzeros(C2, aggIndLen, aggIndLen)
    for I = 1:aggIndLen
        elind1, vibind1 = aggIndices[I]
        elOrder1 = elIndOrder(elind1)
        for J = 1:aggIndLen
            elind2, vibind2 = aggIndices[J]
            elOrder2 = elIndOrder(elind2)
            if I == J
                val = getAggStateEnergy(agg, elind1, vibind1)
                if val != 0.0
                    Ham[I, J] = val
                end
            else
                if elind1 != elind2
                    val = agg.coupling[elOrder1, elOrder2] * franckCondonFactors[I, J]
                    if val != 0.0
                        Ham[I, J] = val
                    end
                end
            end
        end
    end
    E0 = Ham[1, 1]
    for I = 1:aggIndLen
        Ham[I, I] -= E0
    end
    b = GenericBasis([aggIndLen])
    return SparseOperator(b, b, Ham)
end

getAggHamiltonianSparse(
    agg::Aggregate{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamiltonianSparse(
    agg::Aggregate{T,C1,C2},
    nothing,
    nothing
)

getAggHamiltonianSparse(
    agg::Aggregate{T,C1,C2},
    aggIndices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamiltonianSparse(
    agg::Aggregate{T,C1,C2},
    aggIndices,
    nothing
)

"""
    setupAggregate(agg; groundEnergy=true, verbose=false)

Generate all basic data from the [`Aggregate`](@ref). Returns
`aggInds, vibindices, aggIndLen, basis, FCFact, FCProd, Ham, Ham_0, Ham_I`.

"""
function setupAggregate(agg; groundEnergy = true, verbose = false)
    if verbose
        println("aggInds")
        @time aggInds = getIndices(agg)
        println("vibindices")
        @time vibindices = getVibIndices(agg, aggInds)
        aggIndLen = length(aggInds)
        println("aggIndLen: ", aggIndLen)
        println("basis")
        @time basis = GenericBasis([aggIndLen])
        println("FCFact")
        @time FCFact = getFranckCondonFactors(agg, aggInds)
        println("FCProd")
        @time FCProd =
            getFCProd(agg, FCFact, aggInds, vibindices)
        println("Ham")
        @time Ham = getAggHamiltonian(
            agg,
            aggInds,
            FCFact;
            groundEnergy = groundEnergy,
        )
        println("Ham_0")
        @time Ham_0 = getAggHamSysBath(
            agg,
            aggInds;
            groundEnergy = groundEnergy,
        )
        println("Ham_I")
        @time Ham_I = Ham - Ham_0
    else
        aggInds = getIndices(agg)
        vibindices = getVibIndices(agg, aggInds)
        aggIndLen = length(aggInds)
        basis = GenericBasis([aggIndLen])
        FCFact = getFranckCondonFactors(agg, aggInds)
        FCProd = getFCProd(agg, FCFact, aggInds, vibindices)
        Ham = getAggHamiltonian(
            agg,
            aggInds,
            FCFact;
            groundEnergy = groundEnergy,
        )
        Ham_0 = getAggHamSysBath(
            agg,
            aggInds;
            groundEnergy = groundEnergy,
        )
        Ham_I = Ham - Ham_0
    end
    return (aggInds, vibindices, aggIndLen, basis, FCFact, FCProd, Ham, Ham_0, Ham_I)
end
