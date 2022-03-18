
abstract type AbstractAggregateCore end

"""
    AggregateCore{T,C1,C2}(molecules, coupling)
    AggregateCore{T,C1,C2}(molecules)

Mutable stuct which purpose is to store important information about the whole aggregate.

# Arguments
* `molecules`: Vector of molecules ([`Molecule`](@ref)).
* `coupling`: Matrix of couplings ``J_{nm}`` between molecules.
"""
struct AggregateCore{T<:Integer,C1<:ComputableType,C2<:ComputableType} <: AbstractAggregateCore
    molecules::Vector{Molecule{T,C1,C2}}
    coupling::Matrix{C1}
    molCount::Int64
    function AggregateCore{T,C1,C2}(
        molecules::Vector{Molecule{T,C1,C2}},
        coupling::Matrix{C1},
        molCount::Int64
    ) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
        new(molecules, coupling, molCount)
    end
end

AggregateCore(molecules::Vector{Molecule{T,C1,C2}}, coupling::Matrix{C1}) where {T,C1,C2} =
    AggregateCore{T,C1,C2}(molecules, coupling, length(molecules))
AggregateCore(molecules::Vector{Molecule{T,C1,C2}}) where {T,C1,C2} =
    AggregateCore{T,C1,C2}(molecules, zeros(C1, (length(molecules) + 1, length(molecules) + 1)), length(molecules))


Base.:(==)(x::AggregateCore, y::AggregateCore) =
    x.molecules == y.molecules && x.coupling == y.coupling && x.molCount == y.molCount

"""
    getNvib(aggCore)

Get maximum number of vibrational states of each molecule in the [`Aggregate`](@ref).

"""
function getNvib(aggCore::AggregateCore)
    # TODO: add typeofs
    NvibMols = Array{Array{Int64,1},1}(undef, 0)
    for mol in aggCore.molecules
        push!(NvibMols, fill(mol.Nvib, length(mol.modes)))
    end
    return NvibMols
end

"""
    getShifts(aggCore)

Get shifts for every mode on each molecule in the [`Aggregate`](@ref).

"""
function getShifts(aggCore::AggregateCore)
    # TODO: add typeofs
    shifts = Array{Array{Float64,1},1}(undef, 0)
    for mol in aggCore.molecules
        push!(shifts, getMolShifts(mol))
    end
    return shifts
end

"""
    getFrequencies(aggCore)

Get frequencies for every mode on each molecule in the [`Aggregate`](@ref).

"""
function getFrequencies(aggCore::AggregateCore)
    # TODO: add typeofs
    frequencies = Array{Array{Float64,1},1}(undef, 0)
    for mol in aggCore.molecules
        push!(frequencies, getMolFrequencies(mol))
    end
    return frequencies
end

"""
    getAggStateEnergy(agg, aggElState, aggVibState)

Get Hamiltonian of the [`Aggregate`](@ref) in a form of DenseOperator.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggElState`: Aggregate electric state (e.g. [1, 1, 2]).
* `aggVibState`: Aggregate vibrational state (e.g. [[9, 1], [2, 5, 3], [10]]).
"""
function getAggStateEnergy(
    aggCore::AggregateCore,
    aggElState::Vector{U},
    aggVibState::Vector{Vector{U}},
)::Float64 where {U<:Integer}
    # TODO: typeof
    energy = 0.0::Float64
    for mol_i = 1:length(aggCore.molecules)
        molElState = aggElState[mol_i]
        molVibState = aggVibState[mol_i]
        energy += getMolStateEnergy(aggCore.molecules[mol_i], molElState, molVibState)
    end
    return energy
end

function elIndOrder(elInd::Vector{T})::T where {T<:Integer}
    len = length(elInd)
    ind = 1::T
    for i = 1:len
        ind += (elInd[i] - 1) * i
    end
    return ind
end
