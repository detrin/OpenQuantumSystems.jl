
# include("core.jl")

struct Mode{T}
    omega::T
    shift::T
    function Mode(omega::T, shift::T) where {T<:ComputableType}
        new{T}(omega, shift)
    end
end

function franckCondonFactors(size::T, shift::U) where {T<:Integer,U<:ComputableType}
    b = GenericBasis([size + 100])
    shift_op = ShiftOperator(b, shift)
    return shift_op.data[1:size, 1:size]
end

function vibrationalIndices(maxInds::Vector{T}) where {T<:Integer}
    vibInds = Array{Array{T,1},1}(undef, 0)
    indLen = length(maxInds)
    currentInds = fill(1, (indLen))
    ind = 0::T
    push!(vibInds, copy(currentInds))
    while currentInds != maxInds
        ind += 1
        currentInds[indLen] += 1
        indPos = indLen
        while currentInds[indPos] > maxInds[indPos]
            currentInds[indPos] = (currentInds[indPos] - 1) % maxInds[indPos] + 1
            indPos -= 1
            currentInds[indPos] += 1
        end
        push!(vibInds, copy(currentInds))
    end
    return vibInds
end

function electronicIndices(molCount::T; groundState = true) where {T<:Integer}
    vibInds = Array{Array{T,1},1}(undef, 0)
    currentInds = fill(1, (molCount))
    if groundState
        push!(vibInds, copy(currentInds))
    end
    for indPos = 1:molCount
        currentInds[indPos] += 1
        push!(vibInds, copy(currentInds))
        currentInds[indPos] -= 1
    end
    return vibInds
end

mutable struct Molecule{T<:Integer,C1<:ComputableType,C2<:ComputableType}
    modes::Vector{Mode{C1}}
    Nvib::T
    fcFactors::Array{Array{C2,2},1}
    E::Array{C1,1}
    function Molecule{T,C1,C2}(
        modes::Vector{Mode{C1}},
        Nvib::T,
        E::Array{C1,1},
    ) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
        fcFactors = Array{Array{C2,2},1}(undef, 0)
        for mode in modes
            push!(fcFactors, franckCondonFactors(Nvib, mode.shift))
        end
        new(modes, Nvib, fcFactors, E)
    end
end

Molecule(modes::Vector{Mode{C}}, Nvib::T, E::Array{C,1}) where {C,T} =
    Molecule{T,C,C}(modes, Nvib, E)

function getMolStateEnergy(
    mol::Molecule{T,C1,C2},
    molElState::U,
    molVibState::Vector{U},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType,U<:Integer}
    energy = mol.E[molElState]
    for mode_i = 1:length(mol.modes)
        mode = mol.modes[mode_i]
        energy += mode.omega * (molVibState[mode_i] - 1 + 0.5)
    end
    return energy
end

function getMolStateFC(
    mol::Molecule{T,C1,C2},
    molElState1::U,
    molVibState1::Vector{U},
    molElState2::U,
    molVibState2::Vector{U},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType,U<:Integer}
    fc = 1.0::C2
    for mode_i = 1:length(mol.modes)
        fcFac = mol.fcFactors[mode_i]
        if molElState1 == 2 && molElState2 == 1
            fc *= fcFac[molVibState1[mode_i], molVibState2[mode_i]]
        elseif molElState1 == 1 && molElState2 == 2
            fc *= fcFac[molVibState2[mode_i], molVibState1[mode_i]]
        elseif molVibState1[mode_i] != molVibState2[mode_i]
            fc *= 0
        end # if molVibState1[mode_i] == molVibState2[mode_i] -> fc = 1.
    end
    return fc
end
