

"""
    Mode{T}(omega, shift)

Mutable stuct which purpose is to model vibrational LHO mode in [`Molecule`](@ref).

``V(q) = \\hbar \\omega (q - q_0)^2``

# Arguments
* `omega`: The frequency of LHO (``\\omega``).
* `shift`: The shift of the coordinate of LHO (``q_0``).
"""
mutable struct Mode{T}
    omega::T
    shift::T
    function Mode(omega::T, shift::T) where {T<:ComputableType}
        new{T}(omega, shift)
    end
end

"""
    franckCondonFactors(size, shift)

Get Franck-Condon factors for LHO mode calculated using [`ShiftOperator`](@ref).

"""
function franckCondonFactors(size::T, shift::U) where {T<:Integer,U<:ComputableType}
    b = GenericBasis([size + 100])
    shift_op = ShiftOperator(b, shift)
    return shift_op.data[1:size, 1:size]
end

"""
    vibrationalIndices(maxInds)

Get the vibrational indices for all states on [`Molecule`](@ref).

# Arguments
* `maxInds`: Vector of maximum number of vibrational states.
"""
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


"""
electronicIndices(molCount; groundState = true)

Get the electric indices for all states on [`Aggregate`](@ref).

# Arguments
* `molCount`: Number of molecules in [`Aggregate`](@ref).
* `groundState`: Option for allowing the ground electric state in local basis.
"""
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

"""
    Molecule{T,C1,C2}(modes, Nvib, E)

Mutable stuct which purpose is to model a molecule in [`Aggregate`](@ref).

# Arguments
* `modes`: Vector of modes ([`Mode`](@ref)).
* `Nvib`: Maximum number of vibrational states for all modes.
* `E`: Energy of ground and excited state of molecule (HOMO, LUMO).
"""
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

function updateMolecule!(mol::Molecule{T, C1, C2}) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    mol = Molecule{T, C1, C2}(mol.modes, mol.Nvib, mol.E)
end

"""
    updateMolecule(mol)

Get updated molecule with new params (e.g. [`Mode`](@ref)).
"""
function updateMolecule(mol::Molecule{T, C1, C2}) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    Molecule{T, C1, C2}(mol.modes, mol.Nvib, mol.E)
end

"""
    getMolStateEnergy(mol, molElState, molVibState)

Get the energy of the [`Molecule`](@ref) state.

# Arguments
* `mol`: Instance of [`Molecule`](@ref).
* `molElState`: Electric state of the molecule in local basis (i.e. 1 or 2, where 1 is ground state).
* `molVibState`: Vibrational state of the molecule (e.g. [1, 5, 2, 2]).
"""
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

"""
    getMolStateFC(mol, molElState1, molVibState1, molElState2, molVibState2)

Get the energy of the [`Molecule`](@ref) state.

# Arguments
* `mol`: Instance of [`Molecule`](@ref).
* `molElState1`, `molElState2`: Electric state of the molecule in local basis (i.e. 1 or 2, where 1 is ground state).
* `molVibState1`, `molVibState2`: Vibrational state of the molecule (e.g. [1, 5, 2, 2]).
"""
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
