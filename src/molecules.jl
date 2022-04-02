
function convert_units(E::T; from = "1/cm", to="1/fs")::T where T<:AbstractFloat
    c_const = 299792458.0
    pi_const = 3.141592653589793
    E_new = E
    if from == "1/cm" && to == "1/fs"
        c = 2.0*pi_const*c_const*1e-13
        E_new = E / c
    end
    return E_new
end

function convert_units(E_vec::Vector{T}; from = "1/cm", to="1/fs")::Vector{T} where T<:AbstractFloat
    c_const = 299792458.0
    pi_const = 3.141592653589793
    E_new = deepcopy(E_vec)
    if from == "1/cm" && to == "1/fs"
        c = 2.0*pi_const*c_const*1e-13
        E_new = map((E) -> E/c, E_vec)
    end
    return E_new
end

"""
    Mode{T}(omega, shift)

Stuct which purpose is to model vibrational LHO mode in [`Molecule`](@ref).

``V(q) = \\hbar \\omega (q - q_0)^2, \\quad \\hbar = 1``

# Arguments
* `omega`: The frequency of LHO (``\\omega``).
* `shift`: The shift of the coordinate of LHO (``q_0``).
"""
struct Mode{T}
    omega::T
    shift::T
    function Mode{T}(omega::T, shift::T) where {T<:ComputableType}
        new{T}(omega, shift)
    end
end

Mode(omega::T, shift::T) where {T<:ComputableType} = Mode{T}(omega, shift) 

function Mode(;
    omega = 200.0,
    hr_factor = 0.02
) 
    # S = d / (2 d hbar)
    shift = 2 * hr_factor
    return Mode(omega, shift)
end

Base.:(==)(x::Mode, y::Mode) = x.omega == y.omega && x.shift == y.shift

function convert_units(mode::Mode; from = "1/cm", to="1/fs")::Mode
    omega_new = convert_units(mode.omega, from=from, to=to)
    return Mode(omega_new, mode.shift)
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
electronicIndices(molCount)

Get the electric indices for all states on [`Aggregate`](@ref).

# Arguments
* `molCount`: Number of molecules in [`Aggregate`](@ref).
"""
function electronicIndices(molCount::T) where {T<:Integer}
    vibInds = Array{Array{T,1},1}(undef, 0)
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
    Molecule{T,C1,C2}(modes, Nvib, E)

Stuct which purpose is to model a molecule in [`Aggregate`](@ref).

# Arguments
* `modes`: Vector of modes ([`Mode`](@ref)).
* `Nvib`: Maximum number of vibrational states for all modes.
* `E`: Energy of ground and excited state of molecule (HOMO, LUMO).
"""
struct Molecule{T<:Integer,C1<:ComputableType,C2<:ComputableType}
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

Base.:(==)(x::Molecule, y::Molecule) = 
    x.modes == y.modes && x.Nvib == y.Nvib &&
    x.fcFactors == y.fcFactors && x.E == y.E

function convert_units(molecule::Molecule; from = "1/cm", to="1/fs")::Molecule
    modes_new = map((mode) -> convert_units(mode; from=from, to=to), molecule.modes)
    E_new = map((E_vec) -> convert_units(E_vec; from=from, to=to), molecule.E)
    return Molecule(modes_new, molecule.Nvib, E_new)
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

"""
    getMolShifts(mol)

Get shifts of every mode on the [`Molecule`](@ref) state.

# Arguments
* `mol`: Instance of [`Molecule`](@ref).
"""
function getMolShifts(
    mol::Molecule{T,C1,C2},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    shifts = Array{C1,1}(undef, 0)
    for mode in mol.modes
        push!(shifts, mode.shift)
    end
    return shifts
end

"""
    getMolFrequencies(mol)

Get shifts of every mode on the [`Molecule`](@ref) state.

# Arguments
* `mol`: Instance of [`Molecule`](@ref).
"""
function getMolFrequencies(
    mol::Molecule{T,C1,C2},
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
    frquencies = Array{C1,1}(undef, 0)
    for mode in mol.modes
        push!(frquencies, mode.omega)
    end
    return frquencies
end