
abstract type AbstractAggregateOperators end

struct MyExceptionTree <: Exception
    msg::String
end

"""
getAggHamSystemSmall(agg)

Get Hamiltonian of the system, ``H_S`` only for ``\\mathcal{H}_S``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function getAggHamSystemSmall(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_excited,
    groundEnergy::Bool = true,
)
    Ham_sys = zeros(Float64, (aggCore.molCount + 1, aggCore.molCount + 1))

    agg_shifts = getShifts(aggCore)
    agg_frequencies = getFrequencies(aggCore)
    reorganisation_energies = []
    for mol_i = 1:aggCore.molCount
        mol = aggCore.molecules[mol_i]
        reorganisation_energy = 0.0
        for mode_i = 1:length(mol.modes)
            reorganisation_energy +=
                agg_frequencies[mol_i][mode_i] * agg_shifts[mol_i][mode_i]^2 / 2.0
        end
        push!(reorganisation_energies, reorganisation_energy)
    end

    E_agg = zeros(Float64, (2, aggCore.molCount))
    E_agg[1, :] = map(mol -> mol.E[1], aggCore.molecules)
    E_agg[2, :] = map(mol -> mol.E[2], aggCore.molecules)
    for elInd in aggTools.elIndices
        ind = OpenQuantumSystems.elIndOrder(elInd)
        E_state = 0
        for mol_i = 1:aggCore.molCount
            E_state += E_agg[elInd[mol_i], mol_i]
        end
        if vib_basis == :ground_ground
            # add reorganisation energy
            if ind > 1
                E_state += reorganisation_energies[ind-1]
            end
        end
        Ham_sys[ind, ind] = E_state
    end

    Ham_sys[:, :] += aggCore.coupling[:, :]
    if !groundEnergy
        E0 = Ham_sys[1, 1]
        for i = 1:aggTools.bSystemSize
            Ham_sys[i, i] -= E0
        end
    end

    return DenseOperator(aggTools.basisSystem, aggTools.basisSystem, Ham_sys)
end

"""
getAggHamSystemBig(agg)

Get Hamiltonian of the system, ``H_S`` for ``\\mathcal{H}_S \\otimes \\mathcal{H}_B``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function getAggHamSystemBig(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_excited,
    groundEnergy::Bool = true,
)
    Ham_sys = getAggHamSystemSmall(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = groundEnergy)
    Ham = zeros(Float64, (aggTools.bSize, aggTools.bSize))
    for I = 1:aggTools.bSize
        elind1, vibind1 = aggTools.indices[I]
        elOrder1 = elIndOrder(elind1)
        for J = 1:aggTools.bSize
            elind2, vibind2 = aggTools.indices[J]
            elOrder2 = elIndOrder(elind2)

            # Ham[I, J] = Ham_sys.data[elOrder1, elOrder2] * aggTools.FCfactors[I, J]
            if vibind1 == vibind2
                Ham[I, J] = Ham_sys.data[elOrder1, elOrder2] 
            end
        end
    end
    if !groundEnergy
        E0 = Ham[1, 1]
        for I = 1:aggTools.bSize
            Ham[I, I] -= E0
        end
    end

    return DenseOperator(aggTools.basis, aggTools.basis, Ham)
end

"""
getAggHamBathSmall(agg)

Get Hamiltonian of the bath, ``H_B`` only for ``\\mathcal{H}_B``.

"""
function getAggHamBathSmall(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_excited,
    groundEnergy::Bool = true,
)
    Ham_bath = zeros(Float64, (aggTools.bBathSize, aggTools.bBathSize))

    elind = fill(1, (aggCore.molCount + 1))
    Ham_sys = getAggHamSystemSmall(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = true)
    for I = 1:aggTools.bBathSize
        vibind = aggTools.vibIndices[I]
        Ham_bath[I, I] = getAggStateEnergy(aggCore, elind, vibind) - Ham_sys.data[1, 1]
    end

    if !groundEnergy
        E0 = Ham_bath[1, 1]
        for i = 1:aggTools.bBathSize
            Ham_bath[i, i] -= E0
        end
    end
    return DenseOperator(aggTools.basisBath, aggTools.basisBath, Ham_bath)
end

"""
getAggHamBathBig(agg)

Get Hamiltonian of the bath, ``H_B`` only for ``\\mathcal{H}_S \\otimes \\mathcal{H}_B``.

"""
function getAggHamBathBig(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    groundEnergy::Bool = true,
)

    Ham_bath = getAggHamBathSmall(aggCore, aggTools)
    one_sys = OneDenseOperator(aggTools.basisSystem)
    # TODO: optimise types in future
    one_sys = DenseOperator(aggTools.basisSystem, aggTools.basisSystem, real(one_sys.data))
    Ham_bath = tensor(Ham_bath, one_sys)
    if !groundEnergy
        E0 = Ham_bath.data[1, 1]
        for I = 1:aggTools.bSize
            Ham_bath.data[I, I] -= E0
        end
    end
    return DenseOperator(aggTools.basis, aggTools.basis, Ham_bath.data)
end

"""
getAggHamSystemBath(agg, aggTools.indices, 
\tfranckCondonFactors; groundEnergy = false)
getAggHamSystemBath(agg, aggTools.indices; groundEnergy = false)
getAggHamSystemBath(agg; groundEnergy = false)

Get system and bath Hamiltonian of the [`Aggregate`](@ref) in a form of DenseOperator.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggTools.indices`: Aggregate indices generated by [`getIndices`](@ref).
* `franckCondonFactors`: Franck-Condon factors generated by [`getFranckCondonFactors`](@ref).
"""
function getAggHamSystemBath(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::Symbol=:ground_excited,
    groundEnergy::Bool = true,
)
    Ham_B = getAggHamBathBig(aggCore, aggTools; groundEnergy = true)
    Ham_S = getAggHamSystemBig(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = true)
    Ham_0 = Ham_B + Ham_S
    if !groundEnergy
        E0 = Ham_0.data[1, 1]
        for I = 1:aggTools.bSize
            Ham_0.data[I, I] -= E0
        end
    end
    return Ham_0
end


"""
getAggHamInteraction(agg, aggTools.indices, franckCondonFactors)
getAggHamInteraction(agg, aggTools.indices)
getAggHamInteraction(agg)

Get interation Hamiltonian of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggTools.indices`: Aggregate indices generated by [`getIndices`](@ref).
* `franckCondonFactors`: Franck-Condon factors generated by [`getFranckCondonFactors`](@ref).
"""
function getAggHamInteraction(aggCore::AggregateCore, aggTools::AggregateTools; vib_basis::Symbol=:ground_excited)
    if vib_basis ∉ (:ground_ground, :ground_excited)
        throw(ArgumentError("Optional argument vib_basis has to be selected from (:ground_ground, :ground_excited)"))
    end
    Ham_sys = getAggHamSystemSmall(aggCore, aggTools; vib_basis=vib_basis, groundEnergy = true)

    Ham_I = zeros(Float64, (aggTools.bSize, aggTools.bSize))
    agg_shifts = getShifts(aggCore)
    agg_frequencies = getFrequencies(aggCore)

    for I = 1:aggTools.bSize
        elind1, vibind1 = aggTools.indices[I]
        elOrder1 = elIndOrder(elind1)
        if elOrder1 == 1
            continue
        end
        for J = 1:aggTools.bSize
            elind2, vibind2 = aggTools.indices[J]
            elOrder2 = elIndOrder(elind2)
            if elOrder2 == 1 
                continue
            end
            if elOrder1 != elOrder2
                if vib_basis == :ground_excited
                    Ham_I[I, J] = Ham_sys.data[elOrder1, elOrder2] * aggTools.FCfactors[I, J]
                    if vibind1 == vibind2
                        Ham_I[I, J] -= Ham_sys.data[elOrder1, elOrder2] 
                    end
                end
            else
                # permits only one vib change
                # check all moleciles and all modes, if they differ in max one vib state
                # calculate the values of \Delta V
                #=
                diff_num = 0 # in molecules
                mol_j = 0
                mode_j = 0
                for mol_i = 1:aggCore.molCount
                    mol = aggCore.molecules[mol_i]
                    for mode_i = 1:length(mol.modes)
                        diff_vib = abs(vibind1[mol_i][mode_i] - vibind2[mol_i][mode_i])
                        diff_num += diff_vib
                        if diff_vib == 1 && diff_num == 1
                            mol_j = mol_i
                            mode_j = mode_i
                        end
                    end
                end

                if diff_num != 1 || mol_j != elOrder1 - 1
                    continue
                end
                vib_n = vibind1[mol_j][mode_j]
                vib_m = vibind2[mol_j][mode_j]
                coeff = agg_coeffs[mol_j][mode_j]
                Ham_I[I, J] = -coeff * sqrt(min(vib_n, vib_m))
                =#
                # permits any molecule change
                #=
                coeff = 1.
                mol_i = elOrder1 - 1
                mol = aggCore.molecules[mol_i]
                for mode_i = 1:length(mol.modes)
                    diff_vib = abs(vibind1[mol_i][mode_i] - vibind2[mol_i][mode_i])
                    if diff_vib == 1
                        vib_n = vibind1[mol_i][mode_i]
                        vib_m = vibind2[mol_i][mode_i]
                        coeff *= agg_coeffs[mol_i][mode_i] * sqrt(min(vib_n, vib_m))
                    end
                end
                if coeff != 1.
                    Ham_I[I, J] = - coeff 
                end
                =#
                # somehing smells here, eigenvals of Ham are not the same in both vib bath bases
                # total hours spent: 11
                # What are the correct elements of ΔV?
                # Battle not with monsters, lest ye become a monster, and if you gaze into the abyss, the abyss gazes also into you.
                #      --Nietzsche
                if vib_basis == :ground_ground
                    coeff = 0.
                    diff_vib_sum = 0
                    for mol_i = 1:aggCore.molCount
                        mol = aggCore.molecules[mol_i]
                        for mode_i = 1:length(mol.modes)
                            diff_vib = abs(vibind1[mol_i][mode_i] - vibind2[mol_i][mode_i])
                            diff_vib_sum += diff_vib
                            if diff_vib == 1 && mol_i == elOrder1 - 1
                                vib_n = vibind1[mol_i][mode_i]
                                vib_m = vibind2[mol_i][mode_i]
                                coeff += agg_frequencies[mol_i][mode_i] * agg_shifts[mol_i][mode_i] / Float64(min(vib_n, vib_m))^0.5 / 2.0
                                # println(coeff, " ", vibind1, " ", vibind2, " ", mol_i, " ", mode_i, " ", min(vib_n, vib_m))
                            end # agg_frequencies[mol_i][mode_i] * agg_shifts[mol_i][mode_i] 
                        end
                    end
                    if diff_vib_sum == 1
                        Ham_I[I, J] = coeff  # -
                    end
                end
            end
        end
    end

    return DenseOperator(aggTools.basis, aggTools.basis, Ham_I)
end

"""
getAggHamiltonian(agg, aggTools.indices, franckCondonFactors; 
\t, groundEnergy = true)
getAggHamiltonian(agg, aggTools.indices; groundEnergy = true)
getAggHamiltonian(agg; groundEnergy = true)

Get Hamiltonian of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggTools.indices`: Aggregate indices generated by [`getIndices`](@ref).
* `franckCondonFactors`: Franck-Condon factors generated by [`getFranckCondonFactors`](@ref).
"""
function getAggHamiltonian(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    groundEnergy::Bool = true,
    vib_basis::Symbol=:ground_excited
)
    Ham_I = getAggHamInteraction(aggCore, aggTools; vib_basis = vib_basis)
    Ham_0 = getAggHamSystemBath(aggCore, aggTools; groundEnergy = groundEnergy)
    return Ham_0 + Ham_I
end

### Sparse versions
#=
"""
getAggHamiltonianSparse(agg, aggTools.indices)
getAggHamiltonianSparse(agg)

Sparse version of [`getAggHamiltonian`](@ref).

"""
function getAggHamiltonianSparse(
agg::AggregateCore{T,C1,C2},
aggTools.indices::Any,
franckCondonFactors::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
if aggTools.indices === nothing
    aggTools.indices = getIndices(agg)
end
aggTools.bSize = length(aggTools.indices)
aggCore.molCount = length(agg.molecules)
if franckCondonFactors === nothing
    franckCondonFactors = getFranckCondonFactors(agg, aggTools.indices)
end
Ham = spzeros(C2, aggTools.bSize, aggTools.bSize)
for I = 1:aggTools.bSize
    elind1, vibind1 = aggTools.indices[I]
    elOrder1 = elIndOrder(elind1)
    for J = 1:aggTools.bSize
        elind2, vibind2 = aggTools.indices[J]
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
for I = 1:aggTools.bSize
    Ham[I, I] -= E0
end
b = GenericBasis([aggTools.bSize])
return SparseOperator(b, b, Ham)
end

getAggHamiltonianSparse(
agg::AggregateCore{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamiltonianSparse(
agg::AggregateCore{T,C1,C2},
nothing,
nothing
)

getAggHamiltonianSparse(
agg::AggregateCore{T,C1,C2},
aggTools.indices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = getAggHamiltonianSparse(
agg::AggregateCore{T,C1,C2},
aggTools.indices,
nothing
)
=#

# TODO: add docs and proper types
struct AggregateOperators{O_sys,O_bath,O} <: AbstractAggregateOperators
    Ham_sys::O_sys
    Ham_bath::O_bath
    Ham_S::O
    Ham_B::O
    Ham_0::O
    Ham_I::O
    Ham::O
    groundEnergy::Bool
    vib_basis::Symbol
    function AggregateOperators{O_sys,O_bath,O}(
        Ham_sys::O_sys,
        Ham_bath::O_bath,
        Ham_S::O,
        Ham_B::O,
        Ham_0::O,
        Ham_I::O,
        Ham::O,
        groundEnergy::Bool,
        vib_basis::Symbol,
    ) where {O_sys<:Operator,O_bath<:Operator,O<:Operator}
        new(Ham_sys, Ham_bath, Ham_S, Ham_B, Ham_0, Ham_I, Ham, groundEnergy, vib_basis)
    end
end

function AggregateOperators(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    groundEnergy::Bool = true,
    vib_basis::Symbol=:ground_excited
)

    Ham_sys = getAggHamSystemSmall(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = true)
    Ham_bath = getAggHamBathSmall(aggCore, aggTools; groundEnergy = true)
    Ham_S = getAggHamSystemBig(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = groundEnergy)
    Ham_B = getAggHamBathBig(aggCore, aggTools; groundEnergy = groundEnergy)
    Ham_0 = getAggHamSystemBath(aggCore, aggTools; groundEnergy = groundEnergy)
    Ham_I = getAggHamInteraction(aggCore, aggTools; vib_basis = vib_basis)
    Ham = getAggHamiltonian(aggCore, aggTools; groundEnergy = groundEnergy, vib_basis = vib_basis)

    O_sys = typeof(Ham_sys)
    O_bath = typeof(Ham_bath)
    O = typeof(Ham_S)

    return AggregateOperators{O_sys,O_bath,O}(
        Ham_sys,
        Ham_bath,
        Ham_S,
        Ham_B,
        Ham_0,
        Ham_I,
        Ham,
        groundEnergy,
        vib_basis
    )
end

Base.:(==)(x::AggregateOperators, y::AggregateOperators) =
    x.Ham_sys == y.Ham_sys &&
    x.Ham_bath == y.Ham_bath &&
    x.Ham_S == y.Ham_S &&
    x.Ham_B == y.Ham_B &&
    x.Ham_0 == y.Ham_0 &&
    x.Ham_I == y.Ham_I &&
    x.Ham == y.Ham &&
    x.groundEnergy == y.groundEnergy &&
    x.vib_basis == y.vib_basis
