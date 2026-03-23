
abstract type AbstractAggregateOperators end

struct MyExceptionTree <: Exception
    msg::String
end

"""
get_agg_ham_system_small(agg)

Get Hamiltonian of the system, ``H_S`` only for ``\\mathcal{H}_S``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function get_agg_ham_system_small(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::VibBasisLike=GroundGround(),
    groundEnergy::Bool = true,
)
    vb = _to_vib_basis(vib_basis)
    return _get_agg_ham_system_small(aggCore, aggTools, vb; groundEnergy=groundEnergy)
end

function _get_agg_ham_system_small(
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    vb::AbstractVibBasis;
    groundEnergy::Bool = true,
)
    Ham_sys = zeros(Float64, (aggCore.molCount + 1, aggCore.molCount + 1))

    agg_shifts = get_shifts(aggCore)
    agg_frequencies = get_frequencies(aggCore)
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

    E_agg = zeros(Float64, (N_ELECTRONIC_LEVELS, aggCore.molCount))
    E_agg[ELECTRONIC_GROUND, :] = map(mol -> mol.E[ELECTRONIC_GROUND], aggCore.molecules)
    E_agg[ELECTRONIC_EXCITED, :] = map(mol -> mol.E[ELECTRONIC_EXCITED], aggCore.molecules)
    for elInd in aggTools.elIndices
        ind = elIndOrder(elInd)
        E_state = 0
        for mol_i = 1:aggCore.molCount
            E_state += E_agg[elInd[mol_i], mol_i]
        end
        E_state = _apply_reorganisation(E_state, ind, reorganisation_energies, vb)
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

_apply_reorganisation(E_state, ind, reorganisation_energies, ::GroundGround) =
    ind > 1 ? E_state + reorganisation_energies[ind-1] : E_state
_apply_reorganisation(E_state, ind, reorganisation_energies, ::GroundExcited) = E_state

"""
get_agg_ham_system_big(agg)

Get Hamiltonian of the system, ``H_S`` for ``\\mathcal{H}_S \\otimes \\mathcal{H}_B``.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
"""
function get_agg_ham_system_big(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::VibBasisLike=GroundGround(),
    groundEnergy::Bool = true,
)
    Ham_sys = get_agg_ham_system_small(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = groundEnergy)
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
get_agg_ham_bath_small(agg)

Get Hamiltonian of the bath, ``H_B`` only for ``\\mathcal{H}_B``.

"""
function get_agg_ham_bath_small(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::VibBasisLike=GroundGround(),
    groundEnergy::Bool = true,
)
    Ham_bath = zeros(Float64, (aggTools.bBathSize, aggTools.bBathSize))

    elind = fill(1, (aggCore.molCount + 1))
    Ham_sys = get_agg_ham_system_small(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = true)
    for I = 1:aggTools.bBathSize
        vibind = aggTools.vibIndices[I]
        Ham_bath[I, I] = get_agg_state_energy(aggCore, elind, vibind) - Ham_sys.data[1, 1]
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
get_agg_ham_bath_big(agg)

Get Hamiltonian of the bath, ``H_B`` only for ``\\mathcal{H}_S \\otimes \\mathcal{H}_B``.

"""
function get_agg_ham_bath_big(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    groundEnergy::Bool = true,
)

    Ham_bath = get_agg_ham_bath_small(aggCore, aggTools)
    one_sys = OneDenseOperator(aggTools.basisSystem)
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
get_agg_ham_system_bath(agg, aggTools.indices, 
\tfranck_condon_factors; groundEnergy = false)
get_agg_ham_system_bath(agg, aggTools.indices; groundEnergy = false)
get_agg_ham_system_bath(agg; groundEnergy = false)

Get system and bath Hamiltonian of the [`Aggregate`](@ref) in a form of DenseOperator.

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggTools.indices`: Aggregate indices generated by [`get_indices`](@ref).
* `franck_condon_factors`: Franck-Condon factors generated by [`get_franck_condon_factors`](@ref).
"""
function get_agg_ham_system_bath(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    vib_basis::VibBasisLike=GroundGround(),
    groundEnergy::Bool = true,
)
    Ham_B = get_agg_ham_bath_big(aggCore, aggTools; groundEnergy = true)
    Ham_S = get_agg_ham_system_big(aggCore, aggTools; vib_basis = vib_basis, groundEnergy = true)
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
get_agg_ham_interaction(agg, aggTools.indices, franck_condon_factors)
get_agg_ham_interaction(agg, aggTools.indices)
get_agg_ham_interaction(agg)

Get interation Hamiltonian of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggTools.indices`: Aggregate indices generated by [`get_indices`](@ref).
* `franck_condon_factors`: Franck-Condon factors generated by [`get_franck_condon_factors`](@ref).
"""
function _ham_interaction_off_diagonal!(Ham_I, Ham_sys, aggTools, I, J, elOrder1, elOrder2, vibind1, vibind2)
    Ham_I[I, J] = Ham_sys.data[elOrder1, elOrder2] * aggTools.FCfactors[I, J]
    if vibind1 == vibind2
        Ham_I[I, J] -= Ham_sys.data[elOrder1, elOrder2]
    end
end

function _ham_interaction_ground_ground!(Ham_I, I, J, elOrder1, vibind1, vibind2, aggCore, agg_frequencies, agg_shifts)
    coeff = 0.0
    diff_vib_sum = 0
    mol_i = elOrder1 - 1
    mol = aggCore.molecules[mol_i]
    for mode_i = 1:length(mol.modes)
        diff_vib = abs(vibind1[mol_i][mode_i] - vibind2[mol_i][mode_i])
        diff_vib_sum += diff_vib
        if diff_vib == 1
            vib_n = vibind1[mol_i][mode_i]
            vib_m = vibind2[mol_i][mode_i]
            coeff += agg_frequencies[mol_i][mode_i] * agg_shifts[mol_i][mode_i] * Float64(min(vib_n, vib_m))^0.5 / sqrt(2.0)
        end
    end
    if diff_vib_sum == 1
        Ham_I[I, J] = coeff
    end
end

function get_agg_ham_interaction(aggCore::AggregateCore, aggTools::AggregateTools; vib_basis::VibBasisLike=GroundGround())
    vb = _to_vib_basis(vib_basis)
    Ham_sys = get_agg_ham_system_small(aggCore, aggTools; vib_basis=vb, groundEnergy = true)

    Ham_I = zeros(Float64, (aggTools.bSize, aggTools.bSize))
    agg_shifts = get_shifts(aggCore)
    agg_frequencies = get_frequencies(aggCore)

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
                _ham_interaction_off_diagonal_dispatch!(Ham_I, Ham_sys, aggTools, I, J, elOrder1, elOrder2, vibind1, vibind2, vb)
                continue
            end
            _ham_interaction_diagonal_dispatch!(Ham_I, I, J, elOrder1, vibind1, vibind2, aggCore, agg_frequencies, agg_shifts, vb)
        end
    end

    return DenseOperator(aggTools.basis, aggTools.basis, Ham_I)
end

_ham_interaction_off_diagonal_dispatch!(Ham_I, Ham_sys, aggTools, I, J, elOrder1, elOrder2, vibind1, vibind2, ::GroundExcited) =
    _ham_interaction_off_diagonal!(Ham_I, Ham_sys, aggTools, I, J, elOrder1, elOrder2, vibind1, vibind2)
_ham_interaction_off_diagonal_dispatch!(Ham_I, Ham_sys, aggTools, I, J, elOrder1, elOrder2, vibind1, vibind2, ::GroundGround) = nothing

_ham_interaction_diagonal_dispatch!(Ham_I, I, J, elOrder1, vibind1, vibind2, aggCore, agg_frequencies, agg_shifts, ::GroundGround) =
    _ham_interaction_ground_ground!(Ham_I, I, J, elOrder1, vibind1, vibind2, aggCore, agg_frequencies, agg_shifts)
_ham_interaction_diagonal_dispatch!(Ham_I, I, J, elOrder1, vibind1, vibind2, aggCore, agg_frequencies, agg_shifts, ::GroundExcited) = nothing

"""
get_agg_hamiltonian(agg, aggTools.indices, franck_condon_factors; 
\t, groundEnergy = true)
get_agg_hamiltonian(agg, aggTools.indices; groundEnergy = true)
get_agg_hamiltonian(agg; groundEnergy = true)

Get Hamiltonian of the [`Aggregate`](@ref).

# Arguments
* `agg`: Instance of [`Aggregate`](@ref).
* `aggTools.indices`: Aggregate indices generated by [`get_indices`](@ref).
* `franck_condon_factors`: Franck-Condon factors generated by [`get_franck_condon_factors`](@ref).
"""
function get_agg_hamiltonian(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    groundEnergy::Bool = true,
    vib_basis::VibBasisLike=GroundExcited()
)
    Ham_I = get_agg_ham_interaction(aggCore, aggTools; vib_basis = vib_basis)
    Ham_0 = get_agg_ham_system_bath(aggCore, aggTools; groundEnergy = groundEnergy)
    return Ham_0 + Ham_I
end

### Sparse versions
#=
"""
get_agg_hamiltonianSparse(agg, aggTools.indices)
get_agg_hamiltonianSparse(agg)

Sparse version of [`get_agg_hamiltonian`](@ref).

"""
function get_agg_hamiltonianSparse(
agg::AggregateCore{T,C1,C2},
aggTools.indices::Any,
franck_condon_factors::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType}
if aggTools.indices === nothing
    aggTools.indices = get_indices(agg)
end
aggTools.bSize = length(aggTools.indices)
aggCore.molCount = length(agg.molecules)
if franck_condon_factors === nothing
    franck_condon_factors = get_franck_condon_factors(agg, aggTools.indices)
end
Ham = spzeros(C2, aggTools.bSize, aggTools.bSize)
for I = 1:aggTools.bSize
    elind1, vibind1 = aggTools.indices[I]
    elOrder1 = elIndOrder(elind1)
    for J = 1:aggTools.bSize
        elind2, vibind2 = aggTools.indices[J]
        elOrder2 = elIndOrder(elind2)
        if I == J
            val = get_agg_state_energy(agg, elind1, vibind1)
            if val != 0.0
                Ham[I, J] = val
            end
        else
            if elind1 != elind2
                val = agg.coupling[elOrder1, elOrder2] * franck_condon_factors[I, J]
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

get_agg_hamiltonianSparse(
agg::AggregateCore{T,C1,C2}
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = get_agg_hamiltonianSparse(
agg::AggregateCore{T,C1,C2},
nothing,
nothing
)

get_agg_hamiltonianSparse(
agg::AggregateCore{T,C1,C2},
aggTools.indices::Any
) where {T<:Integer,C1<:ComputableType,C2<:ComputableType} = get_agg_hamiltonianSparse(
agg::AggregateCore{T,C1,C2},
aggTools.indices,
nothing
)
=#

struct AggregateOperators{O_sys,O_bath,O,VB<:AbstractVibBasis} <: AbstractAggregateOperators
    Ham_sys::O_sys
    Ham_bath::O_bath
    Ham_S::O
    Ham_B::O
    Ham_0::O
    Ham_I::O
    Ham::O
    groundEnergy::Bool
    vib_basis::VB
    function AggregateOperators{O_sys,O_bath,O,VB}(
        Ham_sys::O_sys,
        Ham_bath::O_bath,
        Ham_S::O,
        Ham_B::O,
        Ham_0::O,
        Ham_I::O,
        Ham::O,
        groundEnergy::Bool,
        vib_basis::VB,
    ) where {O_sys<:Operator,O_bath<:Operator,O<:Operator,VB<:AbstractVibBasis}
        new(Ham_sys, Ham_bath, Ham_S, Ham_B, Ham_0, Ham_I, Ham, groundEnergy, vib_basis)
    end
end

function AggregateOperators(
    aggCore::AggregateCore,
    aggTools::AggregateTools;
    groundEnergy::Bool = true,
    vib_basis::VibBasisLike=GroundExcited()
)
    vb = _to_vib_basis(vib_basis)

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools; vib_basis = vb, groundEnergy = true)
    Ham_bath = get_agg_ham_bath_small(aggCore, aggTools; groundEnergy = true)
    Ham_S = get_agg_ham_system_big(aggCore, aggTools; vib_basis = vb, groundEnergy = groundEnergy)
    Ham_B = get_agg_ham_bath_big(aggCore, aggTools; groundEnergy = groundEnergy)
    Ham_0 = get_agg_ham_system_bath(aggCore, aggTools; groundEnergy = groundEnergy)
    Ham_I = get_agg_ham_interaction(aggCore, aggTools; vib_basis = vb)
    Ham = get_agg_hamiltonian(aggCore, aggTools; groundEnergy = groundEnergy, vib_basis = vb)

    O_sys = typeof(Ham_sys)
    O_bath = typeof(Ham_bath)
    O = typeof(Ham_S)
    VB = typeof(vb)

    return AggregateOperators{O_sys,O_bath,O,VB}(
        Ham_sys,
        Ham_bath,
        Ham_S,
        Ham_B,
        Ham_0,
        Ham_I,
        Ham,
        groundEnergy,
        vb
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
