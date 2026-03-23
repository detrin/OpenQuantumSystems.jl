# Backward-compatible aliases for renamed functions.
# Old camelCase names map to the new snake_case canonical names.
# These aliases are exported so that existing user code continues to work.

# molecules.jl
const franckCondonFactors = franck_condon_factors
const vibrationalIndices = vibrational_indices
const electronicIndices = electronic_indices
const getMolStateEnergy = get_mol_state_energy
const getMolStateFC = get_mol_state_fc
const getMolShifts = get_mol_shifts
const getMolFrequencies = get_mol_frequencies

# aggregateCore.jl
const getNvib = get_nvib
const getShifts = get_shifts
const getFrequencies = get_frequencies
const getAggStateEnergy = get_agg_state_energy

# aggregateTools.jl
const getIndices = get_indices
const getIndicesMap = get_indices_map
const getFranckCondonFactors = get_franck_condon_factors
const getFCproduct = get_fc_product

# aggregateOperators.jl
const getAggHamSystemSmall = get_agg_ham_system_small
const getAggHamSystemBig = get_agg_ham_system_big
const getAggHamBathSmall = get_agg_ham_bath_small
const getAggHamBathBig = get_agg_ham_bath_big
const getAggHamSystemBath = get_agg_ham_system_bath
const getAggHamInteraction = get_agg_ham_interaction
const getAggHamiltonian = get_agg_hamiltonian

# aggregate.jl
const setupAggregate = setup_aggregate
const setupAggregate! = setup_aggregate!

# evolution.jl
const evolutionOperator = evolution_operator
const evolutionOperatorA = evolution_operator_a
const evolutionSuperOperator = evolution_super_operator
const evolutionOperatorArray = evolution_operator_array
const evolutionSuperOperatorArray = evolution_super_operator_array
const evolutionOperatorIterator = evolution_operator_iterator
const evolutionSuperOperatorIterator = evolution_super_operator_iterator

# interaction_picture.jl
const getInteractionHamIPicture = get_interaction_ham_i_picture
const getInteractionHamIPictureA = get_interaction_ham_i_picture_a
