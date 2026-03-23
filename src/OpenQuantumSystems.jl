module OpenQuantumSystems

using Reexport
@reexport using QuantumOpticsBase:
    Basis,
    GenericBasis,
    basis,
    tensor,
    ⊗,
    permutesystems,
    @samebases,

    # states.jl
    StateVector,
    Bra,
    Ket,
    basisstate,
    norm,
    dagger,
    normalize,
    normalize!,

    # operators.jl
    AbstractOperator,
    DataOperator,
    expect,
    variance,
    identityoperator,
    ptrace,
    embed,
    dense,
    tr,
    sparse,

    # operators_dense.jl
    Operator,
    DenseOperator,
    DenseOpType,
    projector,
    dm,

    # operators_sparse.jl
    SparseOperator,
    diagonaloperator,
    SparseOpType,

    # superoperators.jl
    AbstractSuperOperator,
    SuperOperator,
    DenseSuperOperator,
    DenseSuperOpType,
    SparseSuperOperator,
    SparseSuperOpType,
    spre,
    spost,

    # fock.jl
    FockBasis,
    number,
    destroy,
    create,
    fockstate,
    coherentstate,
    coherentstate!,
    displace,

    # spin.jl
    SpinBasis,
    sigmax,
    sigmay,
    sigmaz,
    sigmap,
    sigmam,
    spinup,
    spindown,

    # state_definitions.jl
    randstate,
    randoperator,
    thermalstate,
    coherentthermalstate,
    phase_average,
    passive_state,

    # nlevel.jl
    NLevelBasis,
    transition,
    nlevelstate,

    # metrics.jl
    tracenorm,
    tracenorm_h,
    tracenorm_nh,
    tracedistance,
    tracedistance_h,
    tracedistance_nh,
    entropy_vn,
    entropy_renyi,
    fidelity,
    ptranspose,
    PPT,
    negativity,
    logarithmic_negativity,
    entanglement_entropy,

    # pauli.jl
    PauliBasis,
    PauliTransferMatrix,
    DensePauliTransferMatrix,
    ChiMatrix,
    DenseChiMatrix,
    avg_gate_fidelity

using SparseArrays, LinearAlgebra

export
    # core.jl
    AbstractVibBasis,
    GroundGround,
    GroundExcited,

    # operators_dense.jl
    AnnihilationOperator,
    CreationOperator,
    PositionOperator,
    MomentumOperator,
    ShiftOperator,
    OneDenseOperator,

    # superoperators.jl
    Commutator,

    # metrics.jl
    tracenorm,
    tracenorm_h,
    tracenorm_nh,
    tracedistance,
    tracedistance_h,
    tracedistance_nh,

    # molecules.jl
    Mode,
    convert_units,
    franck_condon_factors,
    vibrational_indices,
    electronic_indices,
    Molecule,
    get_mol_state_energy,
    get_mol_state_fc,
    get_mol_shifts,
    get_mol_frequencies,

    # aggregateCore.jl
    AggregateCore,
    get_nvib,
    get_shifts,
    get_frequencies,
    get_agg_state_energy,

    # aggregateTools.jl
    AggregateTools,
    get_indices,
    get_indices_map,
    get_franck_condon_factors,
    get_fc_product,
    take_el_part,

    # aggregateOperators.jl
    get_agg_ham_system_small,
    get_agg_ham_system_big,
    get_agg_ham_bath_small,
    get_agg_ham_bath_big,
    get_agg_ham_system_bath,
    get_agg_ham_interaction,
    get_agg_hamiltonian,
    AggregateOperators,

    # aggregate.jl
    Aggregate,
    setup_aggregate,
    setup_aggregate!,

    # evolution.jl
    get_tspan,
    evolution_operator,
    evolution_operator_a,
    evolution_super_operator,
    evolution_operator_array,
    evolution_super_operator_array,
    evolution_operator_iterator,
    evolution_super_operator_iterator,
    evolutionExact,
    evolutionExact!,
    evolutionApproximate,
    evolutionApproximate!,
    evolution_exact,
    evolution_approximate,
    evolution_el_part,
    Evolution_SI_exact,
    Evolution_sI_exact,
    Evolution_SS_exact,
    Evolution_sS_exact,

    # schrodinger.jl
    schroedinger,
    schroedinger_dynamic,

    # liouville.jl
    LvN_sI,
    LvN_sS,
    LvN_SI,
    LvN_SS,

    # interaction_picture.jl
    get_interaction_ham_i_picture,
    get_interaction_ham_i_picture_a,

    # master_exact.jl
    QME_SI_exact,
    QME_SS_exact,
    QME_sI_exact,
    QME_sS_exact,

    # trace.jl
    trace_bath,
    trace_bath_slow,
    get_rho_bath,
    trace_bath_part,
    ad,
    correlation_function,

    # initial_state.jl
    thermal_state,
    thermal_state_old,
    thermal_state_composite,
    thermal_state_composite_old,
    ultrafast_laser_excitation,

    # memory_kernel.jl
    MemoryKernel_1_traced,
    MemoryKernel_2_traced,
    MemoryKernel_3_traced,
    MemoryKernel_4_traced,
    MemoryKernel_traced,

    # rate_constant.jl
    M_aabb_W_bath_intp,
    K_aabb_W_bath_intp,
    M_abcd_W_bath_intp,
    K_abcd_W_bath_intp,

    # master_ansatz.jl
    QME_sI_ansatz,
    QME_sI_ansatz_test,
    QME_sI_ansatz_const_int,
    QME_sI_ansatz_const_sch,
    QME_sI_ansatz_linear_sch,
    QME_sI_ansatz_linear2_sch,
    QME_sI_ansatz_upart1_sch,
    QME_sI_ansatz_upart1_int,
    QME_sI_ansatz_upart2_sch,
    QME_sI_ansatz_upart2_int,

    # redfield.jl
    QME_sI_Redfield,

    # master_iterative
    normalize_bath,
    W_1_bath,
    QME_sI_iterative,
    W_1_markov0_bath,
    QME_sI_iterative_markov0,
    W_1_markov1_bath,
    QME_sI_iterative_markov1,

    # postprocessing.jl
    OperatorVector,
    OperatorVectorArray,
    operator_recast,
    interaction_pic_to_schroedinger_pic,
    schroedinger_pic_to_interaction_pic,
    local_st_to_exciton_st,
    exciton_st_to_local_st,
    tspan_cm_to_fs,
    validate_state,
    validate_trajectory,

    # scoring.jl
    get_rmse_in_time,
    compare_rho,
    compare_rho_in_time,

    # simulation_result.jl
    SimulationResult,
    populations,

    # solve.jl
    solve,

    # forster.jl
    BOLTZMANN_CM,
    LineshapeResult,
    absorption_spectrum,
    emission_spectrum,
    spectral_overlap,
    forster_rate,
    forster_rate_matrix,

    # spectral_density.jl
    SpectralDensity,
    spectral_density,
    reorganization_energy,
    lineshape_function,
    lineshape_derivative,
    lineshape_second_derivative,

    # modified_redfield.jl
    exciton_basis,
    modified_redfield_rates,
    modified_redfield_dynamics,

    # dipole.jl
    DIPOLE_COUPLING_FACTOR,
    TransitionDipole,
    dipole_dipole_coupling,
    coupling_from_dipoles,

    # convenience.jl
    setup_dimer,
    setup_trimer,
    setup_linear_chain,

    # deprecated_aliases.jl (backward-compatible camelCase names)
    franckCondonFactors,
    vibrationalIndices,
    electronicIndices,
    getMolStateEnergy,
    getMolStateFC,
    getMolShifts,
    getMolFrequencies,
    getNvib,
    getShifts,
    getFrequencies,
    getAggStateEnergy,
    getIndices,
    getIndicesMap,
    getFranckCondonFactors,
    getFCproduct,
    getAggHamSystemSmall,
    getAggHamSystemBig,
    getAggHamBathSmall,
    getAggHamBathBig,
    getAggHamSystemBath,
    getAggHamInteraction,
    getAggHamiltonian,
    setupAggregate,
    setupAggregate!,
    evolutionOperator,
    evolutionOperatorA,
    evolutionSuperOperator,
    evolutionOperatorArray,
    evolutionSuperOperatorArray,
    evolutionOperatorIterator,
    evolutionSuperOperatorIterator,
    getInteractionHamIPicture,
    getInteractionHamIPictureA

# base/
include("base/core.jl")
include("base/operators_dense.jl")
include("base/superoperators.jl")
include("base/metrics.jl")

# aggregate/
include("aggregate/molecules.jl")
include("aggregate/aggregateCore.jl")
include("aggregate/aggregateTools.jl")
include("aggregate/aggregateOperators.jl")
include("aggregate/aggregate.jl")
include("aggregate/dipole.jl")
include("aggregate/convenience.jl")

# base/ (depends on aggregate types)
include("base/timeevolution_base.jl")

# evolution/
include("evolution/evolution.jl")
include("evolution/schroedinger.jl")
include("evolution/liouville.jl")
include("evolution/interaction_picture.jl")
include("evolution/master_exact.jl")
include("evolution/trace.jl")
include("evolution/memory_kernel.jl")
include("evolution/rate_constant.jl")
include("evolution/master_ansatz.jl")
include("evolution/redfield.jl")
include("evolution/master_iterative.jl")

# spectroscopy/
include("spectroscopy/forster.jl")
include("spectroscopy/spectral_density.jl")
include("spectroscopy/modified_redfield.jl")

# utils/
include("utils/initial_state.jl")
include("utils/postprocessing.jl")
include("utils/scoring.jl")
include("utils/simulation_result.jl")
include("utils/solve.jl")
include("utils/deprecated_aliases.jl")

end

