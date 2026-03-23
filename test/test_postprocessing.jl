using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "postprocessing" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))


    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 3, [2.0, 200.0])
    mol2 = Molecule([mode2], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setup_aggregate(aggCore)
    aggTools = agg.tools
    aggOperators = agg.operators

    Ham_B = agg.operators.Ham_B
    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct

    # recast of OperatorVectorArray
    tspan = get_tspan(0., 0.02, 100)
    W0, rho0, W0_bath = ultrafast_laser_excitation(10., [0.0, 0.3, 0.7], agg)
    _, rho_t_1 = Evolution_sS_exact(W0, tspan, agg)
    _, rho_t_2 = LvN_sS(
        W0,
        tspan,
        agg;
        reltol = 1e-12,
        abstol = 1e-12,
        alg = OrdinaryDiffEq.Tsit5(),
    )

    rho_t_1_recasted = operator_recast(rho_t_1)
    @test typeof(rho_t_1) == OperatorVectorArray
    @test typeof(rho_t_1_recasted) == OperatorVectorArray

    rho_t_2_recasted = operator_recast(rho_t_2)
    @test typeof(rho_t_2) == OperatorVector
    @test typeof(rho_t_2_recasted) == OperatorVectorArray

    # interaction and schrodinger picture
    _, rho_int_t_ref = Evolution_sI_exact(W0, tspan, agg)
    rho_t = interaction_pic_to_schroedinger_pic(rho_int_t_ref, tspan, agg)
    rho_int_t = schroedinger_pic_to_interaction_pic(rho_t, tspan, agg)
    @test D(rho_int_t_ref, rho_int_t) < 1e-14

    _, rho_int_t_ref = LvN_sI(
        W0,
        tspan,
        agg;
        reltol = 1e-12,
        abstol = 1e-12,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    rho_t = interaction_pic_to_schroedinger_pic(rho_int_t_ref, tspan, agg)
    rho_int_t = schroedinger_pic_to_interaction_pic(rho_t, tspan, agg)
    rho_int_t_ref = operator_recast(rho_int_t_ref)
    @test D(rho_int_t_ref, rho_int_t) < 1e-14

    # local and exciton basis
    _, rho_t_ref = Evolution_sS_exact(W0, tspan, agg)
    rho_t_exc = local_st_to_exciton_st(rho_t_ref, agg)
    rho_t = exciton_st_to_local_st(rho_t_exc, agg)
    @test D(rho_t_ref, rho_t) < 1e-14

    _, rho_t_ref = LvN_sS(
        W0,
        tspan,
        agg;
        reltol = 1e-12,
        abstol = 1e-12,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    rho_t_exc = local_st_to_exciton_st(rho_t_ref, agg)
    rho_t = exciton_st_to_local_st(rho_t_exc, agg)
    rho_t_ref = operator_recast(rho_t_ref)
    @test D(rho_t_ref, rho_t) < 1e-14

    # validate_state
    b = GenericBasis([2])
    valid_rho = Operator(b, b, ComplexF64[0.5 0.0; 0.0 0.5])
    @test validate_state(valid_rho) == true

    bad_trace = Operator(b, b, ComplexF64[0.5 0.0; 0.0 0.0])
    @test_logs (:warn, r"trace deviates") validate_state(bad_trace) == false
    @test validate_state(bad_trace) == false

    nan_rho = Operator(b, b, ComplexF64[NaN 0.0; 0.0 0.5])
    @test_logs (:warn, r"Inf or NaN") validate_state(nan_rho) == false
    @test validate_state(nan_rho) == false

    inf_rho = Operator(b, b, ComplexF64[Inf 0.0; 0.0 0.5])
    @test_logs (:warn, r"Inf or NaN") (:warn, r"trace deviates") validate_state(inf_rho) == false
    @test validate_state(inf_rho) == false

    # validate_trajectory with OperatorVector
    valid_vec = [valid_rho, valid_rho]
    @test validate_trajectory(valid_vec) == true

    bad_vec = [valid_rho, nan_rho]
    @test validate_trajectory(bad_vec) == false

    # validate_trajectory with OperatorVectorArray
    valid_array = zeros(ComplexF64, 2, 2, 2)
    valid_array[1, :, :] = valid_rho.data
    valid_array[2, :, :] = valid_rho.data
    @test validate_trajectory(valid_array) == true

    bad_array = zeros(ComplexF64, 2, 2, 2)
    bad_array[1, :, :] = valid_rho.data
    bad_array[2, :, :] = nan_rho.data
    @test validate_trajectory(bad_array) == false
end