using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "scoring" begin

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
    agg = setupAggregate(aggCore)
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

    tspan = get_tspan(0., 0.02, 100)
    W0, rho0, W0_bath = ultrafast_laser_excitation(10., [0.0, 0.3, 0.7], agg)
    _, W_t_1 = Evolution_SS_exact(W0, tspan, agg)
    _, W_t_2 = Evolution_SS_exact(W0, tspan, agg)
    rmse = get_rmse_in_time(W_t_1, W_t_2)
    @test rmse == 0.

    _, W_t_3 = LvN_SS(
        W0,
        tspan,
        agg;
        reltol = 1e-10,
        abstol = 1e-10,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    rmse = get_rmse_in_time(W_t_3, W_t_3)
    @test rmse == 0.
    rmse = get_rmse_in_time(W_t_1, W_t_3)
    @test rmse < 1e-9
    rmse = get_rmse_in_time(W_t_3, W_t_1)
    @test rmse < 1e-9


    tspan = get_tspan(0., 0.02, 100)
    W0, rho0, W0_bath = ultrafast_laser_excitation(10., [0.0, 0.3, 0.7], agg)
    _, rho_t_1 = Evolution_sS_exact(W0, tspan, agg)
    _, rho_t_2 = Evolution_sS_exact(W0, tspan, agg)
    score = compare_rho(rho_t_1, rho_t_2)
    @test sum(score) == 0.

    _, rho_t_3 = LvN_sS(
        W0,
        tspan,
        agg;
        reltol = 1e-12,
        abstol = 1e-12,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    score = compare_rho(rho_t_3, rho_t_3)
    @test sum(score) == 0.

    score = compare_rho(rho_t_1, rho_t_3)
    @test sum(score) < 1e-9

    score = compare_rho(rho_t_3, rho_t_1)
    @test sum(score) < 1e-9

    tspan = get_tspan(0., 0.02, 100)
    _, rho_t_1 = Evolution_sS_exact(W0, tspan, agg)
    _, rho_t_2 = LvN_sS(
        W0,
        tspan,
        agg;
        reltol = 1e-12,
        abstol = 1e-12,
        alg = OrdinaryDiffEq.Tsit5(),
    )

    score_ref = compare_rho(rho_t_2, rho_t_2)
    score_t = compare_rho_in_time(rho_t_2, rho_t_2; smooth_const=1e-9)
    N, M, K = size(rho_t_1)
    score = zeros(Float64, M, K)
    for t_i in 1:N
        score[:, :] += score_t[t_i, :, :]
    end
    @test D(score_ref, score) == 0.

    score_ref = compare_rho(rho_t_1, rho_t_2)
    score_t = compare_rho_in_time(rho_t_1, rho_t_2; smooth_const=1e-9)
    N, M, K = size(rho_t_1)
    score = zeros(Float64, M, K)
    for t_i in 1:N
        score[:, :] += score_t[t_i, :, :]
    end
    @test D(score_ref, score) < 1e-8

    score_ref = compare_rho(rho_t_2, rho_t_1)
    score_t = compare_rho_in_time(rho_t_2, rho_t_1; smooth_const=1e-9)
    N, M, K = size(rho_t_1)
    score = zeros(Float64, M, K)
    for t_i in 1:N
        score[:, :] += score_t[t_i, :, :]
    end
    @test D(score_ref, score) < 1e-8

    score_ref = compare_rho(rho_t_1, rho_t_1)
    score_t = compare_rho_in_time(rho_t_1, rho_t_1; smooth_const=1e-9)
    N, M, K = size(rho_t_1)
    score = zeros(Float64, M, K)
    for t_i in 1:N
        score[:, :] += score_t[t_i, :, :]
    end
    @test D(score_ref, score) == 0.
end