using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

import OrdinaryDiffEq

@testset "liouville" begin

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
    aggOperators = agg.operators
    aggTools = agg.tools

    Ham = agg.operators.Ham
    Ham_0 = agg.operators.Ham_0
    basis = agg.tools.basis

    ket0 = randstate(basis)
    W0 = dm(ket0)
    # tests have to be quick enough
    tspan = [0.0:0.1:1.0;]

    _, rho_int_t = LvN_sI(
        W0,
        tspan,
        agg;
        reltol = 1e-10,
        abstol = 1e-10,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    for t_i = 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        U_0_op = evolutionOperator(Ham_0, t)
        W_int = U_0_op' * W * U_0_op
        rho_int = trace_bath(W_int, aggCore, aggOperators, aggTools)
        @test 1e-7 > D(rho_int, rho_int_t[t_i])
        # println(t_i, " ", D(W_int, W_int_t[t_i]))
    end

    _, rho_t = LvN_sS(
        W0,
        tspan,
        agg;
        reltol = 1e-10,
        abstol = 1e-10,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    for t_i = 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        rho = trace_bath(W, aggCore, aggOperators, aggTools)
        @test 1e-7 > D(rho, rho_t[t_i])
        # println(t_i, " ", D(W_int, W_int_t[t_i]))
    end

    _, W_int_t = LvN_SI(
        W0,
        tspan,
        agg;
        reltol = 1e-10,
        abstol = 1e-10,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    for t_i = 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        U_0_op = evolutionOperator(Ham_0, t)
        W_int = U_0_op' * W * U_0_op
        @test 1e-7 > D(W_int, W_int_t[t_i])
        # println(t_i, " ", D(W_int, W_int_t[t_i]))
    end

    _, W_t = LvN_SS(
        W0,
        tspan,
        agg;
        reltol = 1e-10,
        abstol = 1e-10,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    for t_i = 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        @test 1e-7 > D(W, W_t[t_i])
        # println(t_i, " ", D(rho.data, rho_t[t_i].data))
    end

end # testset
