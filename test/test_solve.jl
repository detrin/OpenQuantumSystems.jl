using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

import OrdinaryDiffEq

@testset "solve" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    mol1 = Molecule([mode1], 3, [2.0, 200.0])
    mol2 = Molecule([mode2], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)

    Ham = agg.operators.Ham
    basis = agg.tools.basis

    ket0 = randstate(basis)
    W0 = dm(ket0)
    tspan = [0.0:0.1:1.0;]

    @testset "lvn_si dispatches to LvN_sI" begin
        result = solve(W0, tspan, agg; method=:lvn_si,
            reltol=1e-10, abstol=1e-10, alg=OrdinaryDiffEq.Tsit5())
        t2, rho2 = LvN_sI(W0, tspan, agg;
            reltol=1e-10, abstol=1e-10, alg=OrdinaryDiffEq.Tsit5())
        @test result isa SimulationResult
        @test result.method == :lvn_si
        @test length(result) == length(rho2)
        for i in 1:length(result)
            @test D(result[i], rho2[i]) < 1e-12
        end
    end

    @testset "lvn_SS dispatches to LvN_SS" begin
        result = solve(W0, tspan, agg; method=:lvn_SS,
            reltol=1e-10, abstol=1e-10, alg=OrdinaryDiffEq.Tsit5())
        t2, W2 = LvN_SS(W0, tspan, agg;
            reltol=1e-10, abstol=1e-10, alg=OrdinaryDiffEq.Tsit5())
        @test result isa SimulationResult
        @test result.method == :lvn_SS
        @test length(result) == length(W2)
        for i in 1:length(result)
            @test D(result[i], W2[i]) < 1e-12
        end
    end

    @testset "unknown method throws" begin
        @test_throws ArgumentError solve(W0, tspan, agg; method=:nonexistent)
    end

    @testset "qme_iterative without args throws" begin
        @test_throws ArgumentError solve(W0, tspan, agg; method=:qme_iterative)
    end

end
