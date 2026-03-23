using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "simulation_result" begin

    Random.seed!(StableRNG(0), 1)

    b = GenericBasis([2])
    op1 = DenseOperator(b, b, ComplexF64[0.6 0.1; 0.1 0.4])
    op2 = DenseOperator(b, b, ComplexF64[0.5 0.2; 0.2 0.5])
    op3 = DenseOperator(b, b, ComplexF64[0.4 0.05; 0.05 0.6])
    states = [op1, op2, op3]
    tspan = [0.0, 0.5, 1.0]

    @testset "constructor and fields" begin
        r = SimulationResult(tspan, states, :lvn_si)
        @test r.tspan == tspan
        @test r.states == states
        @test r.method == :lvn_si
        @test isempty(r.extra)
    end

    @testset "constructor with extra" begin
        extra = Dict{Symbol,Any}(:W_1_bath_t => [1, 2, 3])
        r = SimulationResult(tspan, states, :qme_iterative, extra)
        @test r.extra[:W_1_bath_t] == [1, 2, 3]
    end

    @testset "length" begin
        r = SimulationResult(tspan, states, :lvn_si)
        @test length(r) == 3
    end

    @testset "getindex" begin
        r = SimulationResult(tspan, states, :lvn_si)
        @test r[1] === op1
        @test r[2] === op2
        @test r[3] === op3
    end

    @testset "firstindex and lastindex" begin
        r = SimulationResult(tspan, states, :lvn_si)
        @test firstindex(r) == 1
        @test lastindex(r) == 3
        @test r[end] === op3
    end

    @testset "iterate" begin
        r = SimulationResult(tspan, states, :lvn_si)
        collected = collect(r)
        @test collected == states
    end

    @testset "populations" begin
        r = SimulationResult(tspan, states, :lvn_si)
        pops = populations(r)
        @test size(pops) == (3, 2)
        @test pops[1, 1] ≈ 0.6
        @test pops[1, 2] ≈ 0.4
        @test pops[2, 1] ≈ 0.5
        @test pops[2, 2] ≈ 0.5
        @test pops[3, 1] ≈ 0.4
        @test pops[3, 2] ≈ 0.6
    end

    @testset "populations returns real values" begin
        r = SimulationResult(tspan, states, :lvn_si)
        pops = populations(r)
        @test eltype(pops) == Float64
    end

end
