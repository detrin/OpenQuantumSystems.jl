using Test
using OpenQuantumSystems

@testset "vib_basis types" begin
    @testset "type hierarchy" begin
        @test GroundGround() isa AbstractVibBasis
        @test GroundExcited() isa AbstractVibBasis
        @test GroundGround() == GroundGround()
        @test GroundExcited() == GroundExcited()
        @test GroundGround() != GroundExcited()
    end

    @testset "conversion helper" begin
        @test OpenQuantumSystems._to_vib_basis(:ground_ground) isa GroundGround
        @test OpenQuantumSystems._to_vib_basis(:ground_excited) isa GroundExcited
        @test OpenQuantumSystems._to_vib_basis(GroundGround()) isa GroundGround
        @test OpenQuantumSystems._to_vib_basis(GroundExcited()) isa GroundExcited
        @test_throws ArgumentError OpenQuantumSystems._to_vib_basis(:invalid)
    end

    @testset "symbol-to-type round-trip" begin
        @test OpenQuantumSystems._vib_basis_to_symbol(GroundGround()) === :ground_ground
        @test OpenQuantumSystems._vib_basis_to_symbol(GroundExcited()) === :ground_excited
    end

    @testset "backward compatibility with Symbol API" begin
        Mode_a = Mode(0.2, 1.0)
        Mode_b = Mode(0.3, 0.5)
        mol1 = Molecule([Mode_a], 3, [2.0, 200.0])
        mol2 = Molecule([Mode_b], 3, [3.0, 300.0])
        aggCore = AggregateCore([mol1, mol2])
        aggTools = AggregateTools(aggCore)

        ham_sym = get_agg_ham_system_small(aggCore, aggTools; vib_basis=:ground_ground)
        ham_type = get_agg_ham_system_small(aggCore, aggTools; vib_basis=GroundGround())
        @test ham_sym == ham_type

        ham_sym2 = get_agg_ham_system_small(aggCore, aggTools; vib_basis=:ground_excited)
        ham_type2 = get_agg_ham_system_small(aggCore, aggTools; vib_basis=GroundExcited())
        @test ham_sym2 == ham_type2

        ham_I_sym = get_agg_ham_interaction(aggCore, aggTools; vib_basis=:ground_ground)
        ham_I_type = get_agg_ham_interaction(aggCore, aggTools; vib_basis=GroundGround())
        @test ham_I_sym == ham_I_type

        ham_I_sym2 = get_agg_ham_interaction(aggCore, aggTools; vib_basis=:ground_excited)
        ham_I_type2 = get_agg_ham_interaction(aggCore, aggTools; vib_basis=GroundExcited())
        @test ham_I_sym2 == ham_I_type2
    end

    @testset "AggregateOperators stores type" begin
        Mode_a = Mode(0.2, 1.0)
        mol1 = Molecule([Mode_a], 3, [2.0, 200.0])
        mol2 = Molecule([Mode_a], 3, [3.0, 300.0])
        aggCore = AggregateCore([mol1, mol2])
        aggTools = AggregateTools(aggCore)

        ops_sym = AggregateOperators(aggCore, aggTools; vib_basis=:ground_ground)
        @test ops_sym.vib_basis isa GroundGround

        ops_type = AggregateOperators(aggCore, aggTools; vib_basis=GroundGround())
        @test ops_type.vib_basis isa GroundGround
        @test ops_sym == ops_type

        ops_ge = AggregateOperators(aggCore, aggTools; vib_basis=:ground_excited)
        @test ops_ge.vib_basis isa GroundExcited
    end

    @testset "setup_aggregate with type" begin
        Mode_a = Mode(0.2, 1.0)
        mol1 = Molecule([Mode_a], 3, [2.0, 200.0])
        mol2 = Molecule([Mode_a], 3, [3.0, 300.0])
        aggCore = AggregateCore([mol1, mol2])

        agg_sym = setup_aggregate(aggCore; vib_basis=:ground_ground)
        agg_type = setup_aggregate(aggCore; vib_basis=GroundGround())
        @test agg_sym == agg_type

        agg_ge = setup_aggregate(aggCore; vib_basis=:ground_excited)
        @test agg_ge.operators.vib_basis isa GroundExcited
    end

    @testset "trace_bath with type dispatch" begin
        Mode_a = Mode(0.2, 1.0)
        mol1 = Molecule([Mode_a], 3, [2.0, 200.0])
        mol2 = Molecule([Mode_a], 3, [3.0, 300.0])
        aggCore = AggregateCore([mol1, mol2])
        agg = setup_aggregate(aggCore; vib_basis=:ground_ground)

        W = DenseOperator(agg.tools.basis, agg.tools.basis, ones(agg.tools.bSize, agg.tools.bSize))

        rho_sym = trace_bath(W, aggCore, agg.operators, agg.tools; vib_basis=:ground_ground)
        rho_type = trace_bath(W, aggCore, agg.operators, agg.tools; vib_basis=GroundGround())
        @test rho_sym == rho_type
    end
end
