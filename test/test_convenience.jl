using Test
using OpenQuantumSystems

@testset "convenience" begin

    @testset "setup_dimer" begin
        agg = setup_dimer()
        @test agg isa Aggregate
        @test agg.core.molCount == 2
        @test agg.core.coupling[2, 3] == 100.0
        @test agg.core.coupling[3, 2] == 100.0
        @test agg.core.molecules[1].E == [0.0, 12500.0]
        @test agg.core.molecules[2].E == [0.0, 12700.0]
        @test agg.core.molecules[1].Nvib == 3
        @test length(agg.core.molecules[1].modes) == 1
        @test agg.core.molecules[1].modes[1].omega == 200.0

        agg2 = setup_dimer(E1 = 10000.0, E2 = 11000.0, J = 50.0, mode_omega = 150.0, mode_hr = 0.05, nvib = 4)
        @test agg2.core.molecules[1].E == [0.0, 10000.0]
        @test agg2.core.molecules[2].E == [0.0, 11000.0]
        @test agg2.core.coupling[2, 3] == 50.0
        @test agg2.core.molecules[1].Nvib == 4
        @test agg2.core.molecules[1].modes[1].omega == 150.0
    end

    @testset "setup_trimer" begin
        J = zeros(Float64, 3, 3)
        J[1, 2] = 100.0; J[2, 1] = 100.0
        J[2, 3] = 80.0;  J[3, 2] = 80.0
        agg = setup_trimer(J_matrix = J)
        @test agg isa Aggregate
        @test agg.core.molCount == 3
        @test agg.core.coupling[2, 3] == 100.0
        @test agg.core.coupling[3, 4] == 80.0
        @test agg.core.coupling[4, 3] == 80.0
        @test agg.core.molecules[1].E == [0.0, 12500.0]
        @test agg.core.molecules[3].E == [0.0, 12900.0]

        @test_throws ArgumentError setup_trimer(E = [1.0, 2.0])
        @test_throws ArgumentError setup_trimer(J_matrix = zeros(Float64, 2, 2))
    end

    @testset "setup_linear_chain" begin
        agg = setup_linear_chain()
        @test agg isa Aggregate
        @test agg.core.molCount == 2
        @test agg.core.coupling[2, 3] == 100.0

        energies = [12000.0, 12200.0, 12400.0, 12600.0]
        agg4 = setup_linear_chain(energies = energies, J_nearest = 75.0, nvib = 2)
        @test agg4.core.molCount == 4
        @test agg4.core.coupling[2, 3] == 75.0
        @test agg4.core.coupling[3, 4] == 75.0
        @test agg4.core.coupling[4, 5] == 75.0
        @test agg4.core.coupling[2, 4] == 0.0
        @test agg4.core.molecules[4].E == [0.0, 12600.0]
        @test agg4.core.molecules[1].Nvib == 2

        @test_throws ArgumentError setup_linear_chain(energies = [12000.0])
    end

end
