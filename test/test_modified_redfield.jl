using Test
using OpenQuantumSystems
using LinearAlgebra

@testset "modified_redfield" begin

    @testset "exciton_basis — dimer" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        aggTools = AggregateTools(aggCore)
        energies, coefficients = exciton_basis(aggCore, aggTools)

        @test length(energies) == 2
        @test energies[1] < energies[2]

        # Coefficients are orthonormal
        @test coefficients * coefficients' ≈ I atol=1e-10

        # Eigenvalues: account for reorganization energy shift (ω × shift²/2 per mode)
        reorg = 200.0 * (2 * 0.02)^2 / 2.0
        H_exc = [12500.0+reorg J; J 12700.0+reorg]
        vals_ref = eigvals(Hermitian(H_exc))
        @test energies ≈ vals_ref atol=1e-10
    end

    @testset "exciton_basis — no coupling" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])

        aggTools = AggregateTools(aggCore)
        energies, coefficients = exciton_basis(aggCore, aggTools)

        reorg = 200.0 * (2 * 0.02)^2 / 2.0
        @test energies ≈ [12500.0 + reorg, 12700.0 + reorg] atol=1e-10
        @test abs(coefficients[1, 1]) ≈ 1.0 atol=1e-10
        @test abs(coefficients[2, 2]) ≈ 1.0 atol=1e-10
    end

    @testset "modified_redfield_rates — basic properties" begin
        # Use stronger system-bath coupling for well-converged rates
        mode = Mode(; omega = 200.0, hr_factor = 0.3)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        rates, energies, coefficients = modified_redfield_rates(
            aggCore; T = 300.0, t_max = 200.0
        )

        @test size(rates) == (2, 2)

        # Columns sum to zero (population conservation)
        for j in 1:2
            @test sum(rates[:, j]) ≈ 0.0 atol=1e-8
        end

        # Off-diagonal rates should be positive and finite
        @test rates[1, 2] > 0.0
        @test rates[2, 1] > 0.0
        @test isfinite(rates[1, 2])
        @test isfinite(rates[2, 1])

        # Downhill rate should be larger at finite T
        @test rates[1, 2] > rates[2, 1]
    end

    @testset "modified_redfield_rates — detailed balance trend" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.3)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        T = 300.0
        rates, energies, _ = modified_redfield_rates(aggCore; T = T, t_max = 200.0)

        # Downhill rate (2→1) should be larger than uphill (1→2)
        @test rates[1, 2] > rates[2, 1]

        # Ratio should be < 1 (uphill/downhill < 1)
        ratio = rates[2, 1] / rates[1, 2]
        @test 0.0 < ratio < 1.0
    end

    @testset "modified_redfield_dynamics — population conservation" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.3)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        tspan = collect(0.0:0.01:1.0)
        rho0 = [1.0, 0.0]

        ts, pop = modified_redfield_dynamics(aggCore, rho0, tspan; T = 300.0, t_max = 200.0)

        @test size(pop) == (length(tspan), 2)

        for i in 1:length(tspan)
            @test sum(pop[i, :]) ≈ 1.0 atol=1e-10
        end

        @test pop[1, :] ≈ rho0
    end

    @testset "modified_redfield_dynamics — transfer occurs" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.3)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        tspan = collect(0.0:0.01:5.0)
        rho0 = [0.0, 1.0]

        ts, pop = modified_redfield_dynamics(aggCore, rho0, tspan; T = 300.0, t_max = 200.0)

        # Starting in exciton 2, population should transfer to exciton 1 (downhill)
        @test pop[end, 1] > pop[1, 1]
        @test pop[end, 2] < pop[1, 2]
    end

    @testset "modified_redfield_rates — zero coupling" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])

        rates, _, _ = modified_redfield_rates(aggCore; T = 300.0, t_max = 500.0)

        # Without coupling, excitons are localized — no transfer
        @test abs(rates[1, 2]) < 1e-10
        @test abs(rates[2, 1]) < 1e-10
    end

    @testset "modified_redfield_rates — gamma damping" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        # With damping, integral converges even for weak system-bath coupling
        rates, _, _ = modified_redfield_rates(aggCore; T = 300.0, t_max = 500.0, gamma = 10.0)

        @test size(rates) == (2, 2)
        for j in 1:2
            @test sum(rates[:, j]) ≈ 0.0 atol=1e-8
        end
        @test isfinite(rates[1, 2])
        @test isfinite(rates[2, 1])
    end

end
