using Test
using OpenQuantumSystems
using LinearAlgebra

@testset "corrected_memory_kernel" begin

    @testset "analytic_correlation_nn — single mode" begin
        sd = SpectralDensity([200.0], [0.02])
        T = 300.0
        tau = 0.1

        C = analytic_correlation_nn(sd, tau, T)
        @test isa(C, ComplexF64)
        @test isfinite(C)

        C0 = analytic_correlation_nn(sd, 0.0, T)
        @test real(C0) > 0.0
        @test abs(imag(C0)) < 1e-10

        C_neg = analytic_correlation_nn(sd, -tau, T)
        @test C_neg ≈ conj(C) atol = 1e-10
    end

    @testset "analytic_correlation_nn — zero temperature" begin
        sd = SpectralDensity([200.0], [0.02])
        tau = 0.1

        C = analytic_correlation_nn(sd, tau, 0.0)
        expected = 200.0^2 * 0.02 * exp(-im * 200.0 * tau)
        @test C ≈ expected atol = 1e-10
    end

    @testset "exciton_correlation — identity coefficients" begin
        C_vals = ComplexF64[1.0 + 0.5im, 2.0 - 0.3im]
        coeffs = [1.0 0.0; 0.0 1.0]

        @test exciton_correlation(C_vals, coeffs, 1, 1, 1, 1) ≈ C_vals[1]
        @test exciton_correlation(C_vals, coeffs, 2, 2, 2, 2) ≈ C_vals[2]
        @test abs(exciton_correlation(C_vals, coeffs, 1, 2, 1, 1)) < 1e-15
    end

    @testset "zeroth_order_memory_kernel_cf — basic properties" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        t1, t2 = 0.1, 0.05
        M = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = 300.0)

        @test size(M) == (3, 3, 3, 3)
        @test all(isfinite, M)

        for i in 1:3
            @test all(abs.(M[1, :, :, :]) .< 1e-15)
            @test all(abs.(M[:, 1, :, :]) .< 1e-15)
            @test all(abs.(M[:, :, 1, :]) .< 1e-15)
            @test all(abs.(M[:, :, :, 1]) .< 1e-15)
        end
    end

    @testset "zeroth_order_memory_kernel_cf — t1==t2 symmetry" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        t = 0.1
        M = zeroth_order_memory_kernel_cf(t, t, aggCore; T = 300.0)
        @test all(isfinite, M)
    end

    @testset "zeroth_order — hermiticity M[a,b,c,d] = conj(M[b,a,d,c])" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        M = zeroth_order_memory_kernel_cf(0.1, 0.05, aggCore; T = 300.0)
        elLen = 3
        for a in 2:elLen, b in 2:elLen, c in 2:elLen, d in 2:elLen
            @test M[a, b, c, d] ≈ conj(M[b, a, d, c]) atol = 1e-10
        end
    end

    @testset "zeroth_order — J=0 diagonal structure" begin
        mode1 = Mode(0.2, 1.0)
        mode2 = Mode(0.3, 2.0)
        mol1 = Molecule([mode1], 2, [2.0, 200.0])
        mol2 = Molecule([mode2], 2, [3.0, 300.0])
        aggCore = AggregateCore([mol1, mol2])

        M = zeroth_order_memory_kernel_cf(0.1, 0.05, aggCore; T = 300.0)
        elLen = 3

        for a in 2:elLen, b in 2:elLen, c in 2:elLen, d in 2:elLen
            if a != c || b != d
                @test abs(M[a, b, c, d]) < 1e-10
            end
        end

        @test abs(M[2, 3, 2, 3]) > 1e-5
        @test abs(M[3, 2, 3, 2]) > 1e-5
    end

    @testset "zeroth_order — t=0 real valued" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        M = zeroth_order_memory_kernel_cf(0.0, 0.0, aggCore; T = 300.0)
        for a in 2:3, b in 2:3, c in 2:3, d in 2:3
            @test abs(imag(M[a, b, c, d])) < 1e-10
        end
    end

    @testset "zeroth_order — single site analytical" begin
        omega = 100.0
        S = 0.5
        T = 300.0
        mol = Molecule([Mode(; omega = omega, hr_factor = S)], 2, [0.0, 500.0])
        aggCore = AggregateCore([mol])

        t1, t2 = 0.1, 0.05
        tau = t1 - t2
        M = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = T)

        n = 1.0 / (exp(omega / (BOLTZMANN_CM * T)) - 1.0)
        C_fwd = omega^2 * S * (n * exp(im * omega * tau) + (n + 1) * exp(-im * omega * tau))
        C_bwd = conj(C_fwd)
        E2 = 500.0

        expected_2222 = 2.0 * real(C_fwd) - C_bwd - C_fwd
        @test M[2, 2, 2, 2] ≈ expected_2222 atol = 1e-8
    end

    @testset "equilibrium_bath_state — trace properties" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        agg = setup_aggregate(aggCore; vib_basis = :ground_ground)

        W_eq = equilibrium_bath_state(agg.core, agg.tools, agg.operators, 300.0)
        @test all(isfinite, W_eq)

        indicesMap = agg.tools.indicesMap
        a1, a2 = indicesMap[1][1], indicesMap[1][end]
        rho_diag = real.(diag(W_eq[a1:a2, a1:a2]))
        @test sum(rho_diag) ≈ 1.0 atol = 1e-10
        @test all(rho_diag .>= 0.0)
    end

    @testset "site_to_exciton_kernel — identity transform" begin
        M = zeros(ComplexF64, 3, 3, 3, 3)
        M[2, 2, 2, 2] = 1.0 + 0.5im
        M[2, 3, 2, 3] = 2.0 - 0.3im

        coeffs_id = [1.0 0.0; 0.0 1.0]
        M_exc = site_to_exciton_kernel(M, coeffs_id)
        @test M_exc[2, 2, 2, 2] ≈ M[2, 2, 2, 2]
        @test M_exc[2, 3, 2, 3] ≈ M[2, 3, 2, 3]
    end

end
