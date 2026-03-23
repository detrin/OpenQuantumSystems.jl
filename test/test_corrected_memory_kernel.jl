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

    @testset "first_order — basic properties" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        J = 100.0
        aggCore.coupling[2, 3] = J
        aggCore.coupling[3, 2] = J

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        t1, t2 = 0.05, 0.02
        M = first_order_memory_kernel_cf(t1, t2, aggCore, rho_I;
            T = 300.0, rtol = 1e-4, atol = 1e-6)

        @test size(M) == (3, 3, 3, 3)
        @test all(isfinite, M)

        for i in 1:3
            @test all(abs.(M[1, :, :, :]) .< 1e-12)
            @test all(abs.(M[:, 1, :, :]) .< 1e-12)
            @test all(abs.(M[:, :, 1, :]) .< 1e-12)
            @test all(abs.(M[:, :, :, 1]) .< 1e-12)
        end
    end

    @testset "first_order — hermiticity" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        M = first_order_memory_kernel_cf(0.05, 0.02, aggCore, rho_I;
            T = 300.0, rtol = 1e-4, atol = 1e-6)
        for a in 2:elLen, b in 2:elLen, c in 2:elLen, d in 2:elLen
            @test M[a, b, c, d] ≈ conj(M[b, a, d, c]) atol = 1e-6
        end
    end

    @testset "first_order — t2=0 recovers zeroth order" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        t1 = 0.1
        M0 = zeroth_order_memory_kernel_cf(t1, 0.0, aggCore; T = 300.0)
        M1 = first_order_memory_kernel_cf(t1, 0.0, aggCore, rho_I; T = 300.0)
        @test M1 ≈ M0 atol = 1e-10
    end

    @testset "first_order — J=0 diagonal structure" begin
        mode1 = Mode(; omega = 200.0, hr_factor = 0.02)
        mode2 = Mode(; omega = 300.0, hr_factor = 0.01)
        mol1 = Molecule([mode1], 2, [0.0, 12500.0])
        mol2 = Molecule([mode2], 2, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.55
        rho_I[3, 3] = 0.45

        M = first_order_memory_kernel_cf(0.05, 0.02, aggCore, rho_I;
            T = 300.0, rtol = 1e-4, atol = 1e-6)

        for a in 2:elLen, b in 2:elLen, c in 2:elLen, d in 2:elLen
            if a != c || b != d
                @test abs(M[a, b, c, d]) < 1e-8
            end
        end
    end

    @testset "first_order — differs from zeroth order" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.05)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        t1, t2 = 0.1, 0.05
        M0 = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = 300.0)
        M1 = first_order_memory_kernel_cf(t1, t2, aggCore, rho_I;
            T = 300.0, rtol = 1e-4, atol = 1e-6)

        diff = maximum(abs.(M1 .- M0))
        @test diff > 1e-10
    end

    @testset "corrected_rates_cf — zeroth order" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        rates = corrected_rates_cf(aggCore;
            T = 300.0, t_ref = 0.5, order = 0, rtol = 1e-3, atol = 1e-5)

        @test size(rates) == (2, 2)
        @test all(isfinite, rates)

        for n in 1:2
            @test abs(sum(rates[:, n])) < 1e-6
        end

        @test rates[1, 2] != 0.0
        @test rates[2, 1] != 0.0
    end

    @testset "corrected_rates_cf — first order (monomer)" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.05)
        mol = Molecule([mode], 2, [0.0, 12500.0])
        aggCore = AggregateCore([mol])

        rates = corrected_rates_cf(aggCore;
            T = 300.0, t_ref = 0.05, order = 1, rtol = 1e-2, atol = 1e-4)

        @test size(rates) == (1, 1)
        @test all(isfinite, rates)
        @test abs(rates[1, 1]) < 1e-6
    end

    @testset "corrected_rates_cf — first order differs from zeroth" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.05)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        r0 = corrected_rates_cf(aggCore;
            T = 300.0, t_ref = 0.05, order = 0, rtol = 1e-3, atol = 1e-5)
        r1 = corrected_rates_cf(aggCore;
            T = 300.0, t_ref = 0.05, order = 1, rtol = 1e-2, atol = 1e-4)

        @test size(r1) == (2, 2)
        @test all(isfinite, r1)
        @test maximum(abs.(r1 .- r0)) > 1e-10
    end

    @testset "corrected_qme_rdm — basic" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        rho0 = [1.0, 0.0]
        tspan = collect(0.0:0.01:0.1)

        t_out, pop_final, all_pops = corrected_qme_rdm(aggCore, rho0, tspan;
            T = 300.0, t_ref = 0.05, max_iter = 2, rtol = 1e-2, atol = 1e-4)

        @test length(t_out) == length(tspan)
        @test size(pop_final) == (length(tspan), 2)
        @test length(all_pops) == 2

        @test pop_final[1, 1] ≈ 1.0 atol = 1e-10
        @test pop_final[1, 2] ≈ 0.0 atol = 1e-10

        for i in 1:length(tspan)
            @test sum(pop_final[i, :]) ≈ 1.0 atol = 1e-6
            @test all(pop_final[i, :] .>= 0.0)
        end
    end

    @testset "corrected_qme_rdm — iterations differ" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.05)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        rho0 = [1.0, 0.0]
        tspan = collect(0.0:0.01:0.1)

        _, _, all_pops = corrected_qme_rdm(aggCore, rho0, tspan;
            T = 300.0, t_ref = 0.05, max_iter = 2, rtol = 1e-2, atol = 1e-4)

        diff = maximum(abs.(all_pops[2] .- all_pops[1]))
        @test diff > 1e-10
    end

    @testset "analytic_correlation_nn_high_T — real valued" begin
        sd = SpectralDensity([200.0], [0.02])
        T_high = 50000.0
        tau = 0.1

        C_ht = analytic_correlation_nn_high_T(sd, tau, T_high)
        @test isa(C_ht, Float64)
        @test isfinite(C_ht)

        C_full = analytic_correlation_nn(sd, tau, T_high)
        @test abs(imag(C_full)) / abs(real(C_full)) < 1e-2
        @test abs(C_ht - real(C_full)) / abs(real(C_full)) < 0.005
    end

    @testset "analytic_correlation_nn_high_T — symmetry C(τ)=C(-τ)" begin
        sd = SpectralDensity([200.0, 400.0], [0.02, 0.01])
        T = 5000.0
        tau = 0.1

        C_fwd = analytic_correlation_nn_high_T(sd, tau, T)
        C_bwd = analytic_correlation_nn_high_T(sd, -tau, T)
        @test C_fwd ≈ C_bwd atol = 1e-12
    end

    @testset "is_high_temperature" begin
        sd_low = SpectralDensity([1000.0], [0.1])
        @test is_high_temperature(sd_low, 300.0) == false
        @test is_high_temperature(sd_low, 10000.0) == true

        sd_high = SpectralDensity([10.0], [0.1])
        @test is_high_temperature(sd_high, 300.0) == true

        @test is_high_temperature(sd_low, 0.0) == false
    end

    @testset "zeroth_order_memory_kernel_high_T — basic properties" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        t1, t2 = 0.1, 0.05
        M = zeroth_order_memory_kernel_high_T(t1, t2, aggCore; T = 10000.0)
        @test size(M) == (3, 3, 3, 3)
        @test all(isfinite, M)

        for i in 1:3
            @test all(abs.(M[1, :, :, :]) .< 1e-15)
            @test all(abs.(M[:, 1, :, :]) .< 1e-15)
            @test all(abs.(M[:, :, 1, :]) .< 1e-15)
            @test all(abs.(M[:, :, :, 1]) .< 1e-15)
        end
    end

    @testset "zeroth_order — high-T hermiticity" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        M = zeroth_order_memory_kernel_high_T(0.1, 0.05, aggCore; T = 50000.0)
        for a in 2:3, b in 2:3, c in 2:3, d in 2:3
            @test M[a, b, c, d] ≈ conj(M[b, a, d, c]) atol = 1e-10
        end
    end

    @testset "zeroth_order — high-T converges to full at extreme T" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        T_extreme = 500000.0
        t1, t2 = 0.1, 0.05
        M_full = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = T_extreme)
        M_ht = zeroth_order_memory_kernel_high_T(t1, t2, aggCore; T = T_extreme)

        max_err = maximum(abs.(M_ht .- M_full))
        max_val = maximum(abs.(M_full))
        @test max_err / max_val < 0.1
    end

    @testset "first_order_memory_kernel_high_T — basic properties" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        t1, t2 = 0.05, 0.02
        M = first_order_memory_kernel_high_T(t1, t2, aggCore, rho_I;
            T = 50000.0, rtol = 1e-4, atol = 1e-6)

        @test size(M) == (3, 3, 3, 3)
        @test all(isfinite, M)

        for i in 1:3
            @test all(abs.(M[1, :, :, :]) .< 1e-12)
            @test all(abs.(M[:, 1, :, :]) .< 1e-12)
            @test all(abs.(M[:, :, 1, :]) .< 1e-12)
            @test all(abs.(M[:, :, :, 1]) .< 1e-12)
        end
    end

    @testset "first_order — high-T t2=0 recovers zeroth order" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        t1 = 0.1
        M0 = zeroth_order_memory_kernel_high_T(t1, 0.0, aggCore; T = 50000.0)
        M1 = first_order_memory_kernel_high_T(t1, 0.0, aggCore, rho_I; T = 50000.0)
        @test M1 ≈ M0 atol = 1e-10
    end

    @testset "first_order — high-T matches full at high temperature" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.02)
        mol1 = Molecule([mode], 3, [0.0, 12500.0])
        mol2 = Molecule([mode], 3, [0.0, 12700.0])
        aggCore = AggregateCore([mol1, mol2])
        aggCore.coupling[2, 3] = 100.0
        aggCore.coupling[3, 2] = 100.0

        elLen = 3
        rho_I = zeros(ComplexF64, elLen, elLen)
        rho_I[1, 1] = 1.0
        rho_I[2, 2] = 0.6
        rho_I[3, 3] = 0.4

        T_high = 50000.0
        t1, t2 = 0.05, 0.02
        M_full = first_order_memory_kernel_cf(t1, t2, aggCore, rho_I;
            T = T_high, rtol = 1e-4, atol = 1e-6)
        M_ht = first_order_memory_kernel_high_T(t1, t2, aggCore, rho_I;
            T = T_high, rtol = 1e-4, atol = 1e-6)

        for a in 2:3, b in 2:3, c in 2:3, d in 2:3
            ref = abs(M_full[a,b,c,d])
            ref < 1e-12 && continue
            @test abs(M_ht[a,b,c,d] - M_full[a,b,c,d]) / ref < 0.05
        end
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
