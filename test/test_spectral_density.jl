using Test
using OpenQuantumSystems

@testset "spectral_density" begin

    @testset "SpectralDensity constructor" begin
        sd = SpectralDensity([200.0, 400.0], [0.02, 0.05])
        @test length(sd.frequencies) == 2
        @test sd.frequencies == [200.0, 400.0]
        @test sd.huang_rhys == [0.02, 0.05]

        @test_throws DimensionMismatch SpectralDensity([200.0], [0.02, 0.05])
    end

    @testset "spectral_density from Molecule" begin
        mode1 = Mode(; omega = 200.0, hr_factor = 0.02)
        mode2 = Mode(; omega = 400.0, hr_factor = 0.05)
        mol = Molecule([mode1, mode2], 3, [0.0, 12500.0])

        sd = spectral_density(mol)
        @test sd.frequencies == [200.0, 400.0]
        @test sd.huang_rhys ≈ [0.02, 0.05]
    end

    @testset "reorganization_energy" begin
        sd = SpectralDensity([200.0, 400.0], [0.02, 0.05])
        λ = reorganization_energy(sd)
        @test λ ≈ 200.0 * 0.02 + 400.0 * 0.05
    end

    @testset "lineshape_function at t=0" begin
        sd = SpectralDensity([200.0], [0.1])
        g = lineshape_function(sd, 0.0, 0.0)
        @test g ≈ 0.0 + 0.0im atol=1e-15
    end

    @testset "lineshape_function T=0" begin
        sd = SpectralDensity([200.0], [0.1])
        t = 0.01
        g = lineshape_function(sd, t, 0.0)
        ωt = 200.0 * t
        expected = 0.1 * ((1.0 - cos(ωt)) + im * (sin(ωt) - ωt))
        @test g ≈ expected atol=1e-12
    end

    @testset "lineshape_function T>0" begin
        sd = SpectralDensity([200.0], [0.1])
        t = 0.01
        T = 300.0
        g = lineshape_function(sd, t, T)
        kT = OpenQuantumSystems.BOLTZMANN_CM * T
        n = 1.0 / (exp(200.0 / kT) - 1.0)
        ωt = 200.0 * t
        expected = 0.1 * ((2n + 1) * (1.0 - cos(ωt)) + im * (sin(ωt) - ωt))
        @test g ≈ expected atol=1e-12
    end

    @testset "lineshape_derivative at t=0" begin
        sd = SpectralDensity([200.0], [0.1])
        gd = lineshape_derivative(sd, 0.0, 0.0)
        @test gd ≈ 0.0 + 0.0im atol=1e-15
    end

    @testset "lineshape_second_derivative at t=0" begin
        sd = SpectralDensity([200.0], [0.1])
        gdd = lineshape_second_derivative(sd, 0.0, 0.0)
        @test real(gdd) ≈ 0.1 * 200.0^2 atol=1e-10
        @test imag(gdd) ≈ 0.0 atol=1e-15
    end

    @testset "lineshape numerical derivative consistency" begin
        sd = SpectralDensity([200.0, 400.0], [0.1, 0.05])
        t = 0.005
        T = 300.0
        dt = 1e-7

        gd_analytic = lineshape_derivative(sd, t, T)
        gd_numerical = (lineshape_function(sd, t + dt, T) - lineshape_function(sd, t - dt, T)) / (2dt)
        @test gd_analytic ≈ gd_numerical atol=1e-4

        gdd_analytic = lineshape_second_derivative(sd, t, T)
        gdd_numerical = (lineshape_derivative(sd, t + dt, T) - lineshape_derivative(sd, t - dt, T)) / (2dt)
        @test gdd_analytic ≈ gdd_numerical atol=1e-2
    end

end
