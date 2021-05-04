using Test
using OpenQuantumSystems
using SparseArrays, LinearAlgebra

@testset "metrics" begin

    b1 = SpinBasis(1 // 2)
    b2 = FockBasis(6)

    psi1 = spinup(b1) ⊗ coherentstate(b2, 0.1)
    psi2 = spindown(b1) ⊗ fockstate(b2, 2)

    rho = tensor(psi1, dagger(psi1))
    sigma = tensor(psi2, dagger(psi2))
    sigma_s = Commutator(sigma)
    rho_s = Commutator(rho)

    # tracenorm
    @test tracenorm(0 * rho_s) ≈ 0.0
    @test tracenorm_h(0 * rho_s) ≈ 0.0
    @test tracenorm_nh(0 * rho_s) ≈ 0.0

    @test tracenorm(rho_s) ≈ 26.0 # 1 * 6 * 4
    @test tracenorm_h(rho_s) ≈ 26.0
    @test tracenorm_nh(rho_s) ≈ 26.0

    # tracedistance
    @test tracedistance(rho_s, sigma_s) ≈ 26.0
    @test tracedistance_h(rho_s, sigma_s) ≈ 26.0
    @test tracedistance_nh(rho_s, sigma_s) ≈ 26.0

    @test tracedistance(rho_s, rho_s) ≈ 0.0
    @test tracedistance_h(rho_s, rho_s) ≈ 0.0
    @test tracedistance_nh(rho_s, rho_s) ≈ 0.0

    @test tracedistance(sigma_s, sigma_s) ≈ 0.0
    @test tracedistance_h(sigma_s, sigma_s) ≈ 0.0
    @test tracedistance_nh(sigma_s, sigma_s) ≈ 0.0

    # tracedistance
    @test tracedistance(rho_s, sigma_s) ≈ 26.0
    @test tracedistance(rho_s, rho_s) ≈ 0.0
    @test tracedistance(sigma_s, sigma_s) ≈ 0.0

    @test tracedistance_nh(rho_s, rho_s) ≈ 0.0

end # testset
