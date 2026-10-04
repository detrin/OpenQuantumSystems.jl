using Test
include("helpers_qme.jl")

# With the exact W(s) inside the memory integral, the second-order QME is an
# identity, so every *_exact solver must reproduce the exact dynamics.
@testset "master" begin

    agg, W0 = qme_test_dimer(0.05)
    tspan = QME_TEST_TSPAN
    kw = qme_test_kwargs()

    @testset "QME_sI_exact equals exact reduced dynamics" begin
        _, rho_t = QME_sI_exact(W0, tspan, agg; kw...)
        ref = qme_exact(agg, W0, tspan; picture = :interaction)
        @test qme_relative_error(rho_t, ref) < 1e-3
        @test qme_trace_drift(rho_t) < 1e-10
        @test qme_hermiticity(rho_t) < 1e-10
    end

    @testset "QME_sS_exact equals exact reduced dynamics" begin
        _, rho_t = QME_sS_exact(W0, tspan, agg; kw...)
        ref = qme_exact(agg, W0, tspan; picture = :schrodinger)
        @test qme_relative_error(rho_t, ref) < 1e-3
        @test qme_trace_drift(rho_t) < 1e-10
        @test qme_hermiticity(rho_t) < 1e-10
    end

    # The full-state solvers feed their own history back into the integral, so
    # their accuracy is limited by the solver tolerance (≈ 3e-3 at 1e-4).
    @testset "QME_SI_exact equals exact full dynamics" begin
        _, W_t = QME_SI_exact(W0, tspan, agg; kw...)
        ref = qme_exact(agg, W0, tspan; picture = :interaction, reduced = false)
        @test qme_relative_error(W_t, ref) < 2e-2
        @test qme_trace_drift(W_t) < 1e-10
        @test qme_hermiticity(W_t) < 1e-10
    end

    @testset "QME_SS_exact equals exact full dynamics" begin
        _, W_t = QME_SS_exact(W0, tspan, agg; kw...)
        ref = qme_exact(agg, W0, tspan; picture = :schrodinger, reduced = false)
        @test qme_relative_error(W_t, ref) < 2e-2
        @test qme_trace_drift(W_t) < 1e-10
        @test qme_hermiticity(W_t) < 1e-10
    end

end # testset
