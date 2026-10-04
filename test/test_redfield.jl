using Test
include("helpers_qme.jl")

@testset "redfield" begin

    tspan = QME_TEST_TSPAN
    kw = qme_test_kwargs()

    @testset "no system-bath coupling: no dynamics in interaction picture" begin
        agg, W0 = qme_test_dimer(0.0)
        _, rho_t = QME_sI_Redfield(W0, tspan, agg; kw...)
        @test maximum(norm(rho.data - rho_t[1].data) for rho in rho_t) < 1e-12
    end

    @testset "weak coupling: close to exact, trace and hermiticity preserved" begin
        agg, W0 = qme_test_dimer(0.05)
        _, rho_t = QME_sI_Redfield(W0, tspan, agg; kw...)
        ref = qme_exact(agg, W0, tspan)
        @test qme_relative_error(rho_t, ref) < 0.1
        @test qme_trace_drift(rho_t) < 1e-10
        @test qme_hermiticity(rho_t) < 1e-10
    end

    @testset "default DDE algorithm is constructible" begin
        agg, W0 = qme_test_dimer(0.05)
        (; reltol, abstol, int_reltol, int_abstol) = kw
        _, rho_t = QME_sI_Redfield(W0, tspan[1:2], agg;
            reltol = reltol, abstol = abstol, int_reltol = int_reltol, int_abstol = int_abstol)
        @test length(rho_t) == 2
    end

    @testset "accepts thermal_state_composite initial state" begin
        agg, _ = qme_test_dimer(0.05)
        W0 = thermal_state_composite(300.0, [0.0, 1.0, 0.0], agg)
        _, rho_t = QME_sI_Redfield(W0, tspan[1:2], agg; kw...)
        @test qme_trace_drift(rho_t) < 1e-10
    end

    # Second-order theory: relative error ∝ S (observed ratio 4.1 for S → S/4).
    @testset "error scales linearly with Huang-Rhys factor" begin
        errs = map((0.05, 0.0125)) do S
            agg, W0 = qme_test_dimer(S)
            _, rho_t = QME_sI_Redfield(W0, tspan, agg; kw...)
            qme_relative_error(rho_t, qme_exact(agg, W0, tspan))
        end
        @test 3.0 < errs[1] / errs[2] < 5.0
    end

end # testset
