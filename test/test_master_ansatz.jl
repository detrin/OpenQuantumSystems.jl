using Test
include("helpers_qme.jl")

@testset "master ansatz" begin

    agg, W0 = qme_test_dimer(0.05)
    tspan = QME_TEST_TSPAN
    kw = qme_test_kwargs()
    ref = qme_exact(agg, W0, tspan)

    # :test feeds the exact bath part into the kernel, so the QME is exact.
    @testset ":test ansatz equals exact dynamics" begin
        _, rho_t = QME_sI_ansatz(W0, tspan, agg; ansatz = :test, kw...)
        @test qme_relative_error(rho_t, ref) < 1e-3
    end

    # Observed relative errors at S = 0.05: const_int 0.06, linear2_sch 0.06,
    # upart*_int 0.01; the Schrödinger-picture ansatzes are much worse (0.2–0.33).
    accurate = (:const_int, :linear2_sch, :upart1_int, :upart2_int)
    rough = (:const_sch, :linear_sch, :upart1_sch, :upart2_sch)
    for ansatz in (accurate..., rough...)
        @testset "$ansatz ansatz" begin
            _, rho_t = QME_sI_ansatz(W0, tspan, agg; ansatz = ansatz, kw...)
            @test qme_relative_error(rho_t, ref) < (ansatz in accurate ? 0.1 : 0.5)
            @test qme_trace_drift(rho_t) < 1e-10
            @test qme_hermiticity(rho_t) < 1e-10
        end
    end

    @testset "const_int error scales linearly with Huang-Rhys factor" begin
        agg_w, W0_w = qme_test_dimer(0.0125)
        _, rho_strong = QME_sI_ansatz(W0, tspan, agg; ansatz = :const_int, kw...)
        _, rho_weak = QME_sI_ansatz(W0_w, tspan, agg_w; ansatz = :const_int, kw...)
        ratio = qme_relative_error(rho_strong, ref) /
                qme_relative_error(rho_weak, qme_exact(agg_w, W0_w, tspan))
        @test 3.0 < ratio < 5.0
    end

    @testset "named wrappers match QME_sI_ansatz" begin
        wrappers = (
            QME_sI_ansatz_test => :test,
            QME_sI_ansatz_const_int => :const_int,
            QME_sI_ansatz_const_sch => :const_sch,
            QME_sI_ansatz_linear_sch => :linear_sch,
            QME_sI_ansatz_linear2_sch => :linear2_sch,
            QME_sI_ansatz_upart1_sch => :upart1_sch,
            QME_sI_ansatz_upart1_int => :upart1_int,
            QME_sI_ansatz_upart2_sch => :upart2_sch,
            QME_sI_ansatz_upart2_int => :upart2_int,
        )
        short = tspan[1:2]
        for (f, ansatz) in wrappers
            _, a = f(W0, short, agg; kw...)
            _, b = QME_sI_ansatz(W0, short, agg; ansatz = ansatz, kw...)
            @test a[end].data ≈ b[end].data
        end
    end

    @testset "accepts thermal_state_composite initial state" begin
        W0_th = thermal_state_composite(300.0, [0.0, 1.0, 0.0], agg)
        _, rho_t = QME_sI_ansatz(W0_th, tspan[1:2], agg; ansatz = :const_int, kw...)
        @test qme_trace_drift(rho_t) < 1e-10
    end

    @testset "unknown ansatz throws" begin
        @test_throws ArgumentError QME_sI_ansatz(W0, tspan, agg; ansatz = :nonexistent)
    end

end # testset
