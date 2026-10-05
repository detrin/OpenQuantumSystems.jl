using Test
using OpenQuantumSystems
using LinearAlgebra
import QuadGK

const _OQS = OpenQuantumSystems

function _iqme_dimer(nvib; S = 0.05, J = 50.0)
    mode = Mode(; omega = 200.0, hr_factor = S)
    aggCore = AggregateCore([Molecule([mode], nvib, [0.0, 12500.0]), Molecule([mode], nvib, [0.0, 12700.0])])
    aggCore.coupling[2, 3] = J
    aggCore.coupling[3, 2] = J
    return aggCore
end

# kernel value Σ_g coeff e^{iΩ t1} e^{-Ẑ(t1-t2)} B(t2) of convolution terms, indexed [a,b,c,d]
function _iqme_kernel(terms, t1, t2, i2, N)
    K = zeros(ComplexF64, N, N, N, N)
    for g in terms
        a = mod1(g.out, N); b = div(g.out - 1, N) + 1
        c = mod1(g.inp, N); d = div(g.inp - 1, N) + 1
        B = g.B === nothing ? one(ComplexF64) : g.B[i2]
        K[a, b, c, d] += g.coeff * cis(g.Ω * t1) * exp(-g.Ẑ * (t1 - t2)) * B
    end
    return K
end

@testset "iterative_qme" begin

    @testset "bath correlation functions" begin
        sd = SpectralDensity([200.0, 350.0], [0.05, 0.02])
        C = bath_correlation(sd, 300.0)
        for τ in (-0.02, 0.0, 0.003, 0.05)
            @test C(τ) ≈ analytic_correlation_nn(sd, τ, 300.0) atol = 1e-9
        end

        dl = DrudeLorentz(20.0, 100.0)
        Cdl = bath_correlation(dl, 300.0; n_matsubara = 50)
        # Im C(τ) = -λγ e^{-γτ} exactly
        @test imag(Cdl(0.01)) ≈ -20.0 * 100.0 * exp(-1.0)
        @test Cdl(-0.01) ≈ conj(Cdl(0.01))
        # real part against the spectral integral (tail beyond 4e4 cm⁻¹ < 1 %)
        β = 1 / (_OQS.BOLTZMANN_CM * 300.0)
        τ = 0.01
        Cnum = QuadGK.quadgk(ω -> _OQS.spectral_density(dl, ω) / π * coth(β * ω / 2) * cos(ω * τ),
            collect(0.0:20.0:4e4)...; rtol = 1e-8)[1]
        @test real(Cdl(τ)) ≈ Cnum rtol = 1e-2
        @test_throws ArgumentError bath_correlation(dl, 0.0)
        @test_throws ArgumentError ExponentialCorrelation([1.0], [-1.0])
    end

    @testset "ExcitonModel" begin
        aggCore = _iqme_dimer(3)
        model = ExcitonModel(aggCore; T = 10.0)
        E, c = exciton_basis(aggCore, AggregateTools(aggCore))
        @test model.energies ≈ E
        @test model.coefficients * model.coefficients' ≈ I
        ρ = ComplexF64[0.7 0.1im; -0.1im 0.3]
        @test exciton_to_site(model, site_to_exciton(model, ρ)) ≈ ρ
        bath = bath_correlation(DrudeLorentz(10.0, 100.0), 300.0)
        @test_throws DimensionMismatch ExcitonModel([1.0 0.0; 0.0 2.0], [bath])
    end

    aggCore = _iqme_dimer(3)
    T = 10.0
    model = ExcitonModel(aggCore; T = T)

    @testset "second-order kernel equals zeroth_order_memory_kernel_cf" begin
        terms = _OQS._second_order_terms(model)
        for (t1, t2) in ((0.003, 0.001), (0.02, 0.005))
            M = -_iqme_kernel(terms, t1, t2, 1, 2)   # terms store -M⁰
            Mref = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = T)[2:3, 2:3, 2:3, 2:3]
            @test maximum(abs.(M .- Mref)) < 1e-10 * maximum(abs.(Mref))
        end
    end

    @testset "fourth-order (thesis) kernel equals first_order_memory_kernel_cf" begin
        R = ComplexF64[0.6 0.2+0.1im; 0.2-0.1im 0.4]
        h = 1e-5
        nt = 301
        ρc = zeros(ComplexF64, nt, 2, 2)
        for i in 1:nt
            ρc[i, :, :] .= R
        end
        terms = _OQS._fourth_order_terms(model, ρc, h, :thesis)
        Rfull = zeros(ComplexF64, 3, 3); Rfull[1, 1] = 1; Rfull[2:3, 2:3] .= R
        t1, i2 = 0.004, nt
        t2 = (i2 - 1) * h
        ΔM = _iqme_kernel(terms, t1, t2, i2, 2)
        M0 = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = T)[2:3, 2:3, 2:3, 2:3]
        M1 = first_order_memory_kernel_cf(t1, t2, aggCore, Rfull; T = T, rtol = 1e-9, atol = 1e-11)[2:3, 2:3, 2:3, 2:3]
        ref = M0 .- M1
        @test maximum(abs.(ΔM .- ref)) < 1e-5 * maximum(abs.(ref))
        @test_throws ArgumentError _OQS._fourth_order_terms(model, ρc, h, :unknown)
    end

    @testset "disconnected part of thesis kernel is -M⁰ ln(ρ(τ)/ρ(0))" begin
        λ = 20.0
        bath = bath_correlation(DrudeLorentz(λ, 100.0), 300.0)
        m = ExcitonModel([12500.0+λ 100.0; 100.0 12700.0+λ], [bath, bath])
        dt = 2.5e-4
        born = QME_cf_born(m, site_to_exciton(m, ComplexF64[1 0; 0 0]), 0.1, dt)
        Tth = _OQS._fourth_order_terms(m, born.ρ, dt, :thesis)
        Tpr = _OQS._fourth_order_terms(m, born.ρ, dt, :projected)
        T2 = _OQS._second_order_terms(m)
        i2 = 301; i1 = i2 + 8
        t1, t2 = (i1 - 1) * dt, (i2 - 1) * dt
        D = _iqme_kernel(Tth, t1, t2, i2, 2) - _iqme_kernel(Tpr, t1, t2, i2, 2)
        M0 = -_iqme_kernel(T2, t1, t2, i2, 2)
        for c in 1:2
            pred = -M0[:, :, c, c] * log(real(born.ρ[i2, c, c]) / real(born.ρ[1, c, c]))
            @test norm(D[:, :, c, c] - pred) < 1e-3 * norm(D[:, :, c, c])
        end
    end

    @testset "solvers match finite-basis QME with the same equations" begin
        agg = setup_aggregate(aggCore; vib_basis = GroundGround())
        ops, tools = agg.operators, agg.tools
        W_eq = equilibrium_bath_state(agg.core, tools, ops, T)
        ρsite = zeros(ComplexF64, 3, 3); ρsite[2, 2] = 1.0
        W0 = DenseOperator(tools.basis, tools.basis, ad(ρsite, W_eq, agg.core, tools))
        ρ0 = site_to_exciton(model, ρsite[2:3, 2:3])
        t_max, dt = 0.01, 2e-5
        tspan = collect(range(0, t_max; length = 3))
        kw = (reltol = 1e-8, abstol = 1e-10, int_reltol = 1e-8, int_abstol = 1e-10)
        at(res, t) = res.ρ[round(Int, t / dt) + 1, :, :]
        to_exc(r) = site_to_exciton(model, r.data[2:3, 2:3])

        born = QME_cf_born(model, ρ0, t_max, dt)
        _, ref = QME_sI_ansatz(W0, tspan, agg; ansatz = :const_int, kw...)
        change = maximum(norm(to_exc(r) - ρ0) for r in ref)
        @test maximum(norm(at(born, t) - to_exc(r)) for (t, r) in zip(tspan, ref)) < 1e-5 * change

        red = QME_cf_redfield(model, ρ0, t_max, dt)
        _, refR = QME_sI_Redfield(W0, tspan, agg; kw...)
        @test maximum(norm(at(red, t) - to_exc(r)) for (t, r) in zip(tspan, refR)) < 1e-5 * change

        for res in (born, red)
            @test all(abs(tr(res.ρ[i, :, :]) - 1) < 1e-10 for i in eachindex(res.t))
            @test all(norm(res.ρ[i, :, :] - res.ρ[i, :, :]') < 1e-10 for i in eachindex(res.t))
        end
    end

    # The corrected iterative QME is a fourth-order theory: its error relative to the
    # change of ρ falls as S², the second-order QMEs only as S.
    @testset "accuracy order against exact finite-bath dynamics" begin
        function case(S)
            ac = _iqme_dimer(4; S = S)
            agg = setup_aggregate(ac; vib_basis = GroundGround())
            ops, tools = agg.operators, agg.tools
            W_eq = equilibrium_bath_state(agg.core, tools, ops, T)
            ρsite = zeros(ComplexF64, 3, 3); ρsite[2, 2] = 0.7; ρsite[3, 3] = 0.3
            W0 = ad(ρsite, W_eq, agg.core, tools)
            m = ExcitonModel(ac; T = T)
            ρ0 = site_to_exciton(m, ρsite[2:3, 2:3])
            t_max, dt = 0.05, 2.5e-4
            n = round(Int, t_max / dt)
            exact = map(0:n) do i
                t = i * dt
                U = exp(-im * Matrix(ops.Ham.data) * t); U0 = exp(-im * Matrix(ops.Ham_0.data) * t)
                site_to_exciton(m, trace_bath(U0' * U * W0 * U' * U0, agg.core, ops, tools)[2:3, 2:3])
            end
            change = maximum(norm(x - ρ0) for x in exact)
            err(ρ) = maximum(norm(ρ[i, :, :] - exact[i]) for i in eachindex(exact)) / change
            return (born = err(QME_cf_born(m, ρ0, t_max, dt).ρ),
                    nz = err(QME_cf_iterative(m, ρ0, t_max, dt; correction = :nz).ρ),
                    projected = err(QME_cf_iterative(m, ρ0, t_max, dt; correction = :projected).ρ))
        end
        weak, weaker = case(0.0125), case(0.00625)
        @test weak.nz < 0.3 * weak.born
        @test 0.8 < log2(weak.born / weaker.born) < 1.3
        @test log2(weak.nz / weaker.nz) > 1.7
        @test log2(weak.projected / weaker.projected) > 1.7
    end

end
