using Test
using OpenQuantumSystems
using LinearAlgebra
import QuadGK

# Exact finite-basis memory kernel, built without any package helpers:
# M_abcd(t1, t2) = tr_B <a| [H_I(t1), [H_I(t2), |c><d| ⊗ w_eq]] |b>
function _brute_memory_kernel(H1, H2, w_eq, indicesMap, elLen)
    n = size(H1, 1)
    M = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    for c in 1:elLen, d in 1:elLen
        X = zeros(ComplexF64, n, n)
        X[indicesMap[c], indicesMap[d]] .= w_eq
        Y = H2 * X - X * H2
        Z = H1 * Y - Y * H1
        for a in 1:elLen, b in 1:elLen
            M[a, b, c, d] = tr(Z[indicesMap[a], indicesMap[b]])
        end
    end
    return M
end

# Rotate the electronic part of a full-space operator from site to exciton basis
function _exciton_rotation(indicesMap, coefficients)
    n = indicesMap[end][end]
    nb = length(indicesMap[1])
    N = size(coefficients, 1)
    U = zeros(ComplexF64, n, n)
    U[indicesMap[1], indicesMap[1]] .= Matrix(I, nb, nb)
    for α in 1:N, k in 1:N
        U[indicesMap[α+1], indicesMap[k+1]] .= coefficients[α, k] .* Matrix(I, nb, nb)
    end
    return U
end

# Exact first-order kernel in exciton basis for a constant exciton-basis RDM R:
# M1_abcd = M0_abcd - tr_B <a|[H1,[H2, |c><d| ⊗ Σ_ef (R_ef/R_cd) ∫∫ <c|[H3,[H4, |e><f| ⊗ w]]|d>]]|b>
function _brute_first_order_kernel(agg, w_eq, U, R, t1, t2; rtol = 1e-7)
    ops, im_ = agg.operators, agg.tools.indicesMap
    n = size(U, 1)
    elLen = length(im_)
    HI(t) = U * get_interaction_ham_i_picture(ops.Ham_0.data, ops.Ham_I.data, t) * U'
    dc(H, X) = H * X - X * H
    lift(e, f, B) = (X = zeros(ComplexF64, n, n); X[im_[e], im_[f]] .= B; X)
    Z = Dict{Tuple{Int,Int},Matrix{ComplexF64}}()
    for e in 2:elLen, f in 2:elLen
        Z[(e, f)], _ = QuadGK.quadgk(0.0, t2; rtol = rtol) do t3
            H3 = HI(t3)
            inner, _ = QuadGK.quadgk(t4 -> dc(H3, dc(HI(t4), lift(e, f, w_eq))), 0.0, t3; rtol = rtol)
            inner
        end
    end
    H1, H2 = HI(t1), HI(t2)
    M = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    for c in 2:elLen, d in 2:elLen
        corr = sum(R[e, f] / R[c, d] .* Z[(e, f)][im_[c], im_[d]] for e in 2:elLen, f in 2:elLen)
        Y = dc(H1, dc(H2, lift(c, d, w_eq .- corr)))
        for a in 2:elLen, b in 2:elLen
            M[a, b, c, d] = tr(Y[im_[a], im_[b]])
        end
    end
    return M
end

function _validation_dimer(nvib; S = 0.02, J = 50.0, vib_basis = GroundGround())
    mode = Mode(; omega = 200.0, hr_factor = S)
    aggCore = AggregateCore([
        Molecule([mode], nvib, [0.0, 12500.0]),
        Molecule([mode], nvib, [0.0, 12700.0]),
    ])
    aggCore.coupling[2, 3] = J
    aggCore.coupling[3, 2] = J
    return setup_aggregate(aggCore; vib_basis = vib_basis)
end

@testset "kernel_validation" begin

    @testset "Huang-Rhys convention S = d^2/2" begin
        mode = Mode(; omega = 200.0, hr_factor = 0.05)
        @test mode.shift ≈ sqrt(2 * 0.05)

        mol = Molecule([mode], 3, [0.0, 12500.0])
        sd = spectral_density(mol)
        @test sd.huang_rhys[1] ≈ 0.05

        # <0|ΔV|1> = ω d / sqrt(2) = ω sqrt(S)
        agg = setup_aggregate(AggregateCore([mol]); vib_basis = GroundGround())
        i2 = agg.tools.indicesMap[2]
        @test agg.operators.Ham_I.data[i2[1], i2[2]] ≈ 200.0 * sqrt(0.05)
    end

    @testset "GroundGround ΔV leaves spectator modes untouched" begin
        agg = _validation_dimer(3; S = 0.05)
        H_I = agg.operators.Ham_I.data
        i2, i3 = agg.tools.indicesMap[2], agg.tools.indicesMap[3]
        vib = agg.tools.vibIndices
        for (k, I) in enumerate(i2), (l, K) in enumerate(i2)
            # site-1 excitation: mode of molecule 2 must not change
            if vib[k][2] != vib[l][2]
                @test H_I[I, K] == 0.0
            end
        end
        for (k, I) in enumerate(i3), (l, K) in enumerate(i3)
            if vib[k][1] != vib[l][1]
                @test H_I[I, K] == 0.0
            end
        end
    end

    @testset "GroundGround and GroundExcited spectra agree" begin
        nvib = 6
        gg = _validation_dimer(nvib; S = 0.05, vib_basis = GroundGround())
        ge = _validation_dimer(nvib; S = 0.05, vib_basis = GroundExcited())
        n = nvib^2
        e_gg = eigvals(Hermitian(Matrix(gg.operators.Ham.data)))[n+1:n+4]
        e_ge = eigvals(Hermitian(Matrix(ge.operators.Ham.data)))[n+1:n+4]
        # the two bases differ by a constant reference energy of the excited states
        @test e_gg .- e_gg[1] ≈ e_ge .- e_ge[1] atol = 1e-3
    end

    @testset "trace_bath_part keeps coherences in GroundGround" begin
        agg = _validation_dimer(2)
        W = rand(ComplexF64, 4, 4)
        @test trace_bath_part(W, 2, 3, agg.tools; vib_basis = GroundGround()) ≈ tr(W)
    end

    @testset "MemoryKernel_traced matches brute force" begin
        agg = _validation_dimer(3)
        ops, tools = agg.operators, agg.tools
        W_eq = equilibrium_bath_state(agg.core, tools, ops, 10.0)
        i1 = tools.indicesMap[1]
        H1 = get_interaction_ham_i_picture(ops.Ham_0.data, ops.Ham_I.data, 0.003)
        H2 = get_interaction_ham_i_picture(ops.Ham_0.data, ops.Ham_I.data, 0.001)
        M = MemoryKernel_traced(H1, H2, W_eq, agg.core, ops, tools)
        M_ref = _brute_memory_kernel(H1, H2, W_eq[i1, i1], tools.indicesMap, 3)
        @test M ≈ M_ref atol = 1e-9
    end

    @testset "zeroth_order_memory_kernel_cf matches exact finite-basis kernel" begin
        T = 10.0
        ex = 2:3
        for J in (0.0, 50.0), nvib in (3, 4)
            agg = _validation_dimer(nvib; J = J)
            ops, tools = agg.operators, agg.tools
            W_eq = equilibrium_bath_state(agg.core, tools, ops, T)
            i1 = tools.indicesMap[1]
            for (t1, t2) in ((0.003, 0.001), (0.01, 0.004))
                H1 = get_interaction_ham_i_picture(ops.Ham_0.data, ops.Ham_I.data, t1)
                H2 = get_interaction_ham_i_picture(ops.Ham_0.data, ops.Ham_I.data, t2)
                M_site = _brute_memory_kernel(H1, H2, W_eq[i1, i1], tools.indicesMap, 3)
                _, c = exciton_basis(agg.core, tools)
                M_exact = site_to_exciton_kernel(M_site, c)
                M_cf = zeroth_order_memory_kernel_cf(t1, t2, agg.core; T = T)
                scale = maximum(abs.(M_exact[ex, ex, ex, ex]))
                @test maximum(abs.(M_cf[ex, ex, ex, ex] .- M_exact[ex, ex, ex, ex])) < 1e-9 * scale
            end
        end
    end

    @testset "first_order_memory_kernel_cf matches exact finite-basis kernel" begin
        T = 10.0
        ex = 2:3
        for J in (0.0, 50.0)
            agg = _validation_dimer(3; S = 0.05, J = J)
            ops, tools = agg.operators, agg.tools
            W_eq = equilibrium_bath_state(agg.core, tools, ops, T)
            i1 = tools.indicesMap[1]
            _, c = exciton_basis(agg.core, tools)
            U = _exciton_rotation(tools.indicesMap, c)
            R = zeros(ComplexF64, 3, 3)
            R[1, 1] = 1.0
            R[2:3, 2:3] .= [0.6 0.2+0.1im; 0.2-0.1im 0.4]
            t1, t2 = 0.003, 0.002
            M_exact = _brute_first_order_kernel(agg, W_eq[i1, i1], U, R, t1, t2)
            M_cf = first_order_memory_kernel_cf(t1, t2, agg.core, R; T = T, rtol = 1e-8, atol = 1e-10)
            M0 = zeroth_order_memory_kernel_cf(t1, t2, agg.core; T = T)
            corr_scale = maximum(abs.(M_exact[ex, ex, ex, ex] .- M0[ex, ex, ex, ex]))
            @test corr_scale > 1.0
            @test maximum(abs.(M_cf[ex, ex, ex, ex] .- M_exact[ex, ex, ex, ex])) < 1e-6 * corr_scale
        end
    end

end
