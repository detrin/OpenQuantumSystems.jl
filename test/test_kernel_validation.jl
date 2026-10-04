using Test
using OpenQuantumSystems
using LinearAlgebra

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

end
