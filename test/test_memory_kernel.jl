using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra
import QuantumOpticsBase


@testset "memory kernel" begin

    Random.seed!(0)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    HR = 0.01
    shift = (2.0 * HR)^0.5
    mode1 = Mode(300.0, shift)
    Energy = [12500.0, 12800.0]
    mol1 = Molecule([mode1], 2, Energy)
    mol2 = Molecule([mode1], 2, Energy)
    agg = Aggregate([mol1, mol2])
    agg.coupling[1, 2] = 100
    agg.coupling[2, 1] = 100
    aggInds = getIndices(agg)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds)
    FCProd = getFCProd(agg, FCFact, aggInds, vibindices)
    Ham = getAggHamiltonian(agg, aggInds, FCFact)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg)
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggInds, FCFact)
    Ham_S = Ham - Ham_int

    # T = 300
    # W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, aggInds; diagonalize=false, diagonal=true)
    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    W0 = DenseOperator(basis, basis, data)
    normalize!(W0)
    W0_bath = get_rho_bath(
        W0,
        agg,
        FCProd,
        aggInds,
        vibindices;
        justCopy = true,
    )
    rho0 = trace_bath(W0.data, agg, FCProd, aggInds, vibindices)

    W_ab = take_el_part(W0.data, 1, 1, vibindices)
    W_ab_len = size(W_ab, 1)
    @test 1e-15 > D(W_ab, W0.data[1:W_ab_len, 1:W_ab_len])

    H_lambda, H_S = eigen(Ham_S.data)
    H_Sinv = inv(H_S)
    H_lambda = diagm(H_lambda)

    tspan = [0.0:0.1:0.3;]
    for t_i = 1:length(tspan)
        for s_i = 1:t_i
            t = tspan[t_i]
            s = tspan[s_i]
            # println(t, " ", s)
            H_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
            H_II_s = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, s)

            MemoryKernel_1 = MemoryKernel_1_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            rho = zero(complex(rho0))
            for a = 1:size(rho, 1),
                b = 1:size(rho, 2),
                c = 1:size(rho, 1),
                d = 1:size(rho, 2)

                rho[a, b] += MemoryKernel_1[a, b, c, d] * rho0[c, d]
            end
            MemoryKernel_1_ref = H_II_t * H_II_s * W0.data
            rho_ref = trace_bath(
                MemoryKernel_1_ref,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))

            MemoryKernel_2 = MemoryKernel_2_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            rho = zero(complex(rho0))
            for a = 1:size(rho, 1),
                b = 1:size(rho, 2),
                c = 1:size(rho, 1),
                d = 1:size(rho, 2)

                rho[a, b] += MemoryKernel_2[a, b, c, d] * rho0[c, d]
            end
            MemoryKernel_2_ref = H_II_t * W0.data * H_II_s
            rho_ref = trace_bath(
                MemoryKernel_2_ref,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))

            MemoryKernel_3 = MemoryKernel_3_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            rho = zero(complex(rho0))
            for a = 1:size(rho, 1),
                b = 1:size(rho, 2),
                c = 1:size(rho, 1),
                d = 1:size(rho, 2)

                rho[a, b] += MemoryKernel_3[a, b, c, d] * rho0[c, d]
            end
            MemoryKernel_3_ref = H_II_s * W0.data * H_II_t
            rho_ref = trace_bath(
                MemoryKernel_3_ref,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))

            MemoryKernel_4 = MemoryKernel_4_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            rho = zero(complex(rho0))
            for a = 1:size(rho, 1),
                b = 1:size(rho, 2),
                c = 1:size(rho, 1),
                d = 1:size(rho, 2)

                rho[a, b] += MemoryKernel_4[a, b, c, d] * rho0[c, d]
            end
            MemoryKernel_4_ref = W0.data * H_II_s * H_II_t
            rho_ref = trace_bath(
                MemoryKernel_4_ref,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))
        end
    end

    tspan = [0.0:0.1:0.3;]
    for t_i = 1:length(tspan)
        for s_i = 1:t_i
            t = tspan[t_i]
            s = tspan[s_i]
            H_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
            H_II_s = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, s)

            MemoryKernel = MemoryKernel_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            rho = zero(complex(rho0))
            for a = 1:size(rho, 1),
                b = 1:size(rho, 2),
                c = 1:size(rho, 1),
                d = 1:size(rho, 2)

                rho[a, b] += MemoryKernel[a, b, c, d] * rho0[c, d]
            end
            MemoryKernel_ref = H_II_s * W0.data - W0.data * H_II_s
            MemoryKernel_ref = H_II_t * MemoryKernel_ref - MemoryKernel_ref * H_II_t
            rho_ref = trace_bath(
                MemoryKernel_ref,
                agg,
                FCProd,
                aggInds,
                vibindices
            )
            @test 1e-6 > D(rho, rho_ref)
        end
    end

end # testset
