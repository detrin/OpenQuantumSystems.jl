using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs
import QuantumOpticsBase


@testset "memory kernel" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    # TODO: change to macro
    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 2, [2.0, 200.0])
    mol2 = Molecule([mode2], 2, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)
    aggTools = agg.tools
    aggOperators = agg.operators

    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    # T = 300
    # W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, indices; diagonalize=false, diagonal=true)
    data = Matrix(Hermitian(rand(ComplexF64, indicesLen, indicesLen)))
    W0 = DenseOperator(basis, basis, data)
    normalize!(W0)
    W0_bath = get_rho_bath(
        W0,
        aggCore,
        aggTools;
        justCopy = true,
    )
    rho0 = trace_bath(W0.data, aggCore, aggTools)

    W_ab = take_el_part(W0.data, 1, 1, indicesMap)
    W_ab_len = size(W_ab, 1)
    @test 1e-15 > D(W_ab, W0.data[1:W_ab_len, 1:W_ab_len])

    tspan = [0.0:0.1:0.3;]
    for t_i = 1:length(tspan)
        for s_i = 1:t_i
            t = tspan[t_i]
            s = tspan[s_i]
            # println(t, " ", s)
            H_II_t = getInteractionHamIPictureA(Ham_I.data, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, t)
            H_II_s = getInteractionHamIPictureA(Ham_I.data, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, s)

            MemoryKernel_1 = MemoryKernel_1_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                aggCore,
                aggTools
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
                aggCore, 
                aggTools
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))

            MemoryKernel_2 = MemoryKernel_2_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                aggCore,
                aggTools
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
                aggCore, 
                aggTools
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))

            MemoryKernel_3 = MemoryKernel_3_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                aggCore,
                aggTools
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
                aggCore, 
                aggTools
            )
            @test 1e-6 > D(rho, rho_ref)
            # println(D(rho, rho_ref))

            MemoryKernel_4 = MemoryKernel_4_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                aggCore,
                aggTools
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
                aggCore, 
                aggTools
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
            H_II_t = getInteractionHamIPictureA(Ham_I.data, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, t)
            H_II_s = getInteractionHamIPictureA(Ham_I.data, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, s)

            MemoryKernel = MemoryKernel_traced(
                H_II_t,
                H_II_s,
                W0_bath.data,
                aggCore,
                aggTools
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
                aggCore, 
                aggTools
            )
            @test 1e-6 > D(rho, rho_ref)
        end
    end

end # testset
