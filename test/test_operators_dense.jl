using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "operators dense" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::Array, op2::Array) = abs(norm(op1 - op2))

    b = GenericBasis([3])
    one_op = OneDenseOperator(b)
    one_op_ref = [
        1.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im 1.0+0.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im 1.0+0.0im
    ]
    @test 1e-12 > D(one_op.data, one_op_ref)

end # testset
