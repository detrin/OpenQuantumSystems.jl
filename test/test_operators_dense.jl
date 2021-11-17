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

    annihilation_op = AnnihilationOperator(b)
    annihilation_op_ref = [0.0 1.0 0.0; 0.0 0.0 1.4142135623730951; 0.0 0.0 0.0]
    @test 1e-12 > D(annihilation_op.data, annihilation_op_ref)

    creation_op = CreationOperator(b)
    creation_op_ref = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.4142135623730951 0.0]
    @test 1e-12 > D(creation_op.data, creation_op_ref)

    position_op = PositionOperator(b)
    position_op_ref = [0.0 0.7071067811865475 0.0; 0.7071067811865475 0.0 1.0; 0.0 1.0 0.0]
    @test 1e-12 > D(position_op.data, position_op_ref)

    momentum_op = MomentumOperator(b)
    momentum_op_ref = ComplexF64[0.0 + 0.0im -0.0 - 0.7071067811865475im 0.0 + 0.0im; 0.0 + 0.7071067811865475im 0.0 + 0.0im -0.0 - 1.0im; 0.0 + 0.0im 0.0 + 1.0im 0.0 + 0.0im]
    @test 1e-12 > D(momentum_op.data, momentum_op_ref)

    shift_op = ShiftOperator(b, 1.0)
    shift_op_ref = [0.779728662995649 + 0.0im -0.5431245605675081 + 0.0im 0.3115107121936076 + 0.0im; 0.5431245605675084 + 0.0im 0.33918598898694735 + 0.0im -0.7680941196124976 + 0.0im; 0.31151071219360765 + 0.0im 0.7680941196124975 + 0.0im 0.559457325991298 + 0.0im]
    @test 1e-12 > D(shift_op.data, shift_op_ref)


end # testset
