using Test
using OpenQuantumSystems
using SparseArrays, LinearAlgebra

@testset "superoperators" begin

    b = GenericBasis([2])
    A = DenseOperator(b, b, [1 2; 3 4])
    B = DenseOperator(b, b, [5 6; 7 8])
    C = Commutator(A)
    C_ref = DenseSuperOperator(
        (b, b),
        (b, b),
        [0 2 -3 0; 3 3 0 -3; -2 0 -3 2; 0 -2 3 0],
    )
    @test Commutator(A) == C_ref
    @test Commutator(A) * B == A * B - B * A

end # testset
