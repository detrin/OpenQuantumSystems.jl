using Test
using OpenQuantumSystems
using SparseArrays, LinearAlgebra

@testset "interaction picture" begin

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
    mol1 = Molecule([mode1], 3, [2.0, 200.0])
    mol2 = Molecule([mode2], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)

    Ham = agg.operators.Ham

    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    for t in [0.0, 1.0, 2.0]
        U_op = evolutionOperator(Ham_0, t)
        U_II_t_ref = U_op' * Ham_I * U_op

        U_II_t = getInteractionHamIPicture(Ham_0.data, Ham_I.data, t)
        @test 1e-14 > D(U_II_t_ref.data, U_II_t)
        U_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
        @test 1e-14 > D(U_II_t_ref, U_II_t)
        U_II_t = getInteractionHamIPictureA(Ham_0, Ham_I, t)
        @test 1e-14 > D(U_II_t_ref.data, U_II_t)
        U_II_t = getInteractionHamIPictureA(Ham_I.data, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, t)
        @test 1e-10 > D(U_II_t_ref.data, U_II_t)
    end

end # testset
