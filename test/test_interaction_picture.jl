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

    mode1 = Mode(0.2, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 2, Energy)
    mol2 = Molecule([mode1], 2, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds)
    Ham = getAggHamiltonian(agg, aggInds, FCFact)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg)
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggInds, FCFact)
    Ham_S = Ham - Ham_int

    H_lambda, H_S = eigen(Ham_S.data)
    H_Sinv = inv(H_S)
    H_lambda = diagm(H_lambda)

    t = 0.0
    U_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    U_op = evolutionOperator(Ham_S, t)
    U_II_t_ref = U_op' * Ham_int * U_op
    @test 1e-15 > D(U_II_t_ref, U_II_t)
    U_II_t = getInteractionHamIPictureA(Ham_S, Ham_int, t)
    @test 1e-15 > D(U_II_t_ref.data, U_II_t)
    U_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
    @test 1e-15 > D(U_II_t_ref.data, U_II_t)

    t = 1.0
    U_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    U_op = evolutionOperator(Ham_S, t)
    U_II_t_ref = U_op' * Ham_int * U_op
    @test 1e-15 > D(U_II_t_ref, U_II_t)
    U_II_t = getInteractionHamIPictureA(Ham_S, Ham_int, t)
    @test 1e-15 > D(U_II_t_ref.data, U_II_t)
    U_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
    @test 1e-15 > D(U_II_t_ref.data, U_II_t)

    t = 2.0
    U_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    U_op = evolutionOperator(Ham_S, t)
    U_II_t_ref = U_op' * Ham_int * U_op
    @test 1e-15 > D(U_II_t_ref, U_II_t)
    U_II_t = getInteractionHamIPictureA(Ham_S, Ham_int, t)
    @test 1e-15 > D(U_II_t_ref.data, U_II_t)
    U_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
    @test 1e-15 > D(U_II_t_ref.data, U_II_t)


end # testset
