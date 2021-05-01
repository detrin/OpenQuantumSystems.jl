using Test
using OpenQuantumSystems
using LinearAlgebra
using SparseArrays

@testset "evolution" begin

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(op1::AbstractOperator, op2::AbstractOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))

    mode1 = Mode(0.2, 1.)
    Energy = [0., 200.]
    mol1 = Molecule([mode1], 2, Energy)
    mol2 = Molecule([mode1], 2, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState=false)
    FCFact = getFranckCondonFactors(agg, aggInds; groundState=false)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState=false)

    t = 0.
    U_op = EvolutionOperator(Ham, t)
    U_op_ref = exp(-1im * Ham * t)
    @test 1e-12 > D(U_op, U_op_ref)

    U_sop = EvolutionSuperOperator(Ham, t)
    U_sop_ref = spre(U_op') * spost(U_op)
    @test 1e-12 > D(U_sop, U_sop_ref)

    t = 1.
    U_op = EvolutionOperator(Ham, t)
    U_op_ref = exp(-1im * Ham * t)
    @test 1e-12 > D(U_op, U_op_ref)

    U_sop = EvolutionSuperOperator(Ham, t)
    U_sop_ref = spre(U_op') * spost(U_op)
    @test 1e-12 > D(U_sop, U_sop_ref)

end