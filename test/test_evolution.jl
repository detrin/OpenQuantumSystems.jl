using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra

@testset "evolution" begin

    Random.seed!(0)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))

    mode1 = Mode(0.2, 1.)
    Energy = [0., 200.]
    mol1 = Molecule([mode1], 2, Energy)
    mol2 = Molecule([mode1], 2, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState=false)
    base = GenericBasis([length(aggInds)])
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

    U_op_array = EvolutionOperatorArray(Ham, 0.0, 1.0, 3)
    U_op1 = EvolutionOperator(Ham, 0.0)
    U_op2 = EvolutionOperator(Ham, 0.5)
    U_op3 = EvolutionOperator(Ham, 1.0)
    @test 1e-12 > D(U_op_array[1], U_op1)
    @test 1e-12 > D(U_op_array[2], U_op2)
    @test 1e-12 > D(U_op_array[3], U_op3)

    U_sop_array = EvolutionSuperOperatorArray(Ham, 0.0, 1.0, 3)
    U_sop1 = EvolutionSuperOperator(Ham, 0.0)
    U_sop2 = EvolutionSuperOperator(Ham, 0.5)
    U_sop3 = EvolutionSuperOperator(Ham, 1.0)
    @test 1e-12 > D(U_sop_array[1], U_sop1)
    @test 1e-12 > D(U_sop_array[2], U_sop2)
    @test 1e-12 > D(U_sop_array[3], U_sop3)

    t = 0.
    foreach(EvolutionOperatorIterator(Ham, 0.0, 1.0, 3)) do U_op
        U_op_ref = EvolutionOperator(Ham, t)
        @test 1e-12 > D(U_op, U_op_ref)
        t += 0.5
    end

    t = 0.
    foreach(EvolutionSuperOperatorIterator(Ham, 0.0, 1.0, 3)) do U_sop
        U_sop_ref = EvolutionSuperOperator(Ham, t)
        @test 1e-12 > D(U_sop, U_sop_ref)
        t += 0.5
    end

    ket = randstate(base)
    ket1 = U_op1 * ket
    ket2 = U_op2 * ket
    ket3 = U_op3 * ket
    ket_array = EvolutionExact(ket, Ham, 0.0, 1.0, 3)
    @test 1e-12 > D(ket1, ket_array[1])
    @test 1e-12 > D(ket2, ket_array[2])
    @test 1e-12 > D(ket3, ket_array[3])

    ket_array = EvolutionApproximate(ket, Ham, 0.0, 1.0, 3)
    @test 1e-12 > D(ket1, ket_array[1])
    @test 1e-12 > D(ket2, ket_array[2])
    @test 1e-12 > D(ket3, ket_array[3])

    op = dm(ket)
    op1 = U_sop1 * op
    op2 = U_sop2 * op
    op3 = U_sop3 * op
    op_array = EvolutionExact(op, Ham, 0.0, 1.0, 3)
    @test 1e-12 > D(op1, op_array[1])
    @test 1e-12 > D(op2, op_array[2])
    @test 1e-12 > D(op3, op_array[3])

    op_array = EvolutionApproximate(op, Ham, 0.0, 1.0, 3)
    @test 1e-12 > D(op1, op_array[1])
    @test 1e-12 > D(op2, op_array[2])
    @test 1e-12 > D(op3, op_array[3])
    

end