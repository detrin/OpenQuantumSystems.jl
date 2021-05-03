using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra

@testset "liouville" begin

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

    ket0 = randstate(base)
    rho0 = dm(ket0)
    tspan = [0.:0.1:1.0;]
    
    U_sop_array = evolutionSuperOperatorArray(Ham, tspan)
    op_array = evolutionExact(rho0, Ham, tspan) #liouvilleVonNeumann
    for t_i in 1:length(tspan)
        rho = U_sop_array[t_i] * rho0
        @test 1e-12 > D(rho, op_array[t_i])
    end

end # testset