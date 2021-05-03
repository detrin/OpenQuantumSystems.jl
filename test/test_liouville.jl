using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra

import OrdinaryDiffEq

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
    aggIndsLen = length(aggInds)
    base = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState=false)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState=false)

    ket0 = randstate(base)
    rho0 = dm(ket0)
    # tests have to be quick enough
    tspan = [0.:0.1:1.0;]

    T, rho_t = liouvilleVonNeumann(rho0, Ham, tspan; reltol=1e-10, abstol=1e-10, alg=OrdinaryDiffEq.Tsit5())
    for t_i in 1:length(tspan)
        U_op = evolutionOperator(Ham, tspan[t_i])
        rho = U_op * rho0 * U_op'
        @test 1e-10 > D(rho, rho_t[t_i])
        # println(t_i, " ", D(rho.data, rho_t[t_i].data))
    end

    T, rho_t = liouvilleVonNeumann(rho0, Ham, tspan; reltol=1e-12, abstol=1e-12, alg=OrdinaryDiffEq.Vern7())
    for t_i in 1:length(tspan)
        U_op = evolutionOperator(Ham, tspan[t_i])
        rho = U_op * rho0 * U_op'
        @test 1e-15 > D(rho, rho_t[t_i])
        # println(t_i, " ", D(rho.data, rho_t[t_i].data))
    end
    

end # testset