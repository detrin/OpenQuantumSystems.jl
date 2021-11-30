using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

import OrdinaryDiffEq

@testset "liouville" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

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
    basis = agg.tools.basis

    ket0 = randstate(basis)
    rho0 = dm(ket0)
    # tests have to be quick enough
    tspan = [0.0:0.1:1.0;]

    T, rho_t = LvN_SS(
        rho0,
        tspan,
        agg;
        reltol = 1e-10,
        abstol = 1e-10,
        alg = OrdinaryDiffEq.Tsit5(),
    )
    for t_i = 1:length(tspan)
        U_op = evolutionOperator(Ham, tspan[t_i])
        rho = U_op * rho0 * U_op'
        @test 1e-7 > D(rho, rho_t[t_i])
        # println(t_i, " ", D(rho.data, rho_t[t_i].data))
    end

    # diabled for now for, it takes about 20 s to run this simulation
    #=
    T, rho_t = liouvilleVonNeumann(rho0, Ham, tspan; reltol=1e-12, abstol=1e-12, alg=OrdinaryDiffEq.Vern7())
    for t_i in 1:length(tspan)
        U_op = evolutionOperator(Ham, tspan[t_i])
        rho = U_op * rho0 * U_op'
        @test 1e-14 > D(rho, rho_t[t_i])
        # println(t_i, " ", D(rho.data, rho_t[t_i].data))
    end
    =#


end # testset
