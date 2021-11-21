using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs
import QuantumOpticsBase

import DelayDiffEq

@testset "master" begin

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

    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    basis = agg.tools.basis
    aggIndsLen = agg.tools.bSize

    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)
    # tests have to be quick enough
    tspan = [0.0:0.002:0.02;]

    # TODO: chenge back reltol, abstol after solving QuadGK
    # TODO: increase precision
    T, rho_t = master_int(
        rho0,
        tspan,
        Ham_0,
        Ham_I;
        reltol = 1e-6,
        abstol = 1e-6,
        int_reltol = 1e-8,
        int_abstol = 0.0,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5()),
    )#, alg=OrdinaryDiffEq.Vern7())
    rho_prev = deepcopy(rho0)
    for t_i = 2:length(rho_t)
        t = T[t_i]
        rho_I = rho_t[t_i]
        U_op_S = evolutionOperator(Ham_0, t)
        rho = U_op_S * rho_I * U_op_S'
        U_op = evolutionOperator(Ham, t)
        rho_ref = U_op * rho0 * U_op'
        @test 1e-4 > D(rho_ref, rho)
    end

    # TODO: increase precision
    T, rho_t = master(
        rho0,
        tspan,
        Ham;
        reltol = 1e-6,
        abstol = 1e-6,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5()),
    )#, alg=OrdinaryDiffEq.Vern7())
    rho_prev = deepcopy(rho0)
    for t_i = 2:length(rho_t)
        t = T[t_i]
        rho = rho_t[t_i]
        U_op = evolutionOperator(Ham, t)
        rho_ref = U_op * rho0 * U_op'
        @test 1e-4 > D(rho_ref, rho)
    end



end # testset
