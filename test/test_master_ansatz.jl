using Test
using Random, SparseArrays, LinearAlgebra, StableRNGs
using OpenQuantumSystems
import DelayDiffEq
import QuantumOpticsBase

@testset "master ansatz" begin

    Random.seed!(StableRNG(0), 1)

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
    aggTools = agg.tools
    aggOperators = agg.operators

    Ham_B = agg.operators.Ham_B
    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct

    data = Matrix(Hermitian(rand(ComplexF64, indicesLen, indicesLen)))
    W0 = DenseOperator(basis, basis, data)
    rho0 = trace_bath(W0, aggCore, aggTools)
    W0_bath = get_rho_bath(W0, aggCore, aggTools)
    normalize!(rho0)
    # tests have to be quick enough
    t_max = 0.001
    t_count = 10
    t0 = 0.
    t_step = (t_max - t0) / (t_count)
    tspan = [t0:t_step:t_max;]

    _, W_t = QME_sI_ansatz_const_test(
        W0,
        tspan,
        agg;
        reltol = 1e-3,
        abstol = 1e-3,
        int_reltol = 1e-4,
        int_abstol = 1e-4,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5()),
    )

    _, W_t = QME_sI_ansatz_const(
        W0,
        tspan,
        agg;
        reltol = 1e-3,
        abstol = 1e-3,
        int_reltol = 1e-4,
        int_abstol = 1e-4,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5()),
    )

end # testset

