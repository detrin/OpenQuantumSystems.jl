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
    mol1 = Molecule([mode1], 2, [2.0, 200.0])
    mol2 = Molecule([mode2], 2, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)
    aggTools = agg.tools
    aggOperators = agg.operators

    Ham_sys = agg.operators.Ham_sys 
    Ham_S = agg.operators.Ham_S 
    Ham_B = agg.operators.Ham_B 
    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    data = Matrix(Hermitian(rand(ComplexF64, indicesLen, indicesLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)
    # tests have to be quick enough
    t_max = 0.001
    t_count = 10
    t0 = 0.
    t_step = (t_max - t0) / (t_count)
    tspan = [t0:t_step:t_max;]

    T = 300
    mu_array = [[2, 1]]
    W01 = thermal_state(T, mu_array, aggCore, aggTools, aggOperators; diagonalize=false, diagonal=true)
    rho_traced1 = trace_bath(W01.data, aggCore, aggTools)
    mu_array = [[1, 2]]
    W02 = thermal_state(T, mu_array, aggCore, aggTools, aggOperators; diagonalize=false, diagonal=true)
    rho_traced2 = trace_bath(W02.data, aggCore, aggTools)

    W0 = 1.0*W01 + 0.0*W02 
    W0_bath = get_rho_bath(W0, aggCore, aggTools)
    rho0 = trace_bath(W0, aggCore, aggTools)

    W0 = DenseOperator(W0.basis_l, W0.basis_r, complex(W0.data))
    W0_bath = DenseOperator(W0_bath.basis_l, W0_bath.basis_r, complex(W0_bath.data))
    rho0 = DenseOperator(rho0.basis_l, rho0.basis_r, complex(rho0.data))

    p = (Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, Ham_B, W0, W0_bath, agg, FCProd, indices, indicesMap, ComplexF64, aggCore, aggTools)
    Tspan, rho_t = master_ansatz(
        rho0,
        tspan,
        p;
        reltol = 1.0e-4,
        abstol = 1.0e-4,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5())
    )
    rho_prev = deepcopy(rho0)
    for t_i = 2:length(rho_t)
        t = Tspan[t_i]
        rho_I = rho_t[t_i]
        U_op_S = evolutionOperator(Ham_sys, t)
        rho = U_op_S * rho_I * U_op_S'
    end

end # testset

