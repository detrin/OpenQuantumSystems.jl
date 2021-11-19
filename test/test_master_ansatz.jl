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
    HR = 0.01
    shift = (2.0 * HR)^0.5
    mode1 = Mode(300., shift)
    mol1 = Molecule([mode1], 3, [12500., 12750.])
    mol2 = Molecule([mode1], 3, [12500., 12800.])
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds)
    FCProd = getFCProd(agg, FCFact, aggInds, vibindices)

    Ham_sys = getAggHamSystemSmall(agg)
    Ham_I = getAggHamInteraction(agg, aggInds, FCFact)
    Ham_0 = getAggHamSystemBath(agg, aggInds, FCFact)
    Ham_B = getAggHamBathBig(agg)
    Ham = getAggHamiltonian(agg, aggInds, FCFact)

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
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
    W01 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize=false, diagonal=true)
    rho_traced1 = trace_bath(W01.data, agg, FCProd, aggInds, vibindices)
    mu_array = [[1, 2]]
    W02 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize=false, diagonal=true)
    rho_traced2 = trace_bath(W02.data, agg, FCProd, aggInds, vibindices)

    W0 = 1.0*W01 + 0.0*W02 
    W0_bath = get_rho_bath(W0, agg, FCProd, aggInds, vibindices)
    rho0 = trace_bath(W0, agg, FCProd, aggInds, vibindices)

    W0 = DenseOperator(W0.basis_l, W0.basis_r, complex(W0.data))
    W0_bath = DenseOperator(W0_bath.basis_l, W0_bath.basis_r, complex(W0_bath.data))
    rho0 = DenseOperator(rho0.basis_l, rho0.basis_r, complex(rho0.data))

    p = (Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, Ham_B, W0, W0_bath, agg, FCProd, aggInds, vibindices, ComplexF64)
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

