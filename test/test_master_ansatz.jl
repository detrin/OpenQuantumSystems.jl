using Test
using Random, SparseArrays, LinearAlgebra
using OpenQuantumSystems
import DelayDiffEq
import QuantumOpticsBase

@testset "master ansatz" begin

    Random.seed!(0)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    HR = 0.01
    shift = (2.0 * HR)^0.5
    mode1 = Mode(300., shift)
    mol1 = Molecule([mode1], 3, [12500., 12750.])
    mol2 = Molecule([mode1], 3, [12500., 12800.])
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState = true)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState = true)
    FCProd = getFCProd(agg, FCFact, aggInds, vibindices; groundState = true)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState = true)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundState = false)
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggInds, FCFact; groundState = true)
    Ham_S = Ham - Ham_int

    H_lambda, H_S = eigen(Ham_S.data)
    H_Sinv = inv(H_S)
    H_lambda = diagm(H_lambda)

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
    rho_traced1 = trace_bath(W01.data, agg, FCProd, aggInds, vibindices; groundState=false)
    mu_array = [[1, 2]]
    W02 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize=false, diagonal=true)
    rho_traced2 = trace_bath(W02.data, agg, FCProd, aggInds, vibindices; groundState=false)

    W0 = 1.0*W01 + 0.0*W02 
    W0_bath = get_rho_bath(W0, agg, FCProd, aggInds, vibindices; groundState=true)
    rho0 = trace_bath(W0, agg, FCProd, aggInds, vibindices; groundState=false)

    W0 = DenseOperator(W0.basis_l, W0.basis_r, complex(W0.data))
    W0_bath = DenseOperator(W0_bath.basis_l, W0_bath.basis_r, complex(W0_bath.data))
    rho0 = DenseOperator(rho0.basis_l, rho0.basis_r, complex(rho0.data))

    p = (Ham_S, Ham_int, H_lambda, H_S, H_Sinv, W0, W0_bath, agg, FCProd, aggInds, vibindices, false, ComplexF64)
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

