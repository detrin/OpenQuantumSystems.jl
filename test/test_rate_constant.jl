using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "rate_constant" begin

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

    tspan = get_tspan(0., 0.02, 100)
    W0, rho0, W0_bath = ultrafast_laser_excitation(10., [0.0, 0.3, 0.7], agg)

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    t = 1.
    # TODO: increase accuracy
    ref = ComplexF64[-0.0 - 0.0im -0.0 - 0.0im -0.0 - 0.0im; -0.0 - 0.0im -1650.7572853212348 - 3.7298879616771476e-14im 3167.3936977428143 + 6.824421196019823e-14im; -0.0 - 0.0im 1082.1218357464745 + 1.7493984984155607e-14im -1979.427319937707 - 3.700622271661764e-14im]
    K_ab_ = K_ab(t, p, tmp1, tmp2)
    @test D(ref, K_ab_) < 1e-10
    
    K_ab_ = K_ab(t, W0, W0_bath, agg)
    @test D(ref, K_ab_) < 1e-10

end