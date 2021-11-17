using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra

@testset "trace" begin

    Random.seed!(0)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    mode1 = Mode(0.2, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 3, Energy)
    mol2 = Molecule([mode1], 3, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds)
    Ham = getAggHamiltonian(agg, aggInds, FCFact)
    Ham_S = getAggHamSysBath2(agg, aggInds)
    Ham_int = Ham - Ham_S

    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.38039889598242 + 0.0im 0.28369491536092306 + 0.248783725284891im 0.2196786799502027 + 0.25348757179772047im; 
        0.28369491536092306 - 0.248783725284891im 0.2741822201009351 - 1.3010426069826053e-17im 0.11525284205621783 + 0.23102806588760527im; 
        0.21967867995020277 - 0.2534875717977205im 0.11525284205621782 - 0.23102806588760527im 0.3377238162369547 - 1.3769367590565906e-17im
    ]
    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices)
    println(rho_traced.data)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, agg, FCProd, aggInds, vibindices)
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho0_ad = ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices)
    @test 1e-14 > D(rho0, rho0_ad)

    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.3803988959824186 + 2.6020852139652106e-18im 0.3487316508666893 - 0.03548586492564312im 0.3626658702848767 + 0.0005197749818341463im; 
        0.3487316508666892 + 0.03548586492564312im 0.2560172766886787 - 1.6543652845120398e-17im 0.11949445051795413 + 0.2239460097353317im; 
        0.3626658702848766 - 0.0005197749818341714im 0.11949445051795418 - 0.2239460097353317im 0.33235575587167504 - 4.0115480381963664e-18im
    ]

    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices)
    println(rho_traced.data)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, agg, FCProd, aggInds, vibindices)
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)

        rho_ab = take_el_part(rho.data, a, b, vibindices)
        rho_traced_ab = trace_bath_part(
            rho_ab,
            a,
            b,
            agg,
            FCProd,
            aggInds,
            vibindices
        )
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho0_ad = ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices)
    @test 1e-14 > D(rho0, rho0_ad)

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices)

    t = 0.0
    Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    prod = Ham_II_t * Ham_int * rho_bath
    corr_ref = trace_bath(prod, agg, FCProd, aggInds, vibindices)
    corr = correlation_function(
        t,
        rho_bath,
        Ham_S,
        Ham_int,
        agg,
        FCProd,
        aggInds,
        vibindices
    )
    @test 1e-14 > D(corr, corr_ref)

    t = 1.0
    Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    prod = Ham_II_t * Ham_int * rho_bath
    corr_ref = trace_bath(prod, agg, FCProd, aggInds, vibindices)
    corr = correlation_function(
        t,
        rho_bath,
        Ham_S,
        Ham_int,
        agg,
        FCProd,
        aggInds,
        vibindices
    )
    @test 1e-14 > D(corr, corr_ref)

end
