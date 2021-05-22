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
    aggInds = getIndices(agg; groundState = false)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState = false)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState = false)
    Ham_S = getAggHamSysBath2(agg, aggInds; groundState = false)
    Ham_int = Ham - Ham_S

    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.5373667261645646+0.0im 0.5960193120071318+0.31849769746151213im
        0.5960193120071318-0.31849769746151213im 0.4626332738354356+0.0im
    ]
    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices; groundState = false)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    for a = 1:2, b = 1:2
        rho_traced_ab = trace_bath(rho, a + 1, b + 1, agg, FCProd, aggInds, vibindices)
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = false)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices; groundState = false)
    rho_traced_ref = [1.0 1.0; 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = false)
    rho0_ad =
        ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho0, rho0_ad)

    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.5373667261645646+4.804317299001507e-18im 0.5960193120071318+0.31849769746151213im
        0.5960193120071318-0.31849769746151213im 0.46263327383543557-4.716400113933126e-18im
    ]
    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices; groundState = false)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    for a = 1:2, b = 1:2
        rho_traced_ab = trace_bath(rho, a + 1, b + 1, agg, FCProd, aggInds, vibindices)
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = false)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices; groundState = false)
    rho_traced_ref = [1.0 1.0; 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = false)
    rho0_ad =
        ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices; groundState = false)
    @test 1e-14 > D(rho0, rho0_ad)

    mode1 = Mode(0.2, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 3, Energy)
    mol2 = Molecule([mode1], 3, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState = true)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState = true)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState = true)
    Ham_S = getAggHamSysBath2(agg, aggInds; groundState = true)
    Ham_int = Ham - Ham_S

    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.27983831453187563+0.0im 0.39510751997977556+0.3832496461602654im 0.36717929072360034+0.48380689550795175im
        0.39510751997977556-0.3832496461602654im 0.3757521692494342+0.0im 0.4295074106269881+0.3124835865738737im
        0.36717929072360034-0.48380689550795175im 0.4295074106269881-0.3124835865738737im 0.3444095162186903+0.0im
    ]
    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices; groundState = true)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, agg, FCProd, aggInds, vibindices)
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = true)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices; groundState = true)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = true)
    rho0_ad = ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho0, rho0_ad)

    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.27983831453187474+4.807665611244999e-18im 0.5271823941867131-0.15833182524717854im 0.6013924791789407-0.08495192560698843im
        0.5271823941867131+0.15833182524717854im 0.3757521692494339+5.194043729796511e-18im 0.4295074106269877+0.31248358657387476im
        0.6013924791789407+0.08495192560698844im 0.4295074106269877-0.3124835865738748im 0.34440951621869137+4.376703186460492e-18im
    ]

    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices; groundState = true)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices; groundState = true)
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
            vibindices;
            groundState = true,
        )
        @test 1e-14 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = true)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices; groundState = true)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = true)
    rho0_ad = ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices; groundState = true)
    @test 1e-14 > D(rho0, rho0_ad)

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices; groundState = true)

    t = 0.0
    Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    prod = Ham_II_t * Ham_int * rho_bath
    corr_ref = trace_bath(prod, agg, FCProd, aggInds, vibindices; groundState = true)
    corr = correlation_function(
        t,
        rho_bath,
        Ham_S,
        Ham_int,
        agg,
        FCProd,
        aggInds,
        vibindices;
        groundState = true,
    )
    @test 1e-14 > D(corr, corr_ref)

    t = 1.0
    Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    prod = Ham_II_t * Ham_int * rho_bath
    corr_ref = trace_bath(prod, agg, FCProd, aggInds, vibindices; groundState = true)
    corr = correlation_function(
        t,
        rho_bath,
        Ham_S,
        Ham_int,
        agg,
        FCProd,
        aggInds,
        vibindices;
        groundState = true,
    )
    @test 1e-14 > D(corr, corr_ref)

end
