using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "trace" begin

    Random.seed!(StableRNG(0), 1)

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

    data = Matrix(Hermitian(rand(StableRNG(0), ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[0.30545089951602383 + 0.0im 0.14375197756961597 + 0.21658981002618588im 0.23626964940425038 + 0.1914303091850988im; 0.14375197756961597 - 0.21658981002618588im 0.30379617012935617 + 1.5178830414797062e-18im 0.15887255193342634 + 0.11195793936579877im; 0.23626964940425035 - 0.19143030918509882im 0.15887255193342636 - 0.11195793936579874im 0.4026345547478164 - 1.734723475976807e-18im]
    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, agg, FCProd, aggInds, vibindices)
        @test 1e-7 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho0_ad = ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices)
    @test 1e-7 > D(rho0, rho0_ad)

    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[0.30545089951602455 + 0.0im 0.27487510559103273 + 0.03248766557875096im 0.30184368281947643 - 0.05602820966070799im; 0.2748751055910327 - 0.03248766557875095im 0.2979302015891273 + 3.0612122231274e-19im 0.17734441884490892 + 0.11901706815341803im; 0.3018436828194764 + 0.05602820966070803im 0.177344418844909 - 0.11901706815341806im 0.3955255687091211 + 1.0315865594447475e-17im]
    rho_traced = trace_bath_slow(rho, agg, FCFact, aggInds, vibindices)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced =
        trace_bath_slow(rho.data, agg, FCFact, aggInds, vibindices)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    FCProd = getFCProd(agg, FCFact, aggInds, vibindices)
    rho_traced = trace_bath(rho, agg, FCProd, aggInds, vibindices)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCProd, aggInds, vibindices)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, agg, FCProd, aggInds, vibindices)
        @test 1e-7 > abs(rho_traced_ref[a, b] - rho_traced_ab)

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
        @test 1e-7 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end

    rho_bath = get_rho_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho_traced = trace_bath(rho_bath, agg, FCProd, aggInds, vibindices)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, agg, FCProd, aggInds, vibindices)
    rho0_ad = ad(rho_traced, rho_bath, agg, FCProd, aggInds, vibindices)
    @test 1e-7 > D(rho0, rho0_ad)

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
    @test 1e-7 > D(corr, corr_ref)

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
    @test 1e-7 > D(corr, corr_ref)

end
