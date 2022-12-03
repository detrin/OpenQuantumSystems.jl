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
    aggOperators = agg.operators
    aggTools = agg.tools

    # not dependent on aggregate
    Ham_I = Matrix(rand(StableRNG(1), Float64, agg.tools.bSize, agg.tools.bSize))
    Ham_I = DenseOperator(agg.tools.basis, agg.tools.basis, Ham_I)
    Ham_0 = Matrix(rand(StableRNG(2), Float64, agg.tools.bSize, agg.tools.bSize))
    Ham_0 = DenseOperator(agg.tools.basis, agg.tools.basis, Ham_0)
    Ham = Ham_I + Ham_0

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct


    data = Matrix(Hermitian(rand(StableRNG(0), ComplexF64, indicesLen, indicesLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[0.30545089951602383 + 0.0im 0.14375197756961597 + 0.21658981002618588im 0.09446278210921694 + 0.06885750587885092im; 0.14375197756961597 - 0.21658981002618588im 0.30379617012935617 + 1.5178830414797062e-18im 0.06403549575681955 + 0.0490347030584035im; 0.09446278210921695 - 0.06885750587885092im 0.06403549575681954 - 0.04903470305840349im 0.3601657072334817 - 1.1629184732212025e-18im]
    rho_traced = trace_bath_slow(rho, aggCore, aggTools)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath_slow(rho.data, aggCore, aggTools)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    rho_traced_ref = ComplexF64[0.30545089951602383 + 0.0im 0.2935867848856678 + 0.3293162968124289im 0.31391118201089463 + 0.26706491332993326im; 0.2935867848856678 - 0.3293162968124289im 0.31736579991967256 + 0.0im 0.33780358907521824 + 0.297414426438633im; 0.31391118201089463 - 0.26706491332993326im 0.33780358907521824 - 0.297414426438633im 0.37718330056430377 + 0.0im]
    rho_traced = trace_bath(rho, aggCore, aggOperators, aggTools)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, aggCore, aggOperators, aggTools)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, aggTools; vib_basis=aggOperators.vib_basis)
        @test 1e-7 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end
    
    rho_bath = get_rho_bath(rho0, aggCore, aggOperators, aggTools)
    rho_traced = trace_bath(rho_bath, aggCore, aggOperators, aggTools)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, aggCore, aggOperators, aggTools)
    rho0_ad = ad(rho_traced, rho_bath, aggCore, aggTools)
    @test 1e-7 > D(rho0, rho0_ad)
    
    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[0.2595039421087527 + 1.1102230246251565e-16im 0.6626972429778065 - 0.3034395850675509im 0.12521539202671672 - 0.4920188624519769im; 0.6626972429778069 + 0.3034395850675506im 0.21095069278599016 + 2.2730010622775303e-16im 0.29117337783386915 + 1.2392616289083571im; 0.12521539202671667 + 0.492018862451977im 0.29117337783386926 - 1.239261628908357im 1.3575702817649218 + 4.0047474713266894e-17im]
    rho_traced = trace_bath_slow(rho, aggCore, aggTools)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath_slow(rho.data, aggCore, aggTools)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    rho_traced_ref = ComplexF64[0.25950394210875094 - 1.8422648233666214e-16im 1.018659546379254 + 0.6588686562184847im -0.40327861040729096 + 1.0419870590205613im; 1.0186595463792536 - 0.6588686562184841im -0.2813033800012592 + 7.151701628581004e-17im 0.5154723483602309 - 0.6806166222460325im; -0.4032786104072908 - 1.0419870590205615im 0.515472348360231 + 0.6806166222460321im 1.5605546808339001 + 1.7829708069959374e-16im]
    rho_traced = trace_bath(rho, aggCore, aggOperators, aggTools)
    @test 1e-7 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, aggCore, aggOperators, aggTools)
    @test 1e-7 > D(rho_traced, rho_traced_ref)

    # TODO: fix for ground_excited
    #=
    for a = 1:3, b = 1:3
        rho_traced_ab = trace_bath(rho, a, b, aggTools; vib_basis=:ground_ground)
        @test 1e-7 > abs(rho_traced_ref[a, b] - rho_traced_ab)

        rho_ab = take_el_part(rho.data, a, b, indicesMap)
        rho_traced_ab = trace_bath_part(
            rho_ab,
            a,
            b,
            aggTools;
            vib_basis=:ground_ground
        )
        @test 1e-7 > abs(rho_traced_ref[a, b] - rho_traced_ab)
    end
    =#

    rho_bath = get_rho_bath(rho0, aggCore, aggOperators, aggTools)
    rho_traced = trace_bath(rho_bath, aggCore, aggOperators, aggTools)
    rho_traced_ref = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]
    @test 1e-10 > D(rho_traced_ref, rho_traced.data)

    rho_traced = trace_bath(rho0, aggCore, aggOperators, aggTools)
    rho0_ad = ad(rho_traced, rho_bath, aggCore, aggTools)
    @test 1e-7 > D(rho0, rho0_ad)

    rho_bath = get_rho_bath(rho0, aggCore, aggOperators, aggTools)
    
    t = 0.0
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    prod = Ham_II_t * Ham_I * rho_bath
    corr_ref = trace_bath(prod, aggCore, aggOperators, aggTools)
    corr = correlation_function(
        t,
        rho_bath,
        Ham_0,
        Ham_I,
        aggCore, 
        aggOperators,
        aggTools
    )
    @test 1e-7 > D(corr, corr_ref)

    t = 1.0
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)
    prod = Ham_II_t * Ham_I * rho_bath
    corr_ref = trace_bath(prod, aggCore, aggOperators, aggTools)
    corr = correlation_function(
        t,
        rho_bath,
        Ham_0,
        Ham_I,
        aggCore, 
        aggOperators,
        aggTools
    )
    @test 1e-7 > D(corr, corr_ref)
    
end
