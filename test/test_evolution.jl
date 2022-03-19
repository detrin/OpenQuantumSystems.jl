using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "evolution" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    tspan = get_tspan(0., 1., 10)
    ref = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    @test tspan == ref

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 3, [2.0, 200.0])
    mol2 = Molecule([mode2], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)

    Ham = agg.operators.Ham
    H_lambda, H_S = eigen(Ham.data)
    H_Sinv = inv(H_S)
    H_lambda = diagm(H_lambda)

    t = 0.0
    U_op_ref = exp(-1im * Ham * t)
    U_op = evolutionOperator(Ham, t)
    @test 1e-12 > D(U_op, U_op_ref)
    U_op = evolutionOperatorA(H_lambda, H_S, H_Sinv, t)
    @test 1e-12 > D(U_op, U_op_ref.data)

    U_sop = evolutionSuperOperator(Ham, t)
    U_sop_ref = spre(U_op_ref) * spost(U_op_ref')
    @test 1e-12 > D(U_sop, U_sop_ref)

    t = 1.0
    U_op_ref = exp(-1im * Ham * t)
    U_op = evolutionOperator(Ham, t)
    @test 1e-12 > D(U_op, U_op_ref)
    U_op = evolutionOperatorA(H_lambda, H_S, H_Sinv, t)
    @test 1e-12 > D(U_op, U_op_ref.data)

    U_sop = evolutionSuperOperator(Ham, t)
    U_sop_ref = spre(U_op_ref) * spost(U_op_ref')
    @test 1e-12 > D(U_sop, U_sop_ref)

    tspan = [0.0:0.5:1.0;]
    U_op_array = evolutionOperatorArray(Ham, tspan)
    U_op1 = evolutionOperator(Ham, 0.0)
    U_op2 = evolutionOperator(Ham, 0.5)
    U_op3 = evolutionOperator(Ham, 1.0)
    @test 1e-12 > D(U_op_array[1], U_op1)
    @test 1e-12 > D(U_op_array[2], U_op2)
    @test 1e-11 > D(U_op_array[3], U_op3)

    U_sop_array = evolutionSuperOperatorArray(Ham, tspan)
    U_sop1 = evolutionSuperOperator(Ham, 0.0)
    U_sop2 = evolutionSuperOperator(Ham, 0.5)
    U_sop3 = evolutionSuperOperator(Ham, 1.0)
    @test 1e-11 > D(U_sop_array[1], U_sop1)
    @test 1e-10 > D(U_sop_array[2], U_sop2)
    @test 1e-10 > D(U_sop_array[3], U_sop3)

    t = 0.0
    foreach(evolutionOperatorIterator(Ham, tspan; diagonalize = true, approximate = false)
    ) do U_op
        U_op_ref = evolutionOperator(Ham, t)
        @test 1e-11 > D(U_op, U_op_ref)
        t += 0.5
    end

    t = 0.0
    foreach(
        evolutionOperatorIterator(Ham, tspan; diagonalize = false, approximate = false),
    ) do U_op
        U_op_ref = evolutionOperator(Ham, t)
        @test 1e-12 > D(U_op, U_op_ref)
        t += 0.5
    end

    t = 0.0
    foreach(
        evolutionOperatorIterator(Ham, tspan; diagonalize = false, approximate = true),
    ) do U_op
        U_op_ref = evolutionOperator(Ham, t)
        @test 1e-12 > D(U_op, U_op_ref)
        t += 0.5
    end

    t = 0.0
    foreach(
        evolutionSuperOperatorIterator(Ham, tspan; diagonalize = false, approximate = true)
    ) do U_sop
        U_sop_ref = evolutionSuperOperator(Ham, t)
        @test 1e-10 > D(U_sop, U_sop_ref)
        t += 0.5
    end

    t = 0.0
    foreach(
        evolutionSuperOperatorIterator(Ham, tspan; diagonalize = false, approximate = false),
    ) do U_op
        U_op_ref = evolutionSuperOperator(Ham, t)
        @test 1e-12 > D(U_op, U_op_ref)
        t += 0.5
    end

    t = 0.0
    foreach(
        evolutionSuperOperatorIterator(Ham, tspan; diagonalize = false, approximate = true),
    ) do U_op
        U_op_ref = evolutionSuperOperator(Ham, t)
        @test 1e-12 > D(U_op, U_op_ref)
        t += 0.5
    end

    ket = randstate(agg.tools.basis)
    ket1 = U_op1 * ket
    ket2 = U_op2 * ket
    ket3 = U_op3 * ket
    ket_array = evolutionExact(ket, tspan, Ham)
    @test 1e-12 > D(ket1, ket_array[1])
    @test 1e-12 > D(ket2, ket_array[2])
    @test 1e-12 > D(ket3, ket_array[3])

    ket_array = Array{Array{ComplexF32,1},1}(undef, 0)
    evolutionExact!(ket_array, ket, tspan, Ham)
    @test 1e-7 > D(ket1.data, ket_array[1])
    @test 1e-7 > D(ket2.data, ket_array[2])
    @test 1e-7 > D(ket3.data, ket_array[3])

    ket_array = evolutionApproximate(ket, tspan, Ham)
    @test 1e-12 > D(ket1, ket_array[1])
    @test 1e-12 > D(ket2, ket_array[2])
    @test 1e-12 > D(ket3, ket_array[3])

    ket_array = Array{Array{ComplexF32,1},1}(undef, 0)
    evolutionApproximate!(ket_array, ket, tspan, Ham)
    @test 1e-7 > D(ket1.data, ket_array[1])
    @test 1e-7 > D(ket2.data, ket_array[2])
    @test 1e-7 > D(ket3.data, ket_array[3])

    op = dm(ket)
    op1 = U_sop1 * op
    op2 = U_sop2 * op
    op3 = U_sop3 * op
    op_array = evolutionExact(op, tspan, Ham)
    @test 1e-12 > D(op1, op_array[1])
    @test 1e-12 > D(op2, op_array[2])
    @test 1e-12 > D(op3, op_array[3])

    op_array = Array{Array{ComplexF32,2},1}(undef, 0)
    evolutionExact!(op_array, op, tspan, Ham)
    @test 1e-7 > D(op1.data, op_array[1])
    @test 1e-7 > D(op2.data, op_array[2])
    @test 1e-7 > D(op3.data, op_array[3])

    op_array = evolutionApproximate(op, tspan, Ham)
    @test 1e-12 > D(op1, op_array[1])
    @test 1e-12 > D(op2, op_array[2])
    @test 1e-12 > D(op3, op_array[3])

    op_array = Array{Array{ComplexF32,2},1}(undef, 0)
    evolutionApproximate!(op_array, op, tspan, Ham)
    @test 1e-7 > D(op1.data, op_array[1])
    @test 1e-7 > D(op2.data, op_array[2])
    @test 1e-7 > D(op3.data, op_array[3])

    data = Matrix(Hermitian(rand(ComplexF64, agg.tools.bSize, agg.tools.bSize)))
    rho0 = DenseOperator(agg.tools.basis, agg.tools.basis, data)
    normalize!(rho0)

    tspan = [0.0:0.1:1.0;]
    rho_t_ref = evolutionExact(rho0, tspan, Ham)

    T, rho_t = evolution_exact(rho0, tspan, Ham; diagonalize = false)
    U_op_array = evolutionOperatorArray(Ham, tspan)
    for t_i = 1:length(tspan)
        U_op = U_op_array[t_i]
        rho_ref = U_op * rho0 * U_op'
        @test 1e-12 > D(rho_t_ref[t_i], rho_t[t_i])
    end

    T, rho_t = evolution_exact(rho0, tspan, Ham; diagonalize = true)
    for t_i = 1:length(tspan)
        U_op = U_op_array[t_i]
        rho_ref = U_op * rho0 * U_op'
        @test 1e-12 > D(rho_t_ref[t_i], rho_t[t_i])
    end

    T, rho_t = evolution_approximate(rho0, tspan, Ham; diagonalize = false)
    for t_i = 1:length(tspan)
        U_op = U_op_array[t_i]
        rho_ref = U_op * rho0 * U_op'
        @test 1e-12 > D(rho_t_ref[t_i], rho_t[t_i])
    end

    T, rho_t = evolution_approximate(rho0, tspan, Ham; diagonalize = true)
    for t_i = 1:length(tspan)
        U_op = U_op_array[t_i]
        rho_ref = U_op * rho0 * U_op'
        @test 1e-12 > D(rho_t_ref[t_i], rho_t[t_i])
    end

    t = 1.
    indicesMap = agg.tools.indicesMap
    U_part = evolution_el_part(Ham, t, 2, 2, indicesMap)
    data = take_el_part(Ham.data, 2, 2, indicesMap)
    b = GenericBasis([size(data, 1)])
    data = evolutionOperator(data, t)
    ref = DenseOperator(b, b, data)
    @test 1e-12 > D(ref, U_part)


    Evolution_SI_exact()
end
