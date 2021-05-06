using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra

@testset "trace" begin

    Random.seed!(0)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))

    mode1 = Mode(0.2, 1.)
    Energy = [0., 200.]
    mol1 = Molecule([mode1], 3, Energy)
    mol2 = Molecule([mode1], 3, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState=false)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState=false)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState=false)
    
    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.5373667261645646 + 0.0im 0.5960193120071318 + 0.31849769746151213im; 
        0.5960193120071318 - 0.31849769746151213im 0.4626332738354356 + 0.0im
    ]
    rho_traced = trace_bath(rho, agg, FCFact, aggInds, vibindices; groundState=false)
    @test 1e-15 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCFact, aggInds, vibindices; groundState=false)
    @test 1e-15 > D(rho_traced, rho_traced_ref)

    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.5373667261645646 + 4.804317299001507e-18im 0.5960193120071318 + 0.31849769746151213im; 
        0.5960193120071318 - 0.31849769746151213im 0.46263327383543557 - 4.716400113933126e-18im
    ]
    rho_traced = trace_bath(rho, agg, FCFact, aggInds, vibindices; groundState=false)
    @test 1e-15 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCFact, aggInds, vibindices; groundState=false)
    @test 1e-15 > D(rho_traced, rho_traced_ref)


    mode1 = Mode(0.2, 1.)
    Energy = [0., 200.]
    mol1 = Molecule([mode1], 3, Energy)
    mol2 = Molecule([mode1], 3, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState=true)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState=true)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState=true)
    
    data = Matrix(Hermitian(rand(ComplexF64, aggIndsLen, aggIndsLen)))
    rho0 = DenseOperator(basis, basis, data)
    normalize!(rho0)

    t = 0.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.27983831453187563 + 0.0im 0.39510751997977556 + 0.3832496461602654im 0.36717929072360034 + 0.48380689550795175im; 
        0.39510751997977556 - 0.3832496461602654im 0.3757521692494342 + 0.0im 0.4295074106269881 + 0.3124835865738737im; 
        0.36717929072360034 - 0.48380689550795175im 0.4295074106269881 - 0.3124835865738737im 0.3444095162186903 + 0.0im
    ]
    rho_traced = trace_bath(rho, agg, FCFact, aggInds, vibindices; groundState=true)
    @test 1e-15 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCFact, aggInds, vibindices; groundState=true)
    @test 1e-15 > D(rho_traced, rho_traced_ref)

    t = 1.0
    U_op = evolutionOperator(Ham, t)
    rho = U_op * rho0 * U_op'
    rho_traced_ref = ComplexF64[
        0.27983831453187474 + 4.807665611244999e-18im 0.5271823941867131 - 0.15833182524717854im 0.6013924791789407 - 0.08495192560698843im; 
        0.5271823941867131 + 0.15833182524717854im 0.3757521692494339 + 5.194043729796511e-18im 0.4295074106269877 + 0.31248358657387476im; 
        0.6013924791789407 + 0.08495192560698844im 0.4295074106269877 - 0.3124835865738748im 0.34440951621869137 + 4.376703186460492e-18im
    ]

    rho_traced = trace_bath(rho, agg, FCFact, aggInds, vibindices; groundState=true)
    @test 1e-15 > D(rho_traced.data, rho_traced_ref)
    rho_traced = trace_bath(rho.data, agg, FCFact, aggInds, vibindices; groundState=true)
    @test 1e-15 > D(rho_traced, rho_traced_ref)


end