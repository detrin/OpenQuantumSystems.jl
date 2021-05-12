using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra
import QuantumOpticsBase


@testset "memory kernel" begin

    Random.seed!(0)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))

    mode1 = Mode(0.2, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 2, Energy)
    mol2 = Molecule([mode1], 2, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg; groundState = false)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds; groundState = false)
    FCProd = getFCProd(agg, FCFact, aggInds, vibindices; groundState = false)
    Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState = false)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundState = false)
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggInds, FCFact; groundState = false)
    Ham_S = Ham - Ham_int

    T = 300
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, aggInds; diagonalize=false, diagonal=true)

    W_ab = take_el_part(W0.data, 2, 2, vibindices)
    W_ab_len = size(W_ab, 1)
    @test 1e-15 > D(W_ab, W0.data[1:W_ab_len, 1:W_ab_len])

end # testset
