using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs
import QuantumOpticsBase

@testset "initial state" begin

    Random.seed!(StableRNG(0), 1)

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
    aggInds = getIndices(agg)
    vibindices = getVibIndices(agg, aggInds)
    aggIndsLen = length(aggInds)
    basis = GenericBasis([aggIndsLen])
    FCFact = getFranckCondonFactors(agg, aggInds)
    Ham = getAggHamiltonian(agg, aggInds, FCFact)

    T = 300.0
    mu_array = [[1, 2]]
    W0_ref = [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25023985364402346 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24999994249781035 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24999994249781035 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2497602613603558
    ]
    # resulting in NaNs
    #=
    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_composite_ref, W0)
    =# 

    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = true)
    @test 1e-3 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = true)
    @test 1e-15 > D(W0_composite_ref, W0)

    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        Ham,
        vibindices,
        aggInds;
        diagonalize = false,
        diagonal = true,
    )
    @test 1e-15 > D(W0_composite_ref, W0)

    T = 1e-10
    W0_ref = [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    # resulting in NaNs
    #=
    println(Ham)
    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false)
    println(W0)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_composite_ref, W0)
    

    # W0 = thermal_state(T, mu_array, Ham, aggInds; diagonalize=true)
    # @test 1e-15 > D(W0_ref, W0.data)
    # println(W0.data)
    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        Ham,
        vibindices,
        aggInds;
        diagonalize = false,
        diagonal = true,
    )
    @test 1e-15 > D(W0_composite_ref, W0)
    =#

    T = 300.0
    mu_array = [[1, 2], [2, 1]]
    W0_ref = [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25023985364402346 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24999994249781035 0.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24999994249781035 0.0; 
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24976026136035578
    ]
    # resulting in NaNs
    #=
    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_composite_ref, W0)
    =#

    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = true)
    @test 1e-3 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = true)
    @test 1e-15 > D(W0_composite_ref, W0)

    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        Ham,
        vibindices,
        aggInds;
        diagonalize = false,
        diagonal = true,
    )
    @test 1e-15 > D(W0_composite_ref, W0)

    aggInds = getIndices(agg)
    FCFact = getFranckCondonFactors(agg, aggInds)
    Ham = getAggHamiltonian(agg, aggInds, FCFact)
    T = 300.0
    mu_array = [[1, 1]]
    W0_ref = [
        0.25023985364402346 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.24999994249781035 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.24999994249781035 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.24976026136035578 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    #=
    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = false)
    @test 1e-15 > D(W0_composite_ref, W0)
    =#

    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = true)
    @test 1e-3 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], Ham, vibindices, aggInds; diagonalize = true)
    @test 1e-15 > D(W0_composite_ref, W0)

    W0 = thermal_state(T, mu_array, Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    @test 1e-15 > D(W0_ref, W0.data)
    W01 = thermal_state(T, [[2, 1]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W02 = thermal_state(T, [[1, 2]], Ham, vibindices, aggInds; diagonalize = false, diagonal = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        Ham,
        vibindices,
        aggInds;
        diagonalize = false,
        diagonal = true,
    )
    @test 1e-15 > D(W0_composite_ref, W0)


end # testset
