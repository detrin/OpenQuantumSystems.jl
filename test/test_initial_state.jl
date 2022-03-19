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
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 2, [2.0, 200.0])
    mol2 = Molecule([mode2], 2, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)
    aggTools = agg.tools
    aggOperators = agg.operators

    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct

    T = 300.0
    mu_array = [[1, 1, 2]]
    W0_ref = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25000143755694526 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2499997124870238 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2500002875090083 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24999856244702248]
    W0 = thermal_state(T, mu_array, aggCore, aggTools, aggOperators; diagonalize = true)
    @test 1e-15 > D(W0_ref, W0.data)

    W0 = thermal_state(T, mu_array, agg; diagonalize = true)
    @test 1e-15 > D(W0_ref, W0.data)

    W01 = thermal_state(T, [[1, 2, 1]], aggCore, aggTools, aggOperators; diagonalize = true)
    W02 = thermal_state(T, [[1, 1, 2]], aggCore, aggTools, aggOperators; diagonalize = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], aggCore, aggTools, aggOperators; diagonalize = true)
    @test 1e-15 > D(W0_composite_ref, W0)

    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], agg; diagonalize = true)
    @test 1e-15 > D(W0_composite_ref, W0)

    W0_ref = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25029983141333245 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24993996472850133 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25005986276494624 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2497003410932199]
    W0 = thermal_state(T, mu_array, aggCore, aggTools, aggOperators; diagonalize = false, diagonal = true)
    @test 1e-15 > D(W0_ref, W0.data)

    W01 = thermal_state(T, [[1, 2, 1]], aggCore, aggTools, aggOperators; diagonalize = false, diagonal = true)
    W02 = thermal_state(T, [[1, 1, 2]], aggCore, aggTools, aggOperators; diagonalize = false, diagonal = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    normalize!(W0_composite_ref)
    W0 = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        aggCore,
        aggTools,
        aggOperators;
        diagonalize = false,
        diagonal = true,
    )
    @test 1e-15 > D(W0_composite_ref, W0)

    W01 = thermal_state(T, [[1, 2, 1]], aggCore, aggTools, aggOperators; diagonalize = true)
    W02 = thermal_state(T, [[1, 1, 2]], aggCore, aggTools, aggOperators; diagonalize = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    # normalize!(W0_composite_ref)
    W0 = thermal_state_composite(T, [0.0, 0.8, 0.2], aggCore, aggTools, aggOperators; diagonalize = true)
    @test 1e-15 > D(W0_composite_ref, W0)

    W0_ref = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25029983141333245 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24993996472850133 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25005986276494624 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2497003410932199]
    W0 = thermal_state(T, mu_array, aggCore, aggTools, aggOperators; diagonalize = false, diagonal = true)
    @test 1e-15 > D(W0_ref, W0.data)
    
    W01 = thermal_state(T, [[1, 2, 1]], aggCore, aggTools, aggOperators; diagonalize = false, diagonal = true)
    W02 = thermal_state(T, [[1, 1, 2]], aggCore, aggTools, aggOperators; diagonalize = false, diagonal = true)
    W0_composite_ref = 0.8 * W01 + 0.2 * W02
    # normalize!(W0_composite_ref)
    W0 = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        aggCore,
        aggTools,
        aggOperators;
        diagonalize = false,
        diagonal = true,
    )
    @test 1e-15 > D(W0_composite_ref, W0)

    T = 0.5
    mu_array = [[1, 1, 2]]
    W0_ref = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.7750674418239181 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.06464140010302448 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1479518276234691 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.012339330449588293]
    W0 = thermal_state(T, mu_array, aggCore, aggTools, aggOperators; diagonalize = true)
    @test 1e-15 > D(W0_ref, W0.data)

    W0, rho0, W0_bath = ultrafast_laser_excitation(T, [0.0, 0.8, 0.2], agg; diagonalize = true)
    W0_ref = thermal_state_composite(
        T,
        [0.0, 0.8, 0.2],
        aggCore,
        aggTools,
        aggOperators;
        diagonalize = true
    )
    normalize!(W0_ref)
    W0_ref = DenseOperator(W0_ref.basis_l, W0_ref.basis_r, complex(W0_ref.data))

    rho0_ref = trace_bath(W0_ref, agg.core, agg.tools)
    rho0_ref = DenseOperator(rho0_ref.basis_l, rho0_ref.basis_r, complex(rho0_ref.data))

    W0_bath_ref = get_rho_bath(W0_ref, agg.core, agg.tools)
    W0_bath_ref = DenseOperator(W0_bath_ref.basis_l, W0_bath_ref.basis_r, complex(W0_bath_ref.data))
    @test 1e-15 > D(W0_ref, W0)
    @test 1e-15 > D(rho0_ref, rho0)
    @test 1e-15 > D(W0_bath_ref, W0_bath)

end # testset
