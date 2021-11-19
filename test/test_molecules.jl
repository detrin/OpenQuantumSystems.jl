using Test
using OpenQuantumSystems
using LinearAlgebra
# using Traceur

@testset "molecules" begin

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))

    fc1 = franckCondonFactors(3, 0.0)
    fc2 = [
        1.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im 1.0+0.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im 1.0+0.0im
    ]
    @test 1e-12 > D(fc1, fc2)

    fc1 = franckCondonFactors(3, 1.0)
    fc2 = [
        0.7788007830714049+0.0im -0.5506953149031837+0.0im 0.2753476574515917+0.0im
        0.5506953149031838+0.0im 0.3894003915357024+0.0im -0.5841005873035535+0.0im
        0.27534765745159184+0.0im 0.5841005873035536+0.0im 0.0973500978839258+0.0im
    ]
    @test 1e-12 > D(fc1, fc2)

    fc1 = franckCondonFactors(3, 2.0im)
    fc2 = franckCondonFactors(3, -2.0im)
    fc2 .= conj(fc2)
    @test 1e-12 > D(fc1, fc2)

    @test vibrationalIndices([2]) == [[1], [2]]

    inds = [[1, 1], [1, 2], [2, 1], [2, 2]]
    @test vibrationalIndices([2, 2]) == inds

    inds = [
        [1, 1],
        [1, 2],
        [1, 3],
        [2, 1],
        [2, 2],
        [2, 3],
        [3, 1],
        [3, 2],
        [3, 3],
        [4, 1],
        [4, 2],
        [4, 3],
    ]
    @test vibrationalIndices([4, 3]) == inds
    # @trace vibrationalIndices([4, 3])
    # println("dispatch 2")

    @test electronicIndices(3) == [[1, 1, 1], [2, 1, 1], [1, 2, 1], [1, 1, 2]]

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.4, 1.0)
    Energy = [0.0, 200.0]
    fcFactors = [franckCondonFactors(3, mode1.shift), franckCondonFactors(3, mode2.shift)]
    mol = Molecule([mode1, mode2], 3, Energy)
    @test mol.modes == [mode1, mode2]
    @test mol.Nvib == 3
    @test 1e-12 > D(mol.fcFactors, fcFactors)
    @test mol.E == Energy

    @test getMolStateEnergy(mol, 1, [1, 1]) ≈ 0.3
    @test getMolStateEnergy(mol, 1, [3, 1]) ≈ 0.7
    @test getMolStateEnergy(mol, 2, [2, 2]) ≈ 200.9

    @test getMolStateFC(mol, 1, [1, 1], 1, [1, 1]) == 1.0
    @test getMolStateFC(mol, 1, [2, 1], 1, [1, 2]) == 0.0
    @test getMolStateFC(mol, 1, [1, 1], 2, [1, 1]) ≈ 0.6065306597126334
    @test getMolStateFC(mol, 2, [3, 1], 2, [1, 2]) == 0.0

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.4, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1, mode2], 3, Energy)
    mol1.E = [0.0, 400.0]
    mol2 = updateMolecule(mol1)
    @test mol2.E == [0.0, 400.0]
    mol2.E = [0.0, 200.0]
    updateMolecule!(mol2)
    mol1.E = [0.0, 200.0]
    mol1 = updateMolecule(mol1)
    @test mol1.E == mol2.E

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.4, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1, mode2], 3, Energy)

    @test getMolShifts(mol1) == [1.0, 2.0]
    @test getMolFrequencies(mol1) == [0.2, 0.4]
    
end
