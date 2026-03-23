using Test
using LinearAlgebra
using OpenQuantumSystems

@testset "dipole" begin

    @testset "TransitionDipole constructor" begin
        td = TransitionDipole([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        @test td.position == [0.0, 0.0, 0.0]
        @test td.dipole == [1.0, 0.0, 0.0]

        # accepts any AbstractVector{<:Real}
        td2 = TransitionDipole([0, 5, 0], [0, 1, 0])
        @test td2.position == [0.0, 5.0, 0.0]

        @test_throws ArgumentError TransitionDipole([0.0, 0.0], [1.0, 0.0, 0.0])
        @test_throws ArgumentError TransitionDipole([0.0, 0.0, 0.0], [1.0, 0.0])
    end

    @testset "dipole_dipole_coupling — geometry" begin
        # Parallel dipoles (both along y), R along x.
        # κ = (ŷ·ŷ) - 3(ŷ·x̂)(ŷ·x̂) = 1 - 0 = 1
        # J = DIPOLE_COUPLING_FACTOR × 1 × 1 × 1 / 10³
        td1 = TransitionDipole([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        td2 = TransitionDipole([10.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        J_parallel = dipole_dipole_coupling(td1, td2)
        @test J_parallel ≈ DIPOLE_COUPLING_FACTOR / 1000.0

        # Head-to-tail (both along x), R along x.
        # κ = (x̂·x̂) - 3(x̂·x̂)² = 1 - 3 = -2
        # J = DIPOLE_COUPLING_FACTOR × (-2) / 10³
        td3 = TransitionDipole([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        td4 = TransitionDipole([10.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        J_head_to_tail = dipole_dipole_coupling(td3, td4)
        @test J_head_to_tail ≈ -2 * DIPOLE_COUPLING_FACTOR / 1000.0

        # Head-to-tail coupling is exactly -2× parallel coupling (same R, same |μ|)
        @test J_head_to_tail ≈ -2 * J_parallel

        # Perpendicular dipoles, R perpendicular to both: κ = 0
        # μ1 along x, μ2 along y, R along z
        td5 = TransitionDipole([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        td6 = TransitionDipole([0.0, 0.0, 10.0], [0.0, 1.0, 0.0])
        @test dipole_dipole_coupling(td5, td6) ≈ 0.0 atol=1e-12

        # Symmetry: J(1→2) == J(2→1)
        @test dipole_dipole_coupling(td1, td2) ≈ dipole_dipole_coupling(td2, td1)
        @test dipole_dipole_coupling(td3, td4) ≈ dipole_dipole_coupling(td4, td3)
    end

    @testset "dipole_dipole_coupling — magnitude scaling" begin
        td1 = TransitionDipole([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        td2 = TransitionDipole([10.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        td3 = TransitionDipole([0.0, 0.0, 0.0], [0.0, 2.0, 0.0])
        td4 = TransitionDipole([10.0, 0.0, 0.0], [0.0, 2.0, 0.0])

        # Doubling both dipole magnitudes → 4× coupling
        @test dipole_dipole_coupling(td3, td4) ≈ 4 * dipole_dipole_coupling(td1, td2)

        # Distance scaling: halving R → 8× coupling
        td5 = TransitionDipole([5.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        @test dipole_dipole_coupling(td1, td5) ≈ 8 * dipole_dipole_coupling(td1, td2)
    end

    @testset "dipole_dipole_coupling — zero distance error" begin
        td = TransitionDipole([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        @test_throws ArgumentError dipole_dipole_coupling(td, td)
    end

    @testset "coupling_from_dipoles" begin
        dipoles = [
            TransitionDipole([0.0,  0.0, 0.0], [0.0, 1.0, 0.0]),
            TransitionDipole([10.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            TransitionDipole([20.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ]
        J = coupling_from_dipoles(dipoles)

        # Diagonal is zero
        @test J[1, 1] == 0.0
        @test J[2, 2] == 0.0
        @test J[3, 3] == 0.0

        # Symmetric
        @test J[1, 2] ≈ J[2, 1]
        @test J[1, 3] ≈ J[3, 1]
        @test J[2, 3] ≈ J[3, 2]

        # Nearest-neighbour (R=10) vs next-nearest (R=20): ratio = 8
        @test J[1, 2] ≈ 8 * J[1, 3]

        # Correct values
        @test J[1, 2] ≈ DIPOLE_COUPLING_FACTOR / 1000.0
        @test J[1, 3] ≈ DIPOLE_COUPLING_FACTOR / 8000.0
    end

    @testset "AggregateCore from dipoles" begin
        modes = [Mode(200.0, 0.2)]
        mol1 = Molecule(modes, 3, [0.0, 12500.0])
        mol2 = Molecule(modes, 3, [0.0, 12700.0])
        molecules = [mol1, mol2]

        dipoles = [
            TransitionDipole([0.0,  0.0, 0.0], [0.0, 1.0, 0.0]),
            TransitionDipole([10.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ]
        aggCore = AggregateCore(molecules, dipoles)

        @test aggCore isa AggregateCore
        @test aggCore.molCount == 2

        J_expected = DIPOLE_COUPLING_FACTOR / 1000.0
        # Coupling stored in padded (molCount+1)×(molCount+1) matrix, molecules at 2:end
        @test aggCore.coupling[2, 3] ≈ J_expected
        @test aggCore.coupling[3, 2] ≈ J_expected
        @test aggCore.coupling[1, 1] == 0.0

        # Mismatch error
        @test_throws DimensionMismatch AggregateCore(molecules, [dipoles[1]])
    end

end
