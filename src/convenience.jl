"""
    setup_dimer(; E1=12500., E2=12700., J=100., mode_omega=200., mode_hr=0.02, nvib=3)

Create a fully set-up dimer [`Aggregate`](@ref) with one vibrational mode per molecule.

Returns an [`Aggregate`](@ref) ready for time-evolution calculations.
"""
function setup_dimer(;
    E1::Float64 = 12500.0,
    E2::Float64 = 12700.0,
    J::Float64 = 100.0,
    mode_omega::Float64 = 200.0,
    mode_hr::Float64 = 0.02,
    nvib::Int = 3,
)
    mode1 = Mode(; omega = mode_omega, hr_factor = mode_hr)
    mode2 = Mode(; omega = mode_omega, hr_factor = mode_hr)
    mol1 = Molecule([mode1], nvib, [0.0, E1])
    mol2 = Molecule([mode2], nvib, [0.0, E2])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = J
    aggCore.coupling[3, 2] = J
    return setupAggregate(aggCore)
end

"""
    setup_trimer(; E=[12500., 12700., 12900.], J_matrix=zeros(3,3), mode_omega=200., mode_hr=0.02, nvib=3)

Create a fully set-up trimer [`Aggregate`](@ref) with one vibrational mode per molecule.

`J_matrix` is a 3x3 coupling matrix (only off-diagonal elements are used).
Returns an [`Aggregate`](@ref) ready for time-evolution calculations.
"""
function setup_trimer(;
    E::Vector{Float64} = [12500.0, 12700.0, 12900.0],
    J_matrix::Matrix{Float64} = zeros(Float64, 3, 3),
    mode_omega::Float64 = 200.0,
    mode_hr::Float64 = 0.02,
    nvib::Int = 3,
)
    length(E) == 3 || throw(ArgumentError("E must have exactly 3 elements"))
    size(J_matrix) == (3, 3) || throw(ArgumentError("J_matrix must be 3x3"))

    mols = [
        Molecule([Mode(; omega = mode_omega, hr_factor = mode_hr)], nvib, [0.0, E[i]])
        for i in 1:3
    ]
    aggCore = AggregateCore(mols)
    for i in 1:3, j in 1:3
        if i != j
            aggCore.coupling[i+1, j+1] = J_matrix[i, j]
        end
    end
    return setupAggregate(aggCore)
end

"""
    setup_linear_chain(; energies=[12500., 12700.], J_nearest=100., mode_omega=200., mode_hr=0.02, nvib=3)

Create a fully set-up linear chain [`Aggregate`](@ref) with nearest-neighbour coupling
and one vibrational mode per molecule.

Returns an [`Aggregate`](@ref) ready for time-evolution calculations.
"""
function setup_linear_chain(;
    energies::Vector{Float64} = [12500.0, 12700.0],
    J_nearest::Float64 = 100.0,
    mode_omega::Float64 = 200.0,
    mode_hr::Float64 = 0.02,
    nvib::Int = 3,
)
    N = length(energies)
    N >= 2 || throw(ArgumentError("energies must have at least 2 elements"))

    mols = [
        Molecule([Mode(; omega = mode_omega, hr_factor = mode_hr)], nvib, [0.0, energies[i]])
        for i in 1:N
    ]
    aggCore = AggregateCore(mols)
    for i in 1:(N-1)
        aggCore.coupling[i+1, i+2] = J_nearest
        aggCore.coupling[i+2, i+1] = J_nearest
    end
    return setupAggregate(aggCore)
end
