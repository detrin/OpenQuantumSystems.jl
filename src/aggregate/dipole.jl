
# Prefactor converting D²/Å³ to cm⁻¹ for the point-dipole coupling formula.
# J [cm⁻¹] = DIPOLE_COUPLING_FACTOR × (μ₁ · μ₂ - 3(μ₁·R̂)(μ₂·R̂)) / R³
# with dipole magnitudes in Debye and R in Å.
# Derived from 1/(4πε₀) in SI with unit conversion.
const DIPOLE_COUPLING_FACTOR = 5034.5

"""
    TransitionDipole(position, dipole)

Stores the 3D position and transition dipole vector of a molecule, used to
compute point-dipole couplings between molecules in an aggregate.

# Arguments
* `position`: 3-element vector giving the molecular position in Å.
* `dipole`: 3-element transition dipole vector in Debye (magnitude encoded in norm).
"""
struct TransitionDipole
    position::Vector{Float64}
    dipole::Vector{Float64}
    function TransitionDipole(position::AbstractVector{<:Real}, dipole::AbstractVector{<:Real})
        length(position) == 3 || throw(ArgumentError("position must be a 3-element vector"))
        length(dipole) == 3 || throw(ArgumentError("dipole must be a 3-element vector"))
        new(collect(Float64, position), collect(Float64, dipole))
    end
end

Base.:(==)(a::TransitionDipole, b::TransitionDipole) =
    a.position == b.position && a.dipole == b.dipole

"""
    dipole_dipole_coupling(td1, td2)

Compute the point-dipole coupling (in cm⁻¹) between two transition dipoles.

Uses the formula:
```
J = DIPOLE_COUPLING_FACTOR × (μ₁·μ₂ - 3(μ₁·R̂)(μ₂·R̂)) / R³
```
where R is the inter-molecular distance in Å, μ in Debye, and the prefactor
converts D²/Å³ to cm⁻¹.

# Arguments
* `td1`, `td2`: [`TransitionDipole`](@ref) instances.
"""
function dipole_dipole_coupling(td1::TransitionDipole, td2::TransitionDipole)::Float64
    R_vec = td2.position .- td1.position
    R = norm(R_vec)
    iszero(R) && throw(ArgumentError("Positions of td1 and td2 are identical; R = 0"))
    R_hat = R_vec ./ R
    mu1 = td1.dipole
    mu2 = td2.dipole
    kappa_mu = dot(mu1, mu2) - 3 * dot(mu1, R_hat) * dot(mu2, R_hat)
    return DIPOLE_COUPLING_FACTOR * kappa_mu / R^3
end

"""
    coupling_from_dipoles(dipoles)

Build an N×N coupling matrix (in cm⁻¹) from a vector of [`TransitionDipole`](@ref)s.
Diagonal entries are zero. Off-diagonal entry (i,j) is the point-dipole coupling
between molecule i and molecule j.

# Arguments
* `dipoles`: Vector of [`TransitionDipole`](@ref) instances, one per molecule.
"""
function coupling_from_dipoles(dipoles::Vector{TransitionDipole})::Matrix{Float64}
    n = length(dipoles)
    J = zeros(Float64, n, n)
    for i in 1:n
        for j in (i+1):n
            J[i, j] = dipole_dipole_coupling(dipoles[i], dipoles[j])
            J[j, i] = J[i, j]
        end
    end
    return J
end

"""
    AggregateCore(molecules, dipoles)

Construct an [`AggregateCore`](@ref) where the coupling matrix is computed
automatically from point-dipole interactions.

# Arguments
* `molecules`: Vector of [`Molecule`](@ref) instances.
* `dipoles`: Vector of [`TransitionDipole`](@ref) instances (one per molecule).
"""
function AggregateCore(
    molecules::Vector{Molecule{T,C1,C2}},
    dipoles::Vector{TransitionDipole},
) where {T,C1,C2}
    n = length(molecules)
    length(dipoles) == n || throw(
        DimensionMismatch(
            "Number of dipoles ($(length(dipoles))) must match number of molecules ($n)"
        )
    )
    J = coupling_from_dipoles(dipoles)
    return AggregateCore(molecules, J)
end
