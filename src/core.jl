
const ComputableType = Union{AbstractFloat,Complex}
const version = "0.2.0"

const ELECTRONIC_GROUND = 1
const ELECTRONIC_EXCITED = 2
const N_ELECTRONIC_LEVELS = 2

const _SAFE_DIV_TOL = 1e-10

_safe_div(a, b; tol=_SAFE_DIV_TOL) = abs(b) < tol ? zero(a) : a / b

_safe_inv(x; tol=_SAFE_DIV_TOL) = abs(x) < tol ? zero(x) : 1 / x

abstract type AbstractVibBasis end
struct GroundGround <: AbstractVibBasis end
struct GroundExcited <: AbstractVibBasis end

const VibBasisLike = Union{Symbol, AbstractVibBasis}

function _to_vib_basis(s::Symbol)
    s === :ground_ground && return GroundGround()
    s === :ground_excited && return GroundExcited()
    s === :none && return s
    throw(ArgumentError("vib_basis must be :ground_ground or :ground_excited, got :$s"))
end
_to_vib_basis(v::AbstractVibBasis) = v

_vib_basis_to_symbol(::GroundGround) = :ground_ground
_vib_basis_to_symbol(::GroundExcited) = :ground_excited
_vib_basis_to_symbol(s::Symbol) = s