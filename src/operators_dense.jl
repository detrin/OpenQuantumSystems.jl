
# import Adapt

# Convert data to CuArray with cu(::Operator)
# Adapt.adapt_structure(to, x::Operator) = Operator(x.basis_l, x.basis_r, Adapt.adapt(to, x.data))

"""
    AnnihilationOperator{BL,BR}(basis_l, basis_r)

Dense annihilation operator as a mutable struct.

"""
mutable struct AnnihilationOperator{BL<:Basis,BR<:Basis} <: DataOperator{BL,BR}
    basis_l::BL
    basis_r::BR
    data::Matrix{Float64}
    function AnnihilationOperator{BL,BR}(
        basis_l::BL,
        basis_r::BR,
    ) where {BL<:Basis,BR<:Basis}
        data = zeros(Float64, length(basis_l), length(basis_r))
        for i = 1:length(basis_l)
            for j = 1:length(basis_r)
                if i == j - 1
                    data[i, j] = Float64(j - 1)^0.5
                end
            end
        end
        new(basis_l, basis_r, data)
    end
end

AnnihilationOperator(bl::BL, br::BR) where {BL,BR} = AnnihilationOperator{BL,BR}(bl, br)
AnnihilationOperator(b::Basis) = AnnihilationOperator(b, b)

"""
    CreationOperator{BL,BR}(basis_l, basis_r)

Dense creation operator as a mutable struct.

"""
mutable struct CreationOperator{BL<:Basis,BR<:Basis} <: DataOperator{BL,BR}
    basis_l::BL
    basis_r::BR
    data::Matrix{Float64}
    function CreationOperator{BL,BR}(basis_l::BL, basis_r::BR) where {BL<:Basis,BR<:Basis}
        data = zeros(Float64, length(basis_l), length(basis_r))
        for i = 1:length(basis_l)
            for j = 1:length(basis_r)
                if i == j + 1
                    data[i, j] = Float64(j)^0.5
                end
            end
        end
        new(basis_l, basis_r, data)
    end
end

CreationOperator(bl::BL, br::BR) where {BL,BR} = CreationOperator{BL,BR}(bl, br)
CreationOperator(b::Basis) = CreationOperator(b, b)

"""
    PositionOperator{BL,BR}(basis_l, basis_r)

Dense position operator as a mutable struct without mass and frequency.

"""
mutable struct PositionOperator{BL<:Basis,BR<:Basis} <: DataOperator{BL,BR}
    basis_l::BL
    basis_r::BR
    data::Matrix{Float64}
    function PositionOperator{BL,BR}(
        basis_l::BL,
        basis_r::BR,
    ) where {BL<:Basis,BR<:Basis}
        # data = zeros(Float64, length(basis_l), length(basis_r))
        creation_op = CreationOperator(basis_l, basis_r)
        annihilation_op = AnnihilationOperator(basis_l, basis_r)
        data = (creation_op.data + annihilation_op.data)/sqrt(2.0)
        # new(basis_l, basis_r, data)
        new(basis_l, basis_r, data)
    end
end

PositionOperator(bl::BL, br::BR) where {BL,BR} = PositionOperator{BL,BR}(bl, br)
PositionOperator(b::Basis) = PositionOperator(b, b)

"""
MomentumOperator{BL,BR}(basis_l, basis_r)

Dense momentum operator as a mutable struct without mass and frequency.

"""
mutable struct MomentumOperator{BL<:Basis,BR<:Basis} <: DataOperator{BL,BR}
    basis_l::BL
    basis_r::BR
    data::Matrix{ComplexF64}
    function MomentumOperator{BL,BR}(
        basis_l::BL,
        basis_r::BR,
    ) where {BL<:Basis,BR<:Basis}
        # data = zeros(ComplexF64, length(basis_l), length(basis_r))
        creation_op = CreationOperator(basis_l, basis_r)
        annihilation_op = AnnihilationOperator(basis_l, basis_r)
        data = 1im*(creation_op.data - annihilation_op.data)/sqrt(2.0)
        # new(basis_l, basis_r, data)
        new(basis_l, basis_r, data)
    end
end

MomentumOperator(bl::BL, br::BR) where {BL,BR} = MomentumOperator{BL,BR}(bl, br)
MomentumOperator(b::Basis) = MomentumOperator(b, b)

"""
    ShiftOperator{BL,BR}(basis_l, basis_r, shift)

Dense shift operator as a mutable struct using the definition

``D(\\alpha) = \\exp(\\alpha a^\\dagger - \\alpha^* a)``.

# Arguments
* `basis_l`: Bra basis.
* `basis_r`: Ket basis.
* `shift`: Shift or ``\\alpha`` parameter, can be complex number.
"""
mutable struct ShiftOperator{BL<:Basis,BR<:Basis,T<:ComputableType} <: DataOperator{BL,BR}
    basis_l::BL
    basis_r::BR
    data::Matrix{ComplexF64}
    shift::T
    function ShiftOperator{BL,BR,T}(
        basis_l::BL,
        basis_r::BR,
        shift::T,
    ) where {BL<:Basis,BR<:Basis,T<:ComputableType}
        data = zeros(ComplexF64, length(basis_l), length(basis_r))
        creation_op = CreationOperator(basis_l, basis_r)
        annihilation_op = AnnihilationOperator(basis_l, basis_r)
        data = (shift * creation_op.data - conj(shift) * annihilation_op.data) / 2^0.5
        data = exp(data)
        new(basis_l, basis_r, data, shift)
    end
end

ShiftOperator(bl::BL, br::BR, shift::T) where {BL,BR,T} =
    ShiftOperator{BL,BR,T}(bl, br, shift)
ShiftOperator(b::Basis, shift::ComputableType) = ShiftOperator(b, b, shift)

"""
    OneDenseOperator(basis_l, basis_r)

DenseOperator with ones on the diagonal.

"""
function OneDenseOperator(basis_l::BL, basis_r::BR) where {BL<:Basis,BR<:Basis}
    op = DenseOperator(basis_l, basis_r)
    opLen = size(op.data, 1)
    for i = 1:opLen
        op.data[i, i] = 1.0
    end
    return op
end

OneDenseOperator(b::B) where {B<:Basis} = OneDenseOperator(b, b)
