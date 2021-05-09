
# import Adapt

# Convert data to CuArray with cu(::Operator)
# Adapt.adapt_structure(to, x::Operator) = Operator(x.basis_l, x.basis_r, Adapt.adapt(to, x.data))

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

function OneDenseOperator(bl::BL, br::BR) where {BL<:Basis,BR<:Basis}
    op = DenseOperator(bl, br)
    opLen = size(op.data, 1)
    for i = 1:opLen
        op.data[i, i] = 1.0
    end
    return op
end

OneDenseOperator(b::B) where {B<:Basis} = OneDenseOperator(b, b)
