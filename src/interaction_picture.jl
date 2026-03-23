
function getInteractionHamIPicture(
    Ham_0::AbstractMatrix,
    Ham_I::AbstractMatrix,
    t::AbstractFloat,
)::AbstractMatrix
    U_array = evolutionOperator(Ham_0, t)
    return adjoint(U_array) * Ham_I * U_array
end

function getInteractionHamIPicture(
    Ham_0::Operator,
    Ham_I::Operator,
    t::AbstractFloat,
)::Operator
    U_op = evolutionOperator(Ham_0, t)
    return U_op' * Ham_I * U_op
end

function getInteractionHamIPictureA(
    Ham_0::Operator,
    Ham_I::Operator,
    t::AbstractFloat,
)::AbstractMatrix
    U_op = evolutionOperator(Ham_0, t)
    return (U_op' * Ham_I * U_op).data
end

function getInteractionHamIPictureA(
    Ham_int::AbstractMatrix,
    H_lambda::AbstractMatrix,
    H_S::AbstractMatrix,
    H_Sinv::AbstractMatrix,
    t::AbstractFloat,
)
    U_op = evolutionOperatorA(H_lambda, H_S, H_Sinv, t)
    U_opd = adjoint(U_op)
    return U_opd * Ham_int * U_op
end
