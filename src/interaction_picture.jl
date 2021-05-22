
function getInteractionHamIPicture(
    Ham_S::Operator,
    Ham_int::Operator,
    t::AbstractFloat,
)::Operator
    U_op = evolutionOperator(Ham_S, t)
    return U_op' * Ham_int * U_op
end

function getInteractionHamIPictureA(
    Ham_S::Operator,
    Ham_int::Operator,
    t::AbstractFloat,
)::Array
    U_op = evolutionOperator(Ham_S, t)
    return (U_op' * Ham_int * U_op).data
end

function getInteractionHamIPictureA(
    Ham_int::Array,
    H_lambda::Array,
    H_S::Array,
    H_Sinv::Array,
    t::AbstractFloat,
)
    U_op = evolutionOperatorA(H_lambda, H_S, H_Sinv, t)
    U_opd = adjoint(U_op)
    return U_opd * Ham_int * U_op
end
