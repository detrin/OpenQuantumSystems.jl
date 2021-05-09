
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
