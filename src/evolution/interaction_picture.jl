
function get_interaction_ham_i_picture(
    Ham_0::AbstractMatrix,
    Ham_I::AbstractMatrix,
    t::AbstractFloat,
)::AbstractMatrix
    U_array = evolution_operator(Ham_0, t)
    return adjoint(U_array) * Ham_I * U_array
end

function get_interaction_ham_i_picture(
    Ham_0::Operator,
    Ham_I::Operator,
    t::AbstractFloat,
)::Operator
    U_op = evolution_operator(Ham_0, t)
    return U_op' * Ham_I * U_op
end

function get_interaction_ham_i_picture_a(
    Ham_0::Operator,
    Ham_I::Operator,
    t::AbstractFloat,
)::AbstractMatrix
    U_op = evolution_operator(Ham_0, t)
    return (U_op' * Ham_I * U_op).data
end

function get_interaction_ham_i_picture_a(
    Ham_int::AbstractMatrix,
    H_lambda::AbstractMatrix,
    H_S::AbstractMatrix,
    H_Sinv::AbstractMatrix,
    t::AbstractFloat,
)
    U_op = evolution_operator_a(H_lambda, H_S, H_Sinv, t)
    U_opd = adjoint(U_op)
    return U_opd * Ham_int * U_op
end
