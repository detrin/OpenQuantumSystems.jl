
const OperatorVectorType = Vector{Operator{GenericBasis{Vector{Int64}}, GenericBasis{Vector{Int64}}, Matrix{ComplexF64}}}
const OperatorDataVectorType = Array{ComplexF64, 3}

function get_rmse_in_time(operator_vec1::OperatorVectorType, operator_vec2::OperatorVectorType)
    t_count = length(operator_vec1)
    c = 0
    for t_i in 1:t_count
        op1 = operator_vec1[t_i]
        op2 = operator_vec2[t_i]
        c += norm(op1.data - op2.data)^2
    end
    return sqrt(c / t_count)
end

function get_rmse_in_time(operator_vec1::OperatorVectorType, operator_data_vec2::OperatorDataVectorType)
    t_count = length(operator_vec1)
    c = 0
    for t_i in 1:t_count
        op1 = operator_vec1[t_i]
        c += norm(op1.data - operator_data_vec2[t_i, :, :])^2
    end
    return sqrt(c / t_count)
end

function get_rmse_in_time(operator_data_vec1::OperatorDataVectorType, operator_data_vec2::OperatorDataVectorType)
    t_count = length(operator_vec1)
    c = 0
    for t_i in 1:t_count
        c += norm(operator_data_vec1[t_i, :, :] - operator_data_vec2[t_i, :, :])^2
    end
    return sqrt(c / t_count)
end