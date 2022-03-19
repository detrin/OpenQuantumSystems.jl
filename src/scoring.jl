
########

function get_rmse_in_time(operator_vec1::OperatorVector, operator_vec2::OperatorVector)
    t_count = length(operator_vec1)
    c = 0
    for t_i in 1:t_count
        op1 = operator_vec1[t_i]
        op2 = operator_vec2[t_i]
        c += norm(op1.data - op2.data)^2
    end
    return sqrt(c / t_count)
end

function get_rmse_in_time(operator_vec1::OperatorVector, operator_data_vec2::Array)
    t_count = length(operator_vec1)
    c = 0
    for t_i in 1:t_count
        op1 = operator_vec1[t_i]
        c += norm(op1.data - operator_data_vec2[t_i, :, :])^2
    end
    return sqrt(c / t_count)
end

function get_rmse_in_time(operator_data_vec1::Array, operator_data_vec2::Array)
    t_count, M, K = size(operator_data_vec1)
    c = 0
    for t_i in 1:t_count
        c += norm(operator_data_vec1[t_i, :, :] - operator_data_vec2[t_i, :, :])^2
    end
    return sqrt(c / t_count)
end

#######

function compare_rho(rho::Array, rho_ref::Array; smooth_const=1e-9)
    N, M, K = size(rho)
    rho_sum = zeros(Float64, M, K)

    for t_i in 1:N
        rho_abs = abs.(rho_ref[t_i, :, :]) 
        rho_d = rho_abs + smooth_const*ones(size(rho_abs))
        rho_sum[:, :] += abs.(rho[t_i, :, :] - rho_ref[t_i, :, :]) ./ rho_d
    end
    rho_sum /= N
    return rho_sum
end

function compare_rho(rho::Array, rho_ref::OperatorVector; smooth_const=1e-9)
    rho_ref_array = operator_recast(rho_ref)
    return compare_rho(rho, rho_ref_array, smooth_const=smooth_const)
end

function compare_rho(rho::OperatorVector, rho_ref::Array; smooth_const=1e-9)
    rho_array = operator_recast(rho)
    return compare_rho(rho_array, rho_ref, smooth_const=smooth_const)
end

function compare_rho(rho::OperatorVector, rho_ref::OperatorVector; smooth_const=1e-9)
    rho_array = operator_recast(rho)
    rho_ref_array = operator_recast(rho_ref)
    return compare_rho(rho_array, rho_ref_array, smooth_const=smooth_const)
end

######

function compare_rho_in_time(rho::Array, rho_ref::Array; smooth_const=1e-9)
    N, M, K = size(rho)
    rho_rel = zeros(Float64, N, M, K)

    for t_i in 1:N
        rho_abs = abs.(rho_ref[t_i, :, :]) 
        rho_d = rho_abs + smooth_const*ones(size(rho_abs))
        rho_rel[t_i, :, :] = abs.(rho[t_i, :, :] - rho_ref[t_i, :, :]) ./ rho_d  
    end
    return rho_rel
end

function compare_rho_in_time(rho::Array, rho_ref::OperatorVector; smooth_const=1e-9)
    rho_ref_array = operator_recast(rho_ref)
    return compare_rho_in_time(rho, rho_ref_array; smooth_const=smooth_const)
end

function compare_rho_in_time(rho::OperatorVector, rho_ref::Array; smooth_const=1e-9)
    rho_array = operator_recast(rho)
    return compare_rho_in_time(rho_array, rho_ref; smooth_const=smooth_const)
end

function compare_rho_in_time(rho::OperatorVector, rho_ref::OperatorVector; smooth_const=1e-9)
    rho_ref_array = operator_recast(rho_ref)
    rho_array = operator_recast(rho)
    return compare_rho_in_time(rho_array, rho_ref_array; smooth_const=smooth_const)
end