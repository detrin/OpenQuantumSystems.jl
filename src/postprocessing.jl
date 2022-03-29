
const OperatorVector = Vector{Operator{GenericBasis{Vector{Int64}}, GenericBasis{Vector{Int64}}, Matrix{ComplexF64}}}
const OperatorVectorArray = Array{ComplexF64, 3}

function operator_recast(rho::OperatorVector)::Array
    N = length(rho)
    M, K = size(rho[1].data)
    rho_array = zeros(eltype(rho[1].data), N, M, K)
    for t_i in 1:N
        rho_array[t_i, :, :] = rho[t_i].data
    end
    return rho_array
end

function operator_recast(rho::OperatorVectorArray)::Array
    return rho
end

######

function interaction_pic_to_schroedinger_pic(rho_int::OperatorVectorArray, tspan::Array, agg::Aggregate)
    rho_sch = deepcopy(rho_int)
    Ham_sys = agg.operators.Ham_sys
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham_sys, t)
        rho_sch[t_i, :, :] = U_op.data * rho_int[t_i, :, :] * (U_op').data 
    end
    return rho_sch
end

function interaction_pic_to_schroedinger_pic(rho_int::OperatorVector, tspan::Array, agg::Aggregate)
    rho_array = operator_recast(rho_int)
    return interaction_pic_to_schroedinger_pic(rho_array, tspan, agg)
end

function schroedinger_pic_to_interaction_pic(rho_sch::OperatorVectorArray, tspan::Array, agg::Aggregate)
    rho_int = deepcopy(rho_sch)
    Ham_sys = agg.operators.Ham_sys
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham_sys, t)
        rho_int[t_i, :, :] = (U_op').data * rho_sch[t_i, :, :] * U_op.data 
    end
    return rho_int
end

function schroedinger_pic_to_interaction_pic(rho_sch::OperatorVector, tspan::Array, agg::Aggregate)
    rho_array = operator_recast(rho_sch)
    return schroedinger_pic_to_interaction_pic(rho_array, tspan, agg)
end

#####

function local_st_to_exciton_st(rho_local::OperatorVectorArray, agg::Aggregate)
    N = size(rho_local, 1)
    rho_exciton = deepcopy(rho_local)
    Ham_sys = agg.operators.Ham_sys
    Ham_sys_lambda, Ham_sys_S = eigen(Ham_sys.data)
    Ham_sys_Sinv = inv(Ham_sys_S)
    # Ham_sys_lambda = diagm(Ham_sys_lambda)

    for t_i in 1:N
        rho_exciton[t_i, :, :] = Ham_sys_Sinv * rho_local[t_i, :, :] * Ham_sys_S
    end
    return rho_exciton
end

function local_st_to_exciton_st(rho_local::OperatorVector, agg::Aggregate)
    rho_array = operator_recast(rho_local)
    return local_st_to_exciton_st(rho_array, agg)
end

function exciton_st_to_local_st(rho_exciton::OperatorVectorArray, agg::Aggregate)
    N = size(rho_exciton, 1)
    rho_local = deepcopy(rho_exciton)
    Ham_sys = agg.operators.Ham_sys
    Ham_sys_lambda, Ham_sys_S = eigen(Ham_sys.data)
    Ham_sys_Sinv = inv(Ham_sys_S)
    # Ham_sys_lambda = diagm(Ham_sys_lambda)

    for t_i in 1:N
        rho_local[t_i, :, :] = Ham_sys_S * rho_exciton[t_i, :, :] * Ham_sys_Sinv
    end
    return rho_local
end

function exciton_st_to_local_st(rho_exciton::OperatorVector, agg::Aggregate)
    rho_array = operator_recast(rho_exciton)
    return exciton_st_to_local_st(rho_array, agg)
end

###

function tspan_cm_to_fs(tspan_cm)
    c_const = 299792458.0
    pi_const = 3.141592653589793
    c = 2.0*pi_const*c_const*1e-13
    tspan_fs = map((x) -> x / c, tspan_cm)
    return tspan_fs 
end