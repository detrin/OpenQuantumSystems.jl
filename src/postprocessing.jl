
const OperatorVector = Vector{Operator{GenericBasis{Vector{Int64}}, GenericBasis{Vector{Int64}}, Matrix{ComplexF64}}}
# const OperatorVectorArray = Array{ComplexF64, 3}

function operator_recast(rho::OperatorVector)::Array
    N = length(rho)
    M, K = size(rho[1].data)
    rho_array = zeros(eltype(rho[1].data), N, M, K)
    for t_i in 1:N
        rho_array[t_i, :, :] = rho[t_i].data
    end
    return rho_array
end

function operator_recast(rho::Array)::Array
    return rho
end


######

function interaction_pic_to_schroedinger_pic(rho_int::Array, tspan::Array, agg::Aggregate)
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
    rho_array = operator_recast(rho)
    return interaction_pic_to_schroedinger_pic(rho_array, tspan, agg)
end

function schroedinger_pic_to_interaction_pic(rho_sch::Array, tspan::Array, agg::Aggregate)
    rho_int = deepcopy(rho_sch)
    Ham_sys = agg.operators.Ham_sys
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham_sys, t)
        rho_int[t_i, :, :] = (U_op').data * rho_sch[t_i, :, :] * U_op.data 
    end
    return rho_int
end

function schroedinger_pic_to_interaction_pic(rho_int::OperatorVector, tspan::Array, agg::Aggregate)
    rho_array = operator_recast(rho)
    return schroedinger_pic_to_interaction_pic(rho_array, tspan, agg)
end

#####

function local_st_to_exciton_st(rho::Array)
end

function local_st_to_exciton_st(rho::OperatorVector)
end

