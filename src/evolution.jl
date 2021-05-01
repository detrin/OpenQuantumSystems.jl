
using Continuables

function EvolutionOperator(Hamiltonian::Operator, t::AbstractFloat)::Operator
    return exp(-1im * Hamiltonian * t)
end

function EvolutionSuperOperator(Hamiltonian::Operator, t::AbstractFloat)::SuperOperator
    U = EvolutionOperator(Hamiltonian, t)
    return spre(U') * spost(U)
end

function EvolutionOperatorArray(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer)::Array
    Ham_lambda, Ham_S = eigen(Ham.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t = t0
    t_step = (t1 - t0) / N
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_op_array = [U_op]
    
    for t_i in 2:N
        t = t0 + t_step * t_i
        U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
        U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_op = DenseOperator(base, base, U)
        push!(U_op_array, U_op)
    end
    return U_op_array 
end

function EvolutionSuperOperatorArray(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer)::Array
    Ham_lambda, Ham_S = eigen(Ham.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t = t0
    t_step = (t1 - t0) / N
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_supop = spre(U_op') * spost(U_op)
    U_supop_array = [U_supop]
    
    for t_i in 2:N
        t = t0 + t_step * t_i
        U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
        U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_op = DenseOperator(base, base, U)
        U_supop = spre(U_op') * spost(U_op)
        push!(U_supop_array, U_supop)
    end
    return U_supop_array 
end

EvolutionOperatorIterator(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer) = @cont begin
    Ham_lambda, Ham_S = eigen(Ham.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t = t0
    t_step = (t1 - t0) / N
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    cont(U_op) 

    for t_i in 2:N
        t = t0 + t_step * t_i
        U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
        U_op.data .= Ham_S * diagm(U_diagonal) * inv(Ham_S)
        cont(U_op) 
    end

end

EvolutionSuperOperatorIterator(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer) = @cont begin
    Ham_lambda, Ham_S = eigen(Ham.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t = t0
    t_step = (t1 - t0) / N
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_supop = spre(U_op') * spost(U_op)
    cont(U_supop)

    for t_i in 2:N
        t = t0 + t_step * t_i
        U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
        U_op.data .= Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_supop.data .= spre(U_op') * spost(U_op)
        cont(U_supop)
    end
end