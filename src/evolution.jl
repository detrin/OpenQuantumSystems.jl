
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
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t_step = (t1 - t0) / (N-1)
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t0), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_op_array = [U_op]
    
    for t_i in 1:N-1
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
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t_step = (t1 - t0) / (N-1)
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t0), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_supop = spre(U_op') * spost(U_op)
    U_supop_array = [U_supop]
    
    for t_i in 1:N-1
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
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t_step = (t1 - t0) / (N-1)
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t0), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    cont(U_op) 

    for t_i in 1:N-1
        t = t0 + t_step * t_i
        U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
        U_op.data .= Ham_S * diagm(U_diagonal) * inv(Ham_S)
        cont(U_op) 
    end
end

EvolutionSuperOperatorIterator(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer) = @cont begin
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    t_step = (t1 - t0) / (N-1)
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * t0), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_supop = spre(U_op') * spost(U_op)
    cont(U_supop)

    for t_i in 1:N-1
        t = t0 + t_step * t_i
        U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
        U_op.data .= Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_supop.data .= (spre(U_op')).data * (spost(U_op)).data
        cont(U_supop)
    end
end

function EvolutionExact(ket0::Ket, Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer)
    ket_array = Array{Ket, 1}(undef, 0)
    foreach(EvolutionOperatorIterator(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer
        )) do U_op
        push!(ket_array, U_op * ket0)
    end
    return ket_array
end

function EvolutionExact(op0::Operator, Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer)
    op_array = Array{typeof(op0), 1}(undef, 0)
    foreach(EvolutionOperatorIterator(
        Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer
        )) do U_op
        push!(op_array, U_op' * op0 * U_op)
    end
    return op_array
end

function EvolutionApproximate(ket0::Ket, Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer)
    ket_array = [deepcopy(ket0)]
    t_step = (t1 - t0) / (N-1)
    U_op_step = EvolutionOperator(Hamiltonian, t_step)
    ket = deepcopy(ket0)
    for t_i in 1:N-1
        ket = U_op_step * ket
        push!(ket_array, ket)
    end
    return ket_array
end

function EvolutionApproximate(op0::Operator, Hamiltonian::Operator, t0::AbstractFloat, t1::AbstractFloat, N::Integer)
    op_array = [deepcopy(op0)]
    t_step = (t1 - t0) / (N-1)
    U_op_step = EvolutionOperator(Hamiltonian, t_step)
    U_op_step_d = U_op_step'
    op = deepcopy(op0)
    for t_i in 1:N-1
        op = U_op_step_d * op * U_op_step
        push!(op_array, op)
    end
    return op_array
end