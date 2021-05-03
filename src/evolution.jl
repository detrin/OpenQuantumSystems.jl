
using Continuables

function evolutionOperator(Hamiltonian::Operator, t::AbstractFloat)::Operator
    return exp(-1im * Hamiltonian * t)
end

function evolutionSuperOperator(Hamiltonian::Operator, t::AbstractFloat)::SuperOperator
    U = evolutionOperator(Hamiltonian, t)
    return spre(U') * spost(U)
end

function evolutionOperatorArray(Hamiltonian::Operator, tspan::Array)::Array
    N = length(tspan)
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    Ham_Sinv = inv(Ham_S)
    U_diagonal = zeros(ComplexF64, size(Ham_lambda))
    
    U_ref = evolutionOperator(Hamiltonian, 0.0)
    U_op_array = [deepcopy(U_ref)]
    for t_i in 2:N
        push!(U_op_array, deepcopy(U_ref))
    end

    for t_i in 1:N
        U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
        U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_op_array[t_i].data[:, :] .= U
    end
    return U_op_array 
end

function evolutionSuperOperatorArray(Hamiltonian::Operator, tspan::Array)::Array
    N = length(tspan)
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    U_diagonal = zeros(ComplexF64, size(Ham_lambda))

    UU_ref = evolutionSuperOperator(Hamiltonian, 0.0)
    U_supop_array = [deepcopy(UU_ref)]
    for t_i in 2:N
        push!(U_supop_array, deepcopy(UU_ref))
    end
    
    for t_i in 1:N
        U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
        U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_op = DenseOperator(base, base, U)
        U_supop = spre(U_op') * spost(U_op)
        U_supop_array[t_i] = spre(U_op') * spost(U_op)
    end
    return U_supop_array 
end

evolutionOperatorIterator(Hamiltonian::Operator, tspan::Array) = @cont begin
    N = length(tspan)
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * tspan[1]), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    cont(U_op) 

    for t_i in 2:N
        U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
        U_op.data .= Ham_S * diagm(U_diagonal) * inv(Ham_S)
        cont(U_op) 
    end
end

evolutionSuperOperatorIterator(Hamiltonian::Operator, tspan::Array) = @cont begin
    N = length(tspan)
    Ham_lambda, Ham_S = eigen(Hamiltonian.data)
    base = GenericBasis([size(Ham_lambda, 1)])
    Ham_Sinv = inv(Ham_S)
    
    U_diagonal = zero(Ham_lambda)
    U_diagonal = map(lambda -> exp(-1im * lambda * tspan[1]), Ham_lambda)
    U = Ham_S * diagm(U_diagonal) * inv(Ham_S)
    U_op = DenseOperator(base, base, U)
    U_supop = spre(U_op') * spost(U_op)
    cont(U_supop)

    for t_i in 2:N
        U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
        U_op.data .= Ham_S * diagm(U_diagonal) * inv(Ham_S)
        U_supop.data .= (spre(U_op')).data * (spost(U_op)).data
        cont(U_supop)
    end
end

function evolutionExact(ket0::Ket, Hamiltonian::Operator, tspan::Array)
    N = length(tspan)
    ket_array = Array{Ket, 1}(undef, 0)
    for t_i in 1:N
        push!(ket_array, deepcopy(ket0))
    end
    t_i = 0
    foreach(evolutionOperatorIterator(Hamiltonian, tspan)) do U_op
        t_i += 1
        ket_array[t_i] = U_op * ket0
    end
    return ket_array
end

function evolutionExact!(ket_array::Array{Array{C,1},1}, ket0::Ket, Hamiltonian::Operator, tspan::Array) where C<:ComputableType
    N = length(tspan)
    buffer = zeros(C, size(ket0.data))
    for t_i in 1:N
        push!(ket_array, deepcopy(ket0.data))
    end
    t_i = 0
    foreach(evolutionOperatorIterator(Hamiltonian, tspan)) do U_op
        t_i += 1
        buffer[:] .= convert(Array{C,1}, (U_op * ket0).data)
        ket_array[t_i] = deepcopy(buffer)
    end
    return ket_array
end

function evolutionExact(op0::Operator, Hamiltonian::Operator, tspan::Array)
    N = length(tspan)
    op_array = Array{typeof(op0), 1}(undef, 0)
    for t_i in 1:N
        push!(op_array, deepcopy(op0))
    end
    t_i = 0
    foreach(evolutionOperatorIterator(Hamiltonian, tspan)) do U_op
        t_i += 1
        op_array[t_i] = U_op' * op0 * U_op
    end
    return op_array
end

function evolutionExact!(op_array::Array{Array{C,2},1}, op0::Operator, Hamiltonian::Operator, tspan::Array) where C<:ComputableType
    N = length(tspan)
    buffer = zeros(C, size(op0.data))
    for t_i in 1:N
        push!(op_array, deepcopy(op0.data))
    end
    t_i = 0
    foreach(evolutionOperatorIterator(Hamiltonian, tspan)) do U_op
        t_i += 1
        buffer[:,:] .= convert(Array{C,2}, (U_op' * op0 * U_op).data)
        op_array[t_i] = deepcopy(buffer)
        push!(op_array, )
    end
    return op_array
end

function evolutionApproximate(ket0::Ket, Hamiltonian::Operator, tspan::Array)
    N = length(tspan)
    ket_array = Array{typeof(ket0), 1}(undef, 0)
    for t_i in 1:N
        push!(ket_array, deepcopy(ket0))
    end
    t_step = tspan[2] - tspan[1] 
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    ket = deepcopy(ket0)
    for t_i in 2:N
        ket = U_op_step * ket
        ket_array[t_i] = deepcopy(ket)
    end
    return ket_array
end

function evolutionApproximate!(ket_array::Array{Array{C,1},1}, ket0::Ket, Hamiltonian::Operator, tspan::Array) where C<:ComputableType
    N = length(tspan)
    buffer = zeros(C, size(ket0.data))
    buffer[:] .= convert(Array{C,1}, ket0.data)
    for t_i in 1:N
        push!(ket_array, deepcopy(buffer))
    end
    t_step = tspan[2] - tspan[1] 
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    ket = deepcopy(ket0)
    for t_i in 1:N-1
        ket = U_op_step * ket
        buffer[:] .= convert(Array{C,1}, ket.data)
        ket_array[t_i+1] = deepcopy(buffer)
    end
    return ket_array
end

function evolutionApproximate(op0::Operator, Hamiltonian::Operator, tspan::Array)
    N = length(tspan)
    op_array = Array{typeof(op0), 1}(undef, 0)
    for t_i in 1:N
        push!(op_array, deepcopy(op0))
    end
    t_step = tspan[2] - tspan[1] 
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    U_op_step_d = U_op_step'
    op = deepcopy(op0)
    for t_i in 2:N
        op = U_op_step_d * op * U_op_step
        op_array[t_i] = deepcopy(op)
    end
    return op_array
end

function evolutionApproximate!(op_array::Array{Array{C,2},1}, op0::Operator, Hamiltonian::Operator, tspan::Array) where C<:ComputableType
    N = length(tspan)
    buffer = zeros(C, size(op0.data))
    buffer[:,:] .= convert(Array{C,2}, op0.data)
    for t_i in 1:N
        push!(op_array, deepcopy(buffer))
    end
    t_step = tspan[2] - tspan[1]  
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    U_op_step_d = U_op_step'
    op = deepcopy(op0)
    for t_i in 2:N
        op = U_op_step_d * op * U_op_step
        buffer[:,:] .= convert(Array{C,2}, op.data)
        op_array[t_i] = deepcopy(buffer)
    end
    return op_array
end