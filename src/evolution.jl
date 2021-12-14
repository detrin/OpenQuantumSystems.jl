
# using Continuables
# import ResumableFunctions: @resumable, @yield
using ResumableFunctions

function get_tspan(t0, t1, count)::Vector{Float64}
    t_step = (t1 - t0) / count
    tspan = [t0:t_step:t1;]
    return tspan
end

"""
    evolutionOperator(Hamiltonian, t)

Get evolution operator as Operator using the definition

``U(t) = e^{-i H t / \\hbar}, \\quad \\hbar = 1``.

"""
function evolutionOperator(Hamiltonian::Operator, t::AbstractFloat)::Operator
    return exp(-1im * Hamiltonian * t)
end

function evolutionOperator(Hamiltonian::Array, t::AbstractFloat)::Array
    return exp(-1im * Hamiltonian * t)
end

"""
    evolutionOperatorA(H_lambda, S, Sinv, t)

Get evolution operator as Array using the definition

``U(t) = S e^{-i H_\\lambda t / \\hbar} S^{-1},``

where ``\\hbar = 1``, ``H_{\\lambda,ii} = \\lambda_i`` are eigenvalues of ``H``,
so ``H`` has to be non-singular, otherwise ``H_{\\lambda,ij} = 0, i \\neq j``. 
``S`` is obtained from eigendecomposition of ``H``, for example

```julia
H_lambda, S = eigen(Hamiltonian.data)
Sinv = inv(S)
H_lambda = diagm(H_lambda)
```
and arguments have to be Arrays.

"""
function evolutionOperatorA(H_lambda::Array, S::Array, Sinv::Array, t::AbstractFloat)::Array
    return S * exp(-1im * H_lambda * t) * Sinv
end

"""
    evolutionSuperOperator(Hamiltonian, t)

Get evolution operator as SuperOperator using the definition

``\\mathcal{U}(t) \\:\\cdot\\: = U(t) \\:\\cdot\\: U^\\dagger(t) = e^{-i H t / \\hbar} \\:\\cdot\\: e^{i H t / \\hbar}, \\quad \\hbar = 1``.

"""
function evolutionSuperOperator(Hamiltonian::Operator, t::AbstractFloat)::SuperOperator
    U = evolutionOperator(Hamiltonian, t)
    return spre(U) * spost(U')
end

"""
    evolutionOperatorArray(Hamiltonian, tspan)

Get evolution operators as Vector or Operators using the definition

``U(t) = e^{-i H t / \\hbar}, \\quad \\hbar = 1``.

"""
function evolutionOperatorArray(Hamiltonian::Operator, tspan::Array)::Array
    N = length(tspan)
    Ham_lambda, S = eigen(Hamiltonian.data)
    Sinv = inv(S)
    U_diagonal = zeros(ComplexF64, size(Ham_lambda))

    U_ref = evolutionOperator(Hamiltonian, 0.0)
    U_op_array = [deepcopy(U_ref)]
    for t_i = 2:N
        push!(U_op_array, deepcopy(U_ref))
    end

    for t_i = 1:N
        U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
        U = S * diagm(U_diagonal) * inv(S)
        U_op_array[t_i].data[:, :] .= U
    end
    return U_op_array
end

"""
    evolutionOperatorArray(Hamiltonian, tspan)

Get evolution superoperators as Vector or SuperOperators using the definition

``\\mathcal{U}(t) \\:\\cdot\\: = U(t) \\:\\cdot\\: U^\\dagger(t) = e^{-i H t / \\hbar} \\:\\cdot\\: e^{i H t / \\hbar}, \\quad \\hbar = 1``.

"""
function evolutionSuperOperatorArray(Hamiltonian::Operator, tspan::Array)::Array
    N = length(tspan)
    Ham_lambda, S = eigen(Hamiltonian.data)
    basis = GenericBasis([size(Ham_lambda, 1)])
    Sinv = inv(S)
    U_diagonal = zeros(ComplexF64, size(Ham_lambda))

    UU_ref = evolutionSuperOperator(Hamiltonian, 0.0)
    U_supop_array = [deepcopy(UU_ref)]
    for t_i = 2:N
        push!(U_supop_array, deepcopy(UU_ref))
    end

    for t_i = 1:N
        U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
        U = S * diagm(U_diagonal) * inv(S)
        U_op = DenseOperator(basis, basis, U)
        U_supop_array[t_i] = spre(U_op) * spost(U_op')
    end
    return U_supop_array
end

"""
    evolutionOperatorIterator(Hamiltonian, tspan; diagonalize = true, approximate = false)

Resumable function that returns evolution operator as Operator type at the time t from tspan.
See [`evolutionOperator`](@ref). The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. The `approximate` option asusmes that tspan is made of 
equidistant points, therefore ``U(t)`` has to be calculated for ``t_0`` and ``t_\\text{step}``.

"""
function evolutionOperatorIterator end

@resumable function evolutionOperatorIterator(
    Hamiltonian::Operator,
    tspan::Array;
    diagonalize = true,
    approximate = false,
)
    N = length(tspan)
    t_step = tspan[2] - tspan[1]
    if diagonalize
        Ham_lambda, S = eigen(Hamiltonian.data)
        basis = GenericBasis([size(Ham_lambda, 1)])
        Sinv = inv(S)

        U_diagonal = zero(Ham_lambda)
        U_diagonal = map(lambda -> exp(-1im * lambda * tspan[1]), Ham_lambda)
        U = S * diagm(U_diagonal) * inv(S)
        U_op = DenseOperator(basis, basis, U)
    elseif !approximate
        U_op = evolutionOperator(Hamiltonian, tspan[1])
    else
        U_op = evolutionOperator(Hamiltonian, tspan[1])
        U_step = evolutionOperator(Hamiltonian, t_step)
    end
    @yield U_op

    for t_i = 2:N
        t = tspan[t_i]
        if diagonalize
            U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
            U_op.data .= S * diagm(U_diagonal) * inv(S)
        elseif !approximate
            U_op = evolutionOperator(Hamiltonian, t)
        else
            U_op = U_step * U_op
        end
        @yield U_op
    end
end


"""
    evolutionSuperOperatorIterator(Hamiltonian, tspan; 
    \tdiagonalize = true, approximate = false)

Resumable function that returns evolution operator as Operator type at the time t from tspan.
See [`evolutionSuperOperator`](@ref). The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. The `approximate` option asusmes that tspan is made of 
equidistant points, therefore ``U(t)`` has to be calculated for ``t_0`` and ``t_\\text{step}``.

"""
function evolutionSuperOperatorIterator end

@resumable function evolutionSuperOperatorIterator(
    Hamiltonian::Operator,
    tspan::Array;
    diagonalize = true,
    approximate = false,
)
    N = length(tspan)
    t_step = tspan[2] - tspan[1]
    if diagonalize
        Ham_lambda, S = eigen(Hamiltonian.data)
        basis = GenericBasis([size(Ham_lambda, 1)])
        Sinv = inv(S)

        U_diagonal = zero(Ham_lambda)
        U_diagonal = map(lambda -> exp(-1im * lambda * tspan[1]), Ham_lambda)
        U = S * diagm(U_diagonal) * Sinv
        U_op = DenseOperator(basis, basis, U)
    elseif !approximate
        U_op = evolutionOperator(Hamiltonian, tspan[1])
    else
        U_op = evolutionOperator(Hamiltonian, tspan[1])
        U_step = evolutionOperator(Hamiltonian, t_step)
    end
    U_supop = spre(U_op) * spost(U_op')
    @yield U_supop

    for t_i = 2:N
        t = tspan[t_i]
        if diagonalize
            U_diagonal .= map(lambda -> exp(-1im * lambda * tspan[t_i]), Ham_lambda)
            U_op.data .= S * diagm(U_diagonal) * Sinv
            U_supop.data .= (spre(U_op)).data * (spost(U_op')).data
        elseif !approximate
            U_op = evolutionOperator(Hamiltonian, t)
            U_supop = spre(U_op) * spost(U_op')
        else
            U_op = U_step * U_op
            U_supop = spre(U_op) * spost(U_op')
        end
        @yield U_supop
    end
end


"""
    evolutionExact(ket0, tspan, Hamiltonian; diagonalize = true, approximate = false)

Calculate exact time evolution of the `ket0` state see [`evolutionOperatorIterator`](@ref). 
The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. The `approximate` option asusmes that tspan is made of 
equidistant points, therefore ``U(t)`` has to be calculated for ``t_0`` and ``t_\\text{step}``.

"""
function evolutionExact(
    ket0::Ket,
    tspan::Array,
    Hamiltonian::Operator;
    diagonalize = true,
    approximate = false,
)
    N = length(tspan)
    ket_array = Array{Ket,1}(undef, 0)
    for t_i = 1:N
        push!(ket_array, deepcopy(ket0))
    end
    t_i = 0
    for U_op in evolutionOperatorIterator(
        Hamiltonian,
        tspan;
        diagonalize = diagonalize,
        approximate = approximate,
    )
        t_i += 1
        ket_array[t_i] = U_op * ket0
    end
    return ket_array
end


"""
    evolutionExact!(ket_array, ket0, tspan, Hamiltonian; 
    \tdiagonalize = true, approximate = false)

Calculate exact time evolution of the `ket0` state inplace see [`evolutionOperatorIterator`](@ref). 
The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. The `approximate` option asusmes that tspan is made of 
equidistant points, therefore ``U(t)`` has to be calculated for ``t_0`` and ``t_\\text{step}``.

"""
function evolutionExact!(
    ket_array::Array{Array{C,1},1},
    ket0::Ket,
    tspan::Array,
    Hamiltonian::Operator;
    diagonalize = true,
    approximate = false,
) where {C<:ComputableType}
    N = length(tspan)
    buffer = zeros(C, size(ket0.data))
    for t_i = 1:N
        push!(ket_array, deepcopy(ket0.data))
    end
    t_i = 0
    for U_op in evolutionOperatorIterator(
        Hamiltonian,
        tspan;
        diagonalize = diagonalize,
        approximate = approximate,
    )
        t_i += 1
        buffer[:] .= convert(Array{C,1}, (U_op * ket0).data)
        ket_array[t_i] = deepcopy(buffer)
    end
    return ket_array
end


"""
    evolutionExact(op0, tspan, Hamiltonian; diagonalize = true, approximate = false)

Calculate exact time evolution of the `op0` state see [`evolutionOperatorIterator`](@ref). 
The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. The `approximate` option asusmes that tspan is made of 
equidistant points, therefore ``U(t)`` has to be calculated for ``t_0`` and ``t_\\text{step}``.

"""
function evolutionExact(
    op0::Operator,
    tspan::Array,
    Hamiltonian::Operator;
    diagonalize = true,
    approximate = false,
)
    N = length(tspan)
    op_array = Array{typeof(op0),1}(undef, 0)
    for t_i = 1:N
        push!(op_array, deepcopy(op0))
    end
    t_i = 0
    for U_op in evolutionOperatorIterator(
        Hamiltonian,
        tspan;
        diagonalize = diagonalize,
        approximate = approximate,
    )
        t_i += 1
        op_array[t_i] = U_op * op0 * U_op'
    end
    return op_array
end


"""
    evolutionExact!(op_array, op0, tspan, Hamiltonian; 
    \tdiagonalize = true, approximate = false)

Calculate exact time evolution of the `op0` state inplace see [`evolutionOperatorIterator`](@ref). 
The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. The `approximate` option asusmes that tspan is made of 
equidistant points, therefore ``U(t)`` has to be calculated for ``t_0`` and ``t_\\text{step}``.

"""
function evolutionExact!(
    op_array::Array{Array{C,2},1},
    op0::Operator,
    tspan::Array,
    Hamiltonian::Operator,
) where {C<:ComputableType}
    N = length(tspan)
    buffer = zeros(C, size(op0.data))
    for t_i = 1:N
        push!(op_array, deepcopy(op0.data))
    end
    t_i = 0
    for U_op in evolutionOperatorIterator(Hamiltonian, tspan)
        t_i += 1
        buffer[:, :] .= convert(Array{C,2}, (U_op * op0 * U_op').data)
        op_array[t_i] = deepcopy(buffer)
        push!(op_array)
    end
    return op_array
end


"""
    evolutionApproximate(ket0, tspan, Hamiltonian)

Calculate approximate time evolution of the `ket0` based on ``U(t)`` that is 
calculated for ``t_0`` and ``t_\\text{step}``. See [`evolutionOperator`](@ref).
This method returns Vector of Ket states.

"""
function evolutionApproximate(ket0::Ket, tspan::Array, Hamiltonian::Operator)
    N = length(tspan)
    ket_array = Array{typeof(ket0),1}(undef, 0)
    t_step = tspan[2] - tspan[1]
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    U_op0 = evolutionOperator(Hamiltonian, tspan[1])
    ket = U_op0 * deepcopy(ket0)
    for t_i = 1:N
        push!(ket_array, deepcopy(ket))
    end
    for t_i = 2:N
        ket = U_op_step * ket
        ket_array[t_i] = deepcopy(ket)
    end
    return ket_array
end


"""
    evolutionApproximate!(ket_array, ket0, tspan, Hamiltonian)

Calculate approximate time evolution of the `ket0` inplace based on ``U(t)`` that is 
calculated for ``t_0`` and ``t_\\text{step}``. See [`evolutionOperator`](@ref).
Argument `ket_array` is Vector of Arrays.

"""
function evolutionApproximate!(
    ket_array::Array{Array{C,1},1},
    ket0::Ket,
    tspan::Array,
    Hamiltonian::Operator,
) where {C<:ComputableType}
    N = length(tspan)
    t_step = tspan[2] - tspan[1]
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    U_op0 = evolutionOperator(Hamiltonian, tspan[1])
    ket = U_op0 * deepcopy(ket0)
    buffer = zeros(C, size(ket.data))
    buffer[:] .= convert(Array{C,1}, ket.data)
    for t_i = 1:N
        push!(ket_array, deepcopy(buffer))
    end
    for t_i = 1:N-1
        ket = U_op_step * ket
        buffer[:] .= convert(Array{C,1}, ket.data)
        ket_array[t_i+1] = deepcopy(buffer)
    end
    return ket_array
end


"""
    evolutionApproximate(op0, tspan, Hamiltonian)

Calculate approximate time evolution of the `op0` based on ``U(t)`` that is 
calculated for ``t_0`` and ``t_\\text{step}``. See [`evolutionOperator`](@ref).
This method returns Vector of Operators.

"""
function evolutionApproximate(op0::Operator, tspan::Array, Hamiltonian::Operator)
    N = length(tspan)
    op_array = Array{typeof(op0),1}(undef, 0)
    t_step = tspan[2] - tspan[1]
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    U_op0 = evolutionOperator(Hamiltonian, tspan[1])
    U_op_step_d = U_op_step'
    op = U_op0 * deepcopy(op0) * U_op0'
    for t_i = 1:N
        push!(op_array, deepcopy(op))
    end
    for t_i = 2:N
        op = U_op_step * op * U_op_step_d
        op_array[t_i] = deepcopy(op)
    end
    return op_array
end


"""
    evolutionApproximate!(op_array, op0, tspan, Hamiltonian)

Calculate approximate time evolution of the `op0` inplace based on ``U(t)`` that is 
calculated for ``t_0`` and ``t_\\text{step}``. See [`evolutionOperator`](@ref).
This method returns Vector of Arrays.

"""
function evolutionApproximate!(
    op_array::Array{Array{C,2},1},
    op0::Operator,
    tspan::Array,
    Hamiltonian::Operator,
) where {C<:ComputableType}
    N = length(tspan)
    t_step = tspan[2] - tspan[1]
    U_op_step = evolutionOperator(Hamiltonian, t_step)
    U_op0 = evolutionOperator(Hamiltonian, tspan[1])
    U_op_step_d = U_op_step'
    op = U_op0 * deepcopy(op0) * U_op0'
    buffer = zeros(C, size(op.data))
    buffer[:, :] .= convert(Array{C,2}, op.data)
    for t_i = 1:N
        push!(op_array, deepcopy(buffer))
    end
    for t_i = 2:N
        op = U_op_step * op * U_op_step_d
        buffer[:, :] .= convert(Array{C,2}, op.data)
        op_array[t_i] = deepcopy(buffer)
    end
    return op_array
end

"""
    evolution_exact(rho0, tspan, Hamiltonian; diagonalize = false)

Calculate exact time evolution of the `rho0` inplace based on ``U(t)``. 
The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. See [`evolutionOperator`](@ref).
This method returns tspan and Vector of Operators.

"""
function evolution_exact(
    rho0::T,
    tspan::Array,
    Hamiltonian::U;
    diagonalize = false,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B}}
    N = length(tspan)

    rho = deepcopy(rho0)
    if diagonalize
        Ham_lambda, S = eigen(Hamiltonian.data)
        basis = GenericBasis([size(Ham_lambda, 1)])
        Sinv = inv(S)
        U_diagonal = zero(Ham_lambda)
        U_diagonal = map(lambda -> exp(-1im * lambda * tspan[1]), Ham_lambda)
        U_data = S * diagm(U_diagonal) * Sinv
        U_op = DenseOperator(basis, basis, U_data)
    else
        U_op = evolutionOperator(Hamiltonian, tspan[1])
    end
    rho = U_op * rho0 * U_op'
    rho_t = [rho]

    for t_i = 2:N
        t = tspan[t_i]
        if diagonalize
            U_diagonal .= map(lambda -> exp(-1im * lambda * t), Ham_lambda)
            U_op.data .= S * diagm(U_diagonal) * Sinv
        else
            U_op = evolutionOperator(Hamiltonian, t)
        end
        rho = U_op * rho0 * U_op'
        push!(rho_t, rho)
    end
    return tspan, rho_t
end

"""
    evolution_exact(rho0, tspan, Hamiltonian; diagonalize = false)

Calculate approximate time evolution of the `rho0` based on ``U(t)`` that is 
calculated for ``t_0`` and ``t_\\text{step}``. The `diagonalize` argument decompose Hamiltonian into ``\\lambda_i, S, S^{-1}``
and calculate the exponential using eigenvalue decomposition. See [`evolutionOperator`](@ref).
This method returns tspan and Vector of Operators.

"""
function evolution_approximate(
    rho0::T,
    tspan::Array,
    Ham::U;
    diagonalize = false,
) where {B<:Basis,T<:Operator{B,B},U<:Operator{B,B}}
    N = length(tspan)
    t_step = tspan[2] - tspan[1]
    for t_i = 1:N-1
        if !(tspan[t_i+1] - tspan[t_i] ≈ t_step)
            throw(ErrorException("Points in in tspan have to be equidistant."))
        end
    end

    rho = deepcopy(rho0)
    if diagonalize
        Ham_lambda, S = eigen(Ham.data)
        basis = GenericBasis([size(Ham_lambda, 1)])
        Sinv = inv(S)

        U_diagonal = zero(Ham_lambda)
        U_diagonal = map(lambda -> exp(-1im * lambda * t_step), Ham_lambda)
        U_data = S * diagm(U_diagonal) * Sinv
        U_op_step = DenseOperator(basis, basis, U_data)

        U_diagonal = map(lambda -> exp(-1im * lambda * tspan[1]), Ham_lambda)
        U_data = S * diagm(U_diagonal) * Sinv
        U_op0 = DenseOperator(basis, basis, U_data)
    else
        U_op_step = evolutionOperator(Ham, t_step)
        U_op0 = evolutionOperator(Ham, tspan[1])
    end
    U_op_step_d = U_op_step'
    U_op0_d = U_op0'
    rho = U_op0 * rho0 * U_op0_d
    rho_t = [rho]

    for t_i = 2:N
        t = tspan[t_i]
        rho = U_op_step * rho * U_op_step_d
        push!(rho_t, rho)
    end
    return tspan, rho_t
end

function evolutionOperatorExp(
    Ham::Array, t::AbstractFloat, n::Integer
)::Array
    # c = ones(size(Ham))
    c = one(Ham)
    Len = size(Ham, 1)
    s = c 
    if n > 0
        for k in 1:n
            c = c * (-1.0im/k * Ham * t)
            s = s + c
            # s[:, :] = s + (-1.0im * Ham * t)^k/factorial(big(k))
        end
    end
    s /= abs(det(s))^(1.0/Len)
    s
end

function evolutionOperatorExp(
    Ham::T, t::AbstractFloat, n::Integer
) where {B<:Basis,T<:Operator{B,B}}
    data = evolutionOperatorExp(Ham.data, t, n)
    DenseOperator(Ham.basis_l, Ham.basis_r, data)
end

function evolution_el_part(
    Ham::Array, 
    t::AbstractFloat, 
    a::Integer, 
    b::Integer, 
    indicesMap::IndicesMap
)
    data = take_el_part(Ham, a, b, indicesMap)
    return evolutionOperator(data, t)
end

function evolution_el_part(
    Ham::T, 
    t::AbstractFloat, 
    a::Integer, 
    b::Integer, 
    indicesMap::IndicesMap
) where {B<:Basis,T<:Operator{B,B}}
    data = evolution_el_part(Ham.data, t, a, b, indicesMap)
    b = GenericBasis([size(data, 1)])
    return DenseOperator(b, b, data)
end

function Evolution_SI_exact(
    W0::T,
    tspan::Array,
    agg::Aggregate
) where {B<:Basis,T<:Operator{B,B}}
    W_t_exact = zeros(ComplexF64, length(tspan), agg.tools.bSize, agg.tools.bSize)
    Ham = agg.operators.Ham
    Ham_0 = agg.operators.Ham_0
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        U_0_op = evolutionOperator(Ham_0, t)
        W = U_0_op' * W * U_0_op
        W_t_exact[t_i, :, :] = W.data
    end
    return tspan, W_t_exact
end

function Evolution_sI_exact(
    W0::T,
    tspan::Array,
    agg::Aggregate
) where {B<:Basis,T<:Operator{B,B}}
    elLen = agg.core.molCount
    rho_t_exact = zeros(ComplexF64, length(tspan), elLen+1, elLen+1)

    Ham = agg.operators.Ham
    Ham_0 = agg.operators.Ham_0
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        U_0_op = evolutionOperator(Ham_0, t)
        W = U_0_op' * W * U_0_op
        rho_t_exact[t_i, :, :] = trace_bath(W.data[:, :], agg.core, agg.tools)
    end
    return tspan, rho_t_exact
end

function Evolution_SS_exact(
    W0::T,
    tspan::Array,
    agg::Aggregate
) where {B<:Basis,T<:Operator{B,B}}
    W_t_exact = zeros(ComplexF64, length(tspan), agg.tools.bSize, agg.tools.bSize)
    Ham = agg.operators.Ham
    Ham_0 = agg.operators.Ham_0
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        W_t_exact[t_i, :, :] = W.data
    end
    return tspan, W_t_exact
end

function Evolution_sS_exact(
    W0::T,
    tspan::Array,
    agg::Aggregate
) where {B<:Basis,T<:Operator{B,B}}
    elLen = agg.core.molCount
    rho_t_exact = zeros(ComplexF64, length(tspan), elLen+1, elLen+1)

    Ham = agg.operators.Ham
    Ham_0 = agg.operators.Ham_0
    for t_i in 1:length(tspan)
        t = tspan[t_i]
        U_op = evolutionOperator(Ham, t)
        W = U_op * W0 * U_op'
        rho_t_exact[t_i, :, :] = trace_bath(W.data[:, :], agg.core, agg.tools)
    end
    return tspan, rho_t_exact
end