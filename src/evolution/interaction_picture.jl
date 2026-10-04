
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

# Same as `get_interaction_ham_i_picture`, but reuses a cached eigendecomposition
# of H_0: U_0' H_I U_0 = V (e^{iλt} ∘ V' H_I V ∘ e^{-iλt}) V'.
function _interaction_ham(Ham_0::AbstractMatrix, Ham_I::AbstractMatrix, t::Real)::AbstractMatrix
    ip = _interaction_picture(Ham_0, Ham_I)
    if ip === nothing
        U_array = evolution_operator(Ham_0, t)
        return adjoint(U_array) * Ham_I * U_array
    end
    phase = cis.(ip.eig.values .* t)
    V = ip.eig.vectors
    return V * (phase .* ip.Ham_I_eig .* transpose(conj.(phase))) * V'
end

function _interaction_ham(Ham_0::Operator, Ham_I::Operator, t::Real)::Operator
    H = _interaction_ham(Ham_0.data, Ham_I.data, t)
    return DenseOperator(Ham_I.basis_l, Ham_I.basis_r, H)
end

# H_I in the eigenbasis of H_0, cached next to the eigendecomposition of H_0.
struct _InteractionPicture
    eig::_HermitianEigen
    Ham_I::Matrix
    Ham_I_eig::Matrix{ComplexF64}
end

const _INTERACTION_CACHE = Dict{Tuple{UInt,UInt},_InteractionPicture}()
const _INTERACTION_CACHE_LOCK = ReentrantLock()

function _interaction_picture(Ham_0::AbstractMatrix, Ham_I::AbstractMatrix)
    eig = _hermitian_eigen(Ham_0)
    eig === nothing && return nothing
    key = (objectid(Ham_0), objectid(Ham_I))
    lock(_INTERACTION_CACHE_LOCK) do
        entry = get(_INTERACTION_CACHE, key, nothing)
        entry !== nothing && entry.eig === eig && entry.Ham_I == Ham_I && return entry
        V = eig.vectors
        entry = _InteractionPicture(eig, copy(Ham_I), V' * Ham_I * V)
        length(_INTERACTION_CACHE) >= _EIGEN_CACHE_SIZE && empty!(_INTERACTION_CACHE)
        _INTERACTION_CACHE[key] = entry
        return entry
    end
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
