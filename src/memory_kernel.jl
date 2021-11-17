
"""
    take_el_part(A, a, b, vibindices)

Take electric part specified by electric indices `a` and `b` from the A (type of Array). 

"""
function take_el_part(A::Array, a, b, vibindices)
    a1 = vibindices[a][1]
    a2 = vibindices[a][end]
    b1 = vibindices[b][1]
    b2 = vibindices[b][end]

    return A[a1:a2, b1:b2]
end



"""
    MemoryKernel_1_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, vibindices)

Calculate the first part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(t) \\hat{H}_I^{(I)}(\\tau) W_\\text{bath} \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `vibindices`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_1_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    agg,
    FCProd,
    aggIndices,
    vibindices
)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    elLen += 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_t_tau = H_II_t * H_II_tau

    for a = 1:elLen
        for c = 1:elLen
            H_II_t_tau_ac = take_el_part(H_II_t_tau, a, c, vibindices)
            for d = 1:elLen
                W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                MK_big = H_II_t_tau_ac * W_bath_cd
                MemoryKernel[a, d, c, d] =
                    trace_bath_part(MK_big, a, d, agg, FCProd, aggIndices, vibindices)
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_2_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, vibindices)

Calculate the second part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(t) W_\\text{bath} \\hat{H}_I^{(I)}(\\tau) \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `vibindices`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_2_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    agg,
    FCProd,
    aggIndices,
    vibindices
)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    elLen += 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    for a = 1:elLen
        for b = 1:elLen
            for c = 1:elLen
                H_II_t_ac = take_el_part(H_II_t, a, c, vibindices)
                for d = 1:elLen
                    H_II_tau_db = take_el_part(H_II_tau, d, b, vibindices)
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    MK_big = H_II_t_ac * W_bath_cd * H_II_tau_db
                    MemoryKernel[a, b, c, d] = trace_bath_part(
                        MK_big,
                        a,
                        b,
                        agg,
                        FCProd,
                        aggIndices,
                        vibindices,
                    )
                end
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_3_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, vibindices)

Calculate the third part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(\\tau) W_\\text{bath} \\hat{H}_I^{(I)}(t) \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `vibindices`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_3_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    agg,
    FCProd,
    aggIndices,
    vibindices
)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    elLen += 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    for a = 1:elLen
        for b = 1:elLen
            for c = 1:elLen
                H_II_tau_ac = take_el_part(H_II_tau, a, c, vibindices)
                for d = 1:elLen
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    H_II_t_db = take_el_part(H_II_t, d, b, vibindices)
                    MK_big = H_II_tau_ac * W_bath_cd * H_II_t_db
                    MemoryKernel[a, b, c, d] = trace_bath_part(
                        MK_big,
                        a,
                        b,
                        agg,
                        FCProd,
                        aggIndices,
                        vibindices,
                    )
                end
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_4_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, vibindices)

Calculate the fourth part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ W_\\text{bath} \\hat{H}_I^{(I)}(t) \\hat{H}_I^{(I)}(\\tau) \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `vibindices`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_4_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    agg,
    FCProd,
    aggIndices,
    vibindices
)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    elLen += 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_tau_t = H_II_tau * H_II_t

    for a = 1:elLen
        for b = 1:elLen
            for d = 1:elLen
                W_bath_ad = take_el_part(W_bath, a, d, vibindices)
                H_II_tau_t_db = take_el_part(H_II_tau_t, d, b, vibindices)
                MK_big = W_bath_ad * H_II_tau_t_db
                MemoryKernel[a, b, a, d] =
                    trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, vibindices)

Calculate Memory Kernel with the definition

`` \\mathcal{M}(t, \\tau) = \\operatorname{tr}_B \\{ [ \\hat{H}_I^{(I)}(t), [ \\hat{H}_I^{(I)}(\\tau), W_\\text{bath} ]]\\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `vibindices`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    agg,
    FCProd,
    aggIndices,
    vibindices
)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    elLen += 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_tau_t = H_II_tau * H_II_t
    H_II_t_tau = H_II_t * H_II_tau

    for a = 1:elLen
        for c = 1:elLen
            H_II_t_tau_ac = take_el_part(H_II_t_tau, a, c, vibindices)
            H_II_tau_ac = take_el_part(H_II_tau, a, c, vibindices)
            H_II_t_ac = take_el_part(H_II_t, a, c, vibindices)
            for d = 1:elLen
                W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                for b = 1:elLen
                    H_II_tau_db = take_el_part(H_II_tau, d, b, vibindices)

                    MK_big = -H_II_t_ac * W_bath_cd * H_II_tau_db

                    H_II_t_db = take_el_part(H_II_t, d, b, vibindices)
                    MK_big[:, :] -= H_II_tau_ac * W_bath_cd * H_II_t_db

                    if a == c
                        W_bath_ad = take_el_part(W_bath, a, d, vibindices)
                        H_II_tau_t_db = take_el_part(H_II_tau_t, d, b, vibindices)
                        MK_big[:, :] += W_bath_ad * H_II_tau_t_db
                    end

                    if d == b
                        MK_big[:, :] += H_II_t_tau_ac * W_bath_cd
                    end
                    MemoryKernel[a, b, c, d] = trace_bath_part(
                        MK_big,
                        a,
                        b,
                        agg,
                        FCProd,
                        aggIndices,
                        vibindices,
                    )
                end
            end
        end
    end
    return MemoryKernel
end
