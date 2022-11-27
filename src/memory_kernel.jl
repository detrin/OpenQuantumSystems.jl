


"""
    MemoryKernel_1_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, indicesMap)

Calculate the first part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(t) \\hat{H}_I^{(I)}(\\tau) W_\\text{bath} \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `indicesMap`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_1_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_t_tau = H_II_t * H_II_tau

    for a = 1:elLen
        for c = 1:elLen
            H_II_t_tau_ac = take_el_part(H_II_t_tau, a, c, indicesMap)
            for d = 1:elLen
                W_bath_cd = take_el_part(W_bath, c, d, indicesMap)
                MK_big = H_II_t_tau_ac * W_bath_cd
                MemoryKernel[a, d, c, d] =
                    trace_bath_part(MK_big, a, d, aggTools; vib_basis=aggOperators.vib_basis)
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_2_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, indicesMap)

Calculate the second part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(t) W_\\text{bath} \\hat{H}_I^{(I)}(\\tau) \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `indicesMap`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_2_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    for a = 1:elLen
        for b = 1:elLen
            for c = 1:elLen
                H_II_t_ac = take_el_part(H_II_t, a, c, indicesMap)
                for d = 1:elLen
                    H_II_tau_db = take_el_part(H_II_tau, d, b, indicesMap)
                    W_bath_cd = take_el_part(W_bath, c, d, indicesMap)
                    MK_big = H_II_t_ac * W_bath_cd * H_II_tau_db
                    MemoryKernel[a, b, c, d] = trace_bath_part(
                        MK_big,
                        a,
                        b,
                        aggTools;
                        vib_basis=aggOperators.vib_basis
                    )
                end
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_3_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, indicesMap)

Calculate the third part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ \\hat{H}_I^{(I)}(\\tau) W_\\text{bath} \\hat{H}_I^{(I)}(t) \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `indicesMap`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_3_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    for a = 1:elLen
        for b = 1:elLen
            for c = 1:elLen
                H_II_tau_ac = take_el_part(H_II_tau, a, c, indicesMap)
                for d = 1:elLen
                    W_bath_cd = take_el_part(W_bath, c, d, indicesMap)
                    H_II_t_db = take_el_part(H_II_t, d, b, indicesMap)
                    MK_big = H_II_tau_ac * W_bath_cd * H_II_t_db
                    MemoryKernel[a, b, c, d] = trace_bath_part(
                        MK_big,
                        a,
                        b,
                        aggTools;
                        vib_basis=aggOperators.vib_basis
                    )
                end
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_4_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, indicesMap)

Calculate the fourth part of Memory Kernel with the definition

`` \\mathcal{M}_1(t, \\tau) = \\operatorname{tr}_B \\{ W_\\text{bath} \\hat{H}_I^{(I)}(t) \\hat{H}_I^{(I)}(\\tau) \\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `indicesMap`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_4_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_tau_t = H_II_tau * H_II_t

    for a = 1:elLen
        for b = 1:elLen
            for d = 1:elLen
                W_bath_ad = take_el_part(W_bath, a, d, indicesMap)
                H_II_tau_t_db = take_el_part(H_II_tau_t, d, b, indicesMap)
                MK_big = W_bath_ad * H_II_tau_t_db
                MemoryKernel[a, b, a, d] =
                    trace_bath_part(MK_big, a, b, aggTools; vib_basis=aggOperators.vib_basis)
            end
        end
    end
    return MemoryKernel
end


"""
    MemoryKernel_traced(H_II_t, H_II_tau, W_bath, agg, FCProd, aggIndices, indicesMap)

Calculate Memory Kernel with the definition

`` \\mathcal{M}(t, \\tau) = \\operatorname{tr}_B \\{ [ \\hat{H}_I^{(I)}(t), [ \\hat{H}_I^{(I)}(\\tau), W_\\text{bath} ]]\\}``.

# Arguments
* `H_II_t`: Interaction Hamiltonian in interaction picutre at the time t, ``\\hat{H}_I^{(I)}(t)``.
* `H_II_tau`: Interaction Hamiltonian in interaction picutre at the time `tau`, ``\\hat{H}_I^{(I)}(\\tau)``.
* `W_bath`: Density matrix representing bath part of the density matrix, see [`get_rho_bath`](@ref).
* `agg`: Aggregate of molecules, see [`Aggregate`](@ref).
* `aggIndices`: Aggregate indices, see [`getIndices`](@ref).
* `indicesMap`: Aggregate vibrational indices, see [`getVibIndices`](@ref).
"""
function MemoryKernel_traced(
    H_II_t::Array,
    H_II_tau::Array,
    W_bath::Array,
    aggCore::AggregateCore,
    aggOperators::AggregateOperators,
    aggTools::AggregateTools
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_tau_t = H_II_tau * H_II_t
    H_II_t_tau = H_II_t * H_II_tau

    for a = 1:elLen
        for c = 1:elLen
            H_II_t_tau_ac = take_el_part(H_II_t_tau, a, c, indicesMap)
            H_II_tau_ac = take_el_part(H_II_tau, a, c, indicesMap)
            H_II_t_ac = take_el_part(H_II_t, a, c, indicesMap)
            for d = 1:elLen
                W_bath_cd = take_el_part(W_bath, c, d, indicesMap)
                W_bath_cd = take_el_part(W_bath, c, d, indicesMap)
                for b = 1:elLen
                    H_II_tau_db = take_el_part(H_II_tau, d, b, indicesMap)

                    MK_big = -H_II_t_ac * W_bath_cd * H_II_tau_db

                    H_II_t_db = take_el_part(H_II_t, d, b, indicesMap)
                    MK_big[:, :] -= H_II_tau_ac * W_bath_cd * H_II_t_db

                    if a == c
                        W_bath_ad = take_el_part(W_bath, a, d, indicesMap)
                        H_II_tau_t_db = take_el_part(H_II_tau_t, d, b, indicesMap)
                        MK_big[:, :] += W_bath_ad * H_II_tau_t_db
                    end

                    if d == b
                        MK_big[:, :] += H_II_t_tau_ac * W_bath_cd
                    end
                    MemoryKernel[a, b, c, d] = trace_bath_part(
                        MK_big,
                        a,
                        b,
                        aggTools;
                        vib_basis=aggOperators.vib_basis
                    )
                end
            end
        end
    end
    return MemoryKernel
end
