
function M_aabb(t, s, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    U_0_op = evolutionOperator(Ham_0, s)
    W0_int_s = U_0_op' * W0_bath * U_0_op
    
    elLen = aggCore.molCount+1
    M_aabb_ = zeros(ComplexF64, elLen, elLen)
    for b=1:elLen
        rho_mantis = zeros(Float64, elLen, elLen)
        rho_mantis[b, b] = 1.
        
        tmp1[:, :] = ad(rho_mantis, W0_int_s.data, aggCore, aggTools)
        tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
        tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data
        M_tr = trace_bath(tmp1, aggCore, aggTools)
        for a=1:elLen
            M_aabb_[a, b] = M_tr[a, a]
        end
    end
    return M_aabb_
end

function K_const_ab(t, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_t = getInteractionHamIPicture(Ham_0, Ham_I, t)

    elLen = aggCore.molCount+1
    K_const_ab_ = zeros(ComplexF64, elLen, elLen)
    for b=1:elLen
        rho_mantis = zeros(Float64, elLen, elLen)
        rho_mantis[b, b] = 1.
        
        tmp1[:, :] = ad(rho_mantis, W0_bath.data, aggCore, aggTools)
        tmp2[:, :] = Ham_II_t.data * tmp1 - tmp1 * Ham_II_t.data
        K_const_tr = trace_bath(tmp1, aggCore, aggTools)
        for a=1:elLen
            K_const_ab_[a, b] = K_const_tr[a, a]
        end
    end
    return K_const_ab_
end

function K_ab(t, p, tmp1, tmp2)
    aggCore, aggTools, aggOperators, W0, W0_bath, elementtype = p
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    M_aabb_integrated, err = QuadGK.quadgk(
        s -> M_aabb(t, s, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = 1e-12,
        atol = 1e-12,
    )
    # K_const_ab_ = K_const_ab(t, p, tmp1, tmp2, Ham_II_t)
    # K_ab_ = -elementtype(im) * K_const_ab_ - M_aabb_integrated
    K_ab_ = - M_aabb_integrated
    return K_ab_
end

function W_1_bath_core(t, s, p, tmp1, tmp2)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    # W_1_bath = deepcopy(W0_bath.data)
    indicesMap = aggTools.indicesMap
    Ham = aggOperators.Ham
    
    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    K_ab_s = K_ab(s, p, tmp1, tmp2)
    for b=1:elLen
        U_bb_ = evolution_el_part(Ham.data, s, b, b, indicesMap)
        U_0_bb = evolution_el_part(Ham_0.data, s, b, b, indicesMap)
        U_bb = U_0_bb * U_bb_
        for a=1:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            U_aa_ = evolution_el_part(Ham.data, t-s, a, a, indicesMap)
            U_0_aa = evolution_el_part(Ham_0.data, t-s, a, a, indicesMap)
            U_aa = U_0_aa * U_aa_
            U_ = U_aa * U_bb
            W_1_bath[a1:a2, a1:a2] = W_1_bath[a1:a2, a1:a2] + 
                K_ab_s[a, b] * U_ * W0_bath.data[b1:b2, b1:b2] * adjoint(U_)
        end
    end
    return W_1_bath
end

function W_1_bath(t, p, tmp1, tmp2)
    aggCore, aggTools, aggOperators, W0, W0_bath, _ = p
    
    W_1_diff, err = QuadGK.quadgk(
        s -> W_1_bath_core(t, s, p, tmp1, tmp2),
        0,
        t,
        rtol = 1e-12,
        atol = 1e-12,
    )
    
    W_1_bath = deepcopy(W0_bath.data) + W_1_diff
    return W_1_bath
end

function normalize_bath(W_bath, aggCore, aggTools)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    W_bath_tr = trace_bath(W_bath, agg.core, agg.tools)
    for b=1:elLen
        for a=1:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            W_bath[a1:a2, b1:b2] /= W_bath[a, b]
        end
    end
    return W_bath
end