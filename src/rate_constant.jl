
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
        M_tr = trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
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
        K_const_tr = trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
        for a=1:elLen
            K_const_ab_[a, b] = K_const_tr[a, a]
        end
    end
    return K_const_ab_
end

function K_ab(t, p, tmp1, tmp2; rtol=1e-12, atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, elementtype = p
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    M_aabb_integrated, err = QuadGK.quadgk(
        s -> M_aabb(t, s, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = rtol,
        atol = atol,
    )
    # K_const_ab_ = K_const_ab(t, p, tmp1, tmp2, Ham_II_t)
    # K_ab_ = -elementtype(im) * K_const_ab_ - M_aabb_integrated
    K_ab_ = - M_aabb_integrated
    return K_ab_
end

function K_ab(t, W0, W0_bath, agg::Aggregate; rtol=1e-12, atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))
    return K_ab(t, p, tmp1, tmp2; rtol=rtol, atol=atol)
end
