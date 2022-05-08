
function interpolate_with_tspan(itp, tspan, t)
    t0 = tspan[1]
    t1 = tspan[end]
    N = length(tspan)
    if t < t0
        t_i = 1
    elseif t > t1
        t_i = N
    else
        t_i = (t-t0)/t1 * (N-1) + 1
    end
    return itp[t_i]
end

function M_aabb_W_bath_intp(t, s, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, _ = p
    
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    rho_t = interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    W_bath_s = interpolate_with_tspan(W_0_bath_t_itp, tspan, s)
    
    elLen = aggCore.molCount+1
    M_aabb_ = zeros(ComplexF64, elLen, elLen)
    for b=1:elLen
        rho_mantis = zeros(Float64, elLen, elLen)
        rho_mantis[b, b] = 1.
        
        tmp1[:, :] = ad(rho_mantis, W_bath_s, aggCore, aggTools)
        tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
        tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data
        M_tr = trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
        
        for a=1:elLen
            if rho_t[a, a] != 0.
                M_aabb_[a, b] = M_tr[a, a] / rho_t[a, a]
            end
        end
    end
    return M_aabb_
end

function K_aabb_W_bath_intp(t, p, tmp1, tmp2; rtol=1e-12, atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, elementtype = p
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    M_aabb_integrated, err = QuadGK.quadgk(
        s -> M_aabb_W_bath_intp(t, s, p, tmp1, tmp2, Ham_II_t),
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

function K_aabb_W_bath_intp(t, W0, W0_bath, agg::Aggregate; rtol=1e-12, atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))
    return K_aabb_W_bath_intp(t, p, tmp1, tmp2; rtol=rtol, atol=atol)
end

function M_abcd_W_bath_intp(t, s, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, _ = p
    
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    rho_t = interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    W_bath_s = interpolate_with_tspan(W_0_bath_t_itp, tspan, s)
    
    elLen = aggCore.molCount+1
    M_abcd_ = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    for c=1:elLen, d=1:elLen
        rho_mantis = zeros(Float64, elLen, elLen)
        rho_mantis[c, d] = 1.
        
        tmp1[:, :] = ad(rho_mantis, W_bath_s, aggCore, aggTools)
        tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
        tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data
        M_tr = trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
        
        for a=1:elLen, b=1:elLen
            if rho_t[a, b] != 0.
                M_abcd_[a, b, c, d] = M_tr[a, b] / rho_t[a, b]
            end
        end
    end
    return M_abcd_
end

function K_abcd_W_bath_intp(t, p, tmp1, tmp2; rtol=1e-12, atol=1e-12)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, elementtype = p
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    M_aabb_integrated, err = QuadGK.quadgk(
        s -> M_abcd_W_bath_intp(t, s, p, tmp1, tmp2, Ham_II_t),
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

function K_abcd_W_bath_intp(t, W0, W0_bath, agg::Aggregate; rtol=1e-12, atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))
    return K_abcd_W_bath_intp(t, p, tmp1, tmp2; rtol=rtol, atol=atol)
end