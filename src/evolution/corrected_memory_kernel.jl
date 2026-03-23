import QuadGK

function _wick_second_order_cf(
    c::AbstractMatrix, site_sds::Vector{SpectralDensity}, T::Real, N::Int,
    i1::Int, i2::Int, i3::Int, i4::Int, i5::Int, i6::Int, i7::Int, i8::Int,
    ta::Real, tb::Real, tc::Real, td::Real,
)::ComplexF64
    val = zero(ComplexF64)
    for n in 1:N, m in 1:N
        K_nnmm = c[i1,n]*c[i2,n]*c[i3,n]*c[i4,n] * c[i5,m]*c[i6,m]*c[i7,m]*c[i8,m]
        K_nmnm = c[i1,n]*c[i2,n]*c[i5,n]*c[i6,n] * c[i3,m]*c[i4,m]*c[i7,m]*c[i8,m]
        K_nmmn = c[i1,n]*c[i2,n]*c[i7,n]*c[i8,n] * c[i3,m]*c[i4,m]*c[i5,m]*c[i6,m]
        if abs(K_nnmm) + abs(K_nmnm) + abs(K_nmmn) < 1e-30
            continue
        end
        Cn_ab = analytic_correlation_nn(site_sds[n], ta - tb, T)
        Cm_cd = analytic_correlation_nn(site_sds[m], tc - td, T)
        Cn_ac = analytic_correlation_nn(site_sds[n], ta - tc, T)
        Cm_bd = analytic_correlation_nn(site_sds[m], tb - td, T)
        Cn_ad = analytic_correlation_nn(site_sds[n], ta - td, T)
        Cm_bc = analytic_correlation_nn(site_sds[m], tb - tc, T)
        val += K_nnmm * Cn_ab * Cm_cd + K_nmnm * Cn_ac * Cm_bd + K_nmmn * Cn_ad * Cm_bc
    end
    return val
end

function _first_order_integrand!(
    result::Array{ComplexF64, 4},
    t1::Real, t2::Real, t3::Real, t4::Real,
    E::Vector{Float64},
    rho_I::AbstractMatrix,
    elLen::Int,
    wk::Function,
)
    for a in 2:elLen, b in 2:elLen, ci in 2:elLen, d in 2:elLen
        ai, bi, cci, di = a-1, b-1, ci-1, d-1
        rho_cd = rho_I[ci, d]
        abs(rho_cd) < _SAFE_DIV_TOL && continue
        val = zero(ComplexF64)

        # --- Outer term 1: +H(t1)H(t2)O ---

        # 1a: +, efg, delta_db
        if d == b
            for e in 2:elLen, f in 2:elLen, g in 2:elLen
                ei, fi, gi = e-1, f-1, g-1
                rd = rho_I[f, d]; abs(rd) < _SAFE_DIV_TOL && continue
                ph = exp(im*((E[a]-E[e])*t1 + (E[e]-E[f])*t2 + (E[f]-E[g])*t3 + (E[g]-E[ci])*t4))
                val += ph * (rho_cd/rd) * wk(ai,ei, ei,fi, fi,gi, gi,cci, t1,t2,t3,t4)
            end
        end

        # 1b: -, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[f, b]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[e])*t1 + (E[e]-E[f])*t2 + (E[f]-E[ci])*t3 + (E[d]-E[b])*t4))
            val -= ph * (rho_cd/rd) * wk(di,bi, ai,ei, ei,fi, fi,cci, t4,t1,t2,t3)
        end

        # 1c: -, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[f, b]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[e])*t1 + (E[e]-E[f])*t2 + (E[d]-E[b])*t3 + (E[f]-E[ci])*t4))
            val -= ph * (rho_cd/rd) * wk(di,bi, ai,ei, ei,fi, fi,cci, t3,t1,t2,t4)
        end

        # 1d: +, eg (f->c)
        for e in 2:elLen, g in 2:elLen
            ei, gi = e-1, g-1
            rd = rho_I[ci, b]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[e])*t1 + (E[e]-E[ci])*t2 + (E[g]-E[b])*t3 + (E[d]-E[g])*t4))
            val += ph * (rho_cd/rd) * wk(di,gi, gi,bi, ai,ei, ei,cci, t4,t3,t1,t2)
        end

        # --- Outer term 2: -H(t1)OH(t2) ---

        # 2a: -, eg (f->d)
        for e in 2:elLen, g in 2:elLen
            ei, gi = e-1, g-1
            rd = rho_I[e, d]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[e])*t1 + (E[d]-E[b])*t2 + (E[e]-E[g])*t3 + (E[g]-E[ci])*t4))
            val -= ph * (rho_cd/rd) * wk(di,bi, ai,ei, ei,gi, gi,cci, t2,t1,t3,t4)
        end

        # 2b: +, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[e, f]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[e])*t1 + (E[f]-E[b])*t2 + (E[e]-E[ci])*t3 + (E[d]-E[f])*t4))
            val += ph * (rho_cd/rd) * wk(di,fi, fi,bi, ai,ei, ei,cci, t4,t2,t1,t3)
        end

        # 2c: +, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[e, f]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[e])*t1 + (E[f]-E[b])*t2 + (E[d]-E[f])*t3 + (E[e]-E[ci])*t4))
            val += ph * (rho_cd/rd) * wk(di,fi, fi,bi, ai,ei, ei,cci, t3,t2,t1,t4)
        end

        # 2d: -, fg (e->c)
        for f in 2:elLen, g in 2:elLen
            fi, gi = f-1, g-1
            rd = rho_I[ci, f]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[a]-E[ci])*t1 + (E[f]-E[b])*t2 + (E[g]-E[f])*t3 + (E[d]-E[g])*t4))
            val -= ph * (rho_cd/rd) * wk(di,gi, gi,fi, fi,bi, ai,cci, t4,t3,t2,t1)
        end

        # --- Outer term 3: -H(t2)OH(t1) ---

        # 3a: -, eg (f->d)
        for e in 2:elLen, g in 2:elLen
            ei, gi = e-1, g-1
            rd = rho_I[e, d]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[d]-E[b])*t1 + (E[a]-E[e])*t2 + (E[e]-E[g])*t3 + (E[g]-E[ci])*t4))
            val -= ph * (rho_cd/rd) * wk(di,bi, ai,ei, ei,gi, gi,cci, t1,t2,t3,t4)
        end

        # 3b: +, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[e, f]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[f]-E[b])*t1 + (E[a]-E[e])*t2 + (E[e]-E[ci])*t3 + (E[d]-E[f])*t4))
            val += ph * (rho_cd/rd) * wk(di,fi, fi,bi, ai,ei, ei,cci, t4,t1,t2,t3)
        end

        # 3c: +, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[e, f]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[f]-E[b])*t1 + (E[a]-E[e])*t2 + (E[d]-E[f])*t3 + (E[e]-E[ci])*t4))
            val += ph * (rho_cd/rd) * wk(di,fi, fi,bi, ai,ei, ei,cci, t3,t1,t2,t4)
        end

        # 3d: -, fg (e->c)
        for f in 2:elLen, g in 2:elLen
            fi, gi = f-1, g-1
            rd = rho_I[ci, f]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[f]-E[b])*t1 + (E[a]-E[ci])*t2 + (E[g]-E[f])*t3 + (E[d]-E[g])*t4))
            val -= ph * (rho_cd/rd) * wk(di,gi, gi,fi, fi,bi, ai,cci, t4,t3,t1,t2)
        end

        # --- Outer term 4: +OH(t2)H(t1) ---

        # 4a: +, fg (e->d)
        for f in 2:elLen, g in 2:elLen
            fi, gi = f-1, g-1
            rd = rho_I[a, d]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[f]-E[b])*t1 + (E[d]-E[f])*t2 + (E[a]-E[g])*t3 + (E[g]-E[ci])*t4))
            val += ph * (rho_cd/rd) * wk(di,fi, fi,bi, ai,gi, gi,cci, t2,t1,t3,t4)
        end

        # 4b: -, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[a, e]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[f]-E[b])*t1 + (E[e]-E[f])*t2 + (E[a]-E[ci])*t3 + (E[d]-E[e])*t4))
            val -= ph * (rho_cd/rd) * wk(di,ei, ei,fi, fi,bi, ai,cci, t4,t2,t1,t3)
        end

        # 4c: -, ef
        for e in 2:elLen, f in 2:elLen
            ei, fi = e-1, f-1
            rd = rho_I[a, e]; abs(rd) < _SAFE_DIV_TOL && continue
            ph = exp(im*((E[f]-E[b])*t1 + (E[e]-E[f])*t2 + (E[d]-E[e])*t3 + (E[a]-E[ci])*t4))
            val -= ph * (rho_cd/rd) * wk(di,ei, ei,fi, fi,bi, ai,cci, t3,t2,t1,t4)
        end

        # 4d: +, efg, delta_ac
        if a == ci
            for e in 2:elLen, f in 2:elLen, g in 2:elLen
                ei, fi, gi = e-1, f-1, g-1
                rd = rho_I[ci, e]; abs(rd) < _SAFE_DIV_TOL && continue
                ph = exp(im*((E[f]-E[b])*t1 + (E[e]-E[f])*t2 + (E[g]-E[e])*t3 + (E[d]-E[g])*t4))
                val += ph * (rho_cd/rd) * wk(di,gi, gi,ei, ei,fi, fi,bi, t4,t3,t2,t1)
            end
        end

        result[a, b, ci, d] = val
    end
end

function first_order_memory_kernel_cf(
    t1::Real, t2::Real,
    aggCore::AggregateCore,
    rho_I::AbstractMatrix;
    T::Real = 300.0,
    rtol::Real = 1e-6,
    atol::Real = 1e-8,
)
    M0 = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = T)
    if t2 <= 0.0
        return M0
    end

    aggTools = AggregateTools(aggCore)
    energies, coefficients = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount
    elLen = N + 1

    site_sds = [spectral_density(aggCore.molecules[k]) for k in 1:N]

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools;
        vib_basis = GroundGround(), groundEnergy = true)
    E = zeros(Float64, elLen)
    E[1] = real(Ham_sys.data[1, 1])
    E[2:end] .= energies

    cc = coefficients

    wk(i1,i2,i3,i4,i5,i6,i7,i8, ta,tb,tc,td) =
        _wick_second_order_cf(cc, site_sds, T, N, i1,i2,i3,i4,i5,i6,i7,i8, ta,tb,tc,td)

    correction, _ = QuadGK.quadgk(0.0, t2; rtol=rtol, atol=atol) do t3
        inner, _ = QuadGK.quadgk(0.0, t3; rtol=rtol, atol=atol) do t4
            result = zeros(ComplexF64, elLen, elLen, elLen, elLen)
            _first_order_integrand!(result, t1, t2, t3, t4, E, rho_I, elLen, wk)
            return result
        end
        return inner
    end

    return M0 .- correction
end

function equilibrium_bath_state(
    aggCore::AggregateCore,
    aggTools::AggregateTools,
    aggOperators::AggregateOperators,
    T::Real,
)
    indicesMap = aggTools.indicesMap
    elLen = aggCore.molCount + 1
    bSize = aggTools.bSize

    a1, a2 = indicesMap[1][1], indicesMap[1][end]
    H_0_g = real.(aggOperators.Ham_0.data[a1:a2, a1:a2])
    rho_B = exp(-H_0_g / (BOLTZMANN_CM * T))
    rho_B ./= LinearAlgebra.tr(rho_B)

    W_eq = zeros(ComplexF64, bSize, bSize)
    for a in 1:elLen, b in 1:elLen
        ai1, ai2 = indicesMap[a][1], indicesMap[a][end]
        bi1, bi2 = indicesMap[b][1], indicesMap[b][end]
        W_eq[ai1:ai2, bi1:bi2] .= rho_B
    end
    return W_eq
end

function zeroth_order_memory_kernel_cf(
    t1::Real, t2::Real,
    aggCore::AggregateCore;
    T::Real = 300.0,
)
    aggTools = AggregateTools(aggCore)
    energies, coefficients = exciton_basis(aggCore, aggTools)
    N = aggCore.molCount
    elLen = N + 1

    site_sds = [spectral_density(aggCore.molecules[k]) for k in 1:N]

    tau = t1 - t2
    C_fwd = ComplexF64[analytic_correlation_nn(sd, tau, T) for sd in site_sds]
    C_bwd = ComplexF64[analytic_correlation_nn(sd, -tau, T) for sd in site_sds]

    Ham_sys = get_agg_ham_system_small(aggCore, aggTools;
        vib_basis = GroundGround(), groundEnergy = true)
    E = zeros(Float64, elLen)
    E[1] = real(Ham_sys.data[1, 1])
    E[2:end] .= energies

    c = coefficients
    M = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    for a in 2:elLen, b in 2:elLen, ci in 2:elLen, d in 2:elLen
        ai, bi, cci, di = a - 1, b - 1, ci - 1, d - 1
        val = zero(ComplexF64)

        if d == b
            for e in 2:elLen
                ei = e - 1
                phase = exp(im * ((E[a] - E[e]) * t1 + (E[e] - E[ci]) * t2))
                val += phase * exciton_correlation(C_fwd, c, ai, ei, ei, cci)
            end
        end

        if a == ci
            for e in 2:elLen
                ei = e - 1
                phase = exp(im * ((E[d] - E[e]) * t2 + (E[e] - E[b]) * t1))
                val += phase * exciton_correlation(C_bwd, c, di, ei, ei, bi)
            end
        end

        phase3 = exp(im * ((E[a] - E[ci]) * t1 + (E[d] - E[b]) * t2))
        val -= phase3 * exciton_correlation(C_bwd, c, di, bi, ai, cci)

        phase4 = exp(im * ((E[a] - E[ci]) * t2 + (E[d] - E[b]) * t1))
        val -= phase4 * exciton_correlation(C_fwd, c, di, bi, ai, cci)

        M[a, b, ci, d] = val
    end

    return M
end

function site_to_exciton_kernel(
    M_site::Array{ComplexF64, 4},
    coefficients::AbstractMatrix,
)
    elLen = size(M_site, 1)
    N = elLen - 1

    T_mat = zeros(Float64, elLen, elLen)
    T_mat[1, 1] = 1.0
    T_mat[2:end, 2:end] .= coefficients

    M_exc = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    for α in 1:elLen, β in 1:elLen, γ in 1:elLen, δ in 1:elLen
        s = zero(ComplexF64)
        for a in 1:elLen, b in 1:elLen, cc in 1:elLen, d in 1:elLen
            s += T_mat[α, a] * T_mat[β, b] * T_mat[γ, cc] * T_mat[δ, d] * M_site[a, b, cc, d]
        end
        M_exc[α, β, γ, δ] = s
    end
    return M_exc
end
