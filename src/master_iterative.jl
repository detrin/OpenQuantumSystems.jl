using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq, Interpolations


function normalize_bath(W_bath, aggCore, aggTools, aggOperators)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    W_bath_tr = trace_bath(W_bath, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
    for b=1:elLen
        for a=1:elLen
            a1 = indicesMap[a][1]
            a2 = indicesMap[a][end]
            b1 = indicesMap[b][1]
            b2 = indicesMap[b][end]
            if abs(W_bath_tr[a, b]) >= _SAFE_DIV_TOL
                W_bath[a1:a2, b1:b2] /= W_bath_tr[a, b]
            end
        end
    end
    return W_bath
end


function W_abcd_1_bath_core(t, t1, t2, p, tmp1, tmp2)
    (; aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp) = p

    # W_1_bath = deepcopy(W0_bath.data)
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I

    Ham_II_t1 = get_interaction_ham_i_picture(Ham_0, Ham_I, t1)
    Ham_II_t2 = get_interaction_ham_i_picture(Ham_0, Ham_I, t2)

    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    rho_t2 = interpolate_with_tspan(rho_0_int_t_itp, tspan, t2)
    W_bath_t2 = interpolate_with_tspan(W_0_bath_t_itp, tspan, t2)
    rho_t1 = interpolate_with_tspan(rho_0_int_t_itp, tspan, t1)
    rho_t1[:, :] = map(_safe_inv, rho_t1)

    tmp1[:, :] = ad(rho_t2, W_bath_t2, aggCore, aggTools)
    tmp1[:, :] = W_bath_t2[:, :]
    tmp2[:, :] = Ham_II_t2.data * tmp1 - tmp1 * Ham_II_t2.data
    tmp1[:, :] = Ham_II_t1.data * tmp2 - tmp2 * Ham_II_t1.data
    W_1_bath[:, :] = ad(rho_t1, tmp1, aggCore, aggTools)
    #=
    for a=2:elLen, b=1:elLen
        a1 = indicesMap[a][1]
        a2 = indicesMap[a][end]
        b1 = indicesMap[b][1]
        b2 = indicesMap[b][end]
        if rho_t1[a, b] != 0
            W_1_bath[a1:a2, b1:b2] = W_1_bath[a1:a2, b1:b2]  / rho_t1[a, b]
        end
    end
    =#

    return W_1_bath
end

function W_1_bath(
        t, p, tmp1, tmp2;
        normalize=false,
        W_1_rtol=1e-12, W_1_atol=1e-12
    )
    (; aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp) = p

    W_1_diff, err = QuadGK.quadgk(
        t1 -> begin
            W_1_diff_t1, err = QuadGK.quadgk(
                t2 -> W_abcd_1_bath_core(t, t1, t2, p, tmp1, tmp2),
                0,
                t1,
                rtol = W_1_rtol,
                atol = W_1_atol,
            );
            return W_1_diff_t1
        end,
        0,
        t,
        rtol = W_1_rtol,
        atol = W_1_atol,
    )
    rho_t = interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    W_bath_t = interpolate_with_tspan(W_0_bath_t_itp, tspan, t)
    W_1_bath = deepcopy(W_bath_t) - W_1_diff
    if normalize
        W_1_bath = normalize_bath(W_1_bath, aggCore, aggTools, aggOperators)
    end
    return W_1_bath
end

function W_1_bath(t, W0, W0_bath, agg::Aggregate; W_1_rtol=1e-12, W_1_atol=1e-12, K_rtol=1e-12, K_atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (aggCore=agg.core, aggTools=agg.tools, aggOperators=agg.operators, W0=W0, W0_bath=W0_bath, elementtype=eltype(W0))
    return W_1_bath(t, p, tmp1, tmp2; W_1_rtol=W_1_rtol, W_1_atol=W_1_atol, K_rtol=K_rtol, K_atol=K_atol)
end

function QME_sI_iterative(
    W0::T,
    rho_0_int_t,
    W_0_bath_t,
    tspan::AbstractVector,
    agg::Aggregate;
    method::Symbol = :default,
    normalize=false,
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    W_1_rtol::AbstractFloat = 1e-12,
    W_1_atol::AbstractFloat = 1e-12,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    method ∈ (:default, :markov0, :markov1) ||
        throw(ArgumentError("method must be :default, :markov0, or :markov1, got :$method"))
    setup = _setup_delayed_integration(W0, tspan, agg)
    (; rho0, history_fun, tmp1, tmp2, tspan_, x0, state, dstate) = setup
    W0_bath = get_rho_bath(W0, agg.core, agg.operators, agg.tools; vib_basis=agg.operators.vib_basis)

    # Calculate and interpolate rho_0_int_t, W_0_bath_t, W_1_bath_t
    if ndims(rho_0_int_t) == 1 && rho_0_int_t[1] isa Operator
        rho_0_int_t = operator_recast(rho_0_int_t)
    end
    if ndims(rho_0_int_t) == 3
        rho_0_int_t = [rho_0_int_t[i, :, :] for i=1:size(rho_0_int_t, 1)]
    end
    rho_0_int_t_itp = Interpolations.interpolate(
        rho_0_int_t,
        Interpolations.BSpline(Interpolations.Linear())
    )

    if ndims(W_0_bath_t) == 1 && W_0_bath_t[1] isa Operator
        W_0_bath_t = operator_recast(W_0_bath_t)
    end
    if ndims(W_0_bath_t) == 3
        W_0_bath_t = [W_0_bath_t[i, :, :] for i=1:size(W_0_bath_t, 1)]
    end
    W_0_bath_t_itp = Interpolations.interpolate(
        W_0_bath_t,
        Interpolations.BSpline(Interpolations.Linear())
    )
    W_1_bath_fn = method == :markov0 ? W_1_markov0_bath :
                  method == :markov1 ? W_1_markov1_bath :
                  W_1_bath
    p = (aggCore=agg.core, aggTools=agg.tools, aggOperators=agg.operators, W0=W0, W0_bath=W0_bath, tspan=tspan, rho_0_int_t_itp=rho_0_int_t_itp, W_0_bath_t_itp=W_0_bath_t_itp, elementtype=eltype(W0))
    elLen = agg.core.molCount+1
    W_1_bath_t = []
    for t_i=1:length(tspan)
        t = tspan[t_i]
        W_1_bath_ = W_1_bath_fn(
            t, p, tmp1, tmp2;
            normalize=normalize,
            W_1_rtol=W_1_rtol, W_1_atol=W_1_atol
        )
        push!(W_1_bath_t, W_1_bath_)
    end
    W_1_bath_itp = Interpolations.interpolate(W_1_bath_t, Interpolations.BSpline(Interpolations.Linear()))

    p = (aggCore=agg.core, aggTools=agg.tools, aggOperators=agg.operators, W0=W0, W0_bath=W0_bath, rho_0_int_t_itp=rho_0_int_t_itp, W_1_bath_itp=W_1_bath_itp, tspan=tspan, elementtype=eltype(W0))
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_iterative(
        t,
        rho,
        drho,
        history_fun,
        tmp1,
        tmp2,
        p,
        int_reltol,
        int_abstol,
    )
    tspan, rho_int_1_t = integrate_delayed(
        tspan_,
        dmaster_,
        history_fun,
        x0,
        state,
        dstate,
        fout;
        p = p,
        reltol = reltol,
        abstol = abstol,
        alg = alg,
        kwargs...,
    )
    return tspan, rho_int_1_t, W_1_bath_t
end

function QME_sI_iterative_markov0(W0, rho_0_int_t, W_0_bath_t, tspan, agg; kwargs...)
    return QME_sI_iterative(W0, rho_0_int_t, W_0_bath_t, tspan, agg; method=:markov0, kwargs...)
end

function QME_sI_iterative_markov1(W0, rho_0_int_t, W_0_bath_t, tspan, agg; kwargs...)
    return QME_sI_iterative(W0, rho_0_int_t, W_0_bath_t, tspan, agg; method=:markov1, kwargs...)
end

function dQME_sI_iterative(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp1::AbstractMatrix,
    tmp2::AbstractMatrix,
    p,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
) where {B<:Basis,T<:Operator{B,B}}
    (; aggCore, aggTools, aggOperators, W0, elementtype) = p

    Ham_II_t = get_interaction_ham_i_picture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
    # println("dQME_sI_iterative", " ", t)
    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_iterative(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_iterative(t, s, h, p, tmp1, tmp2, Ham_II_t)
    (; aggCore, aggTools, aggOperators, W0, W0_bath, rho_0_int_t_itp, W_1_bath_itp, tspan) = p

    rho = h(p, s)

    rho = _to_matrix(rho)

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = get_interaction_ham_i_picture(Ham_0, Ham_I, s)
    W_bath_s = interpolate_with_tspan(W_1_bath_itp, tspan, s)

    tmp1[:, :] = ad(rho, W_bath_s, aggCore, aggTools)
    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggOperators, aggTools; vib_basis=aggOperators.vib_basis)
end

function W_abcd_1_markov0_bath_core(t, t1, t2, p, tmp1, tmp2)
    (; aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp) = p

    # W_1_bath = deepcopy(W0_bath.data)
    # elLen = aggCore.molCount+1
    # indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I

    Ham_II_t1 = get_interaction_ham_i_picture(Ham_0, Ham_I, t1)
    Ham_II_t2 = get_interaction_ham_i_picture(Ham_0, Ham_I, t2)

    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    rho_t2 = interpolate_with_tspan(rho_0_int_t_itp, tspan, t2)
    W_bath_t2 = interpolate_with_tspan(W_0_bath_t_itp, tspan, t2)

    tmp1[:, :] = ad(rho_t2, W_bath_t2, aggCore, aggTools)
    # tmp1[:, :] = W_bath_t2[:, :]
    tmp2[:, :] = Ham_II_t2.data * tmp1 - tmp1 * Ham_II_t2.data
    W_1_bath[:, :] = Ham_II_t1.data * tmp2 - tmp2 * Ham_II_t1.data

    return W_1_bath
end

function W_1_markov0_bath(
        t, p, tmp1, tmp2;
        normalize=false,
        W_1_rtol=1e-12, W_1_atol=1e-12
    )
    (; aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp) = p
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap

    W_1_diff, err = QuadGK.quadgk(
        t1 -> begin
            W_1_diff_t1, err = QuadGK.quadgk(
                t2 -> W_abcd_1_markov0_bath_core(t, t1, t2, p, tmp1, tmp2),
                0,
                t1,
                rtol = W_1_rtol,
                atol = W_1_atol,
            );
            return W_1_diff_t1
        end,
        0,
        t,
        rtol = W_1_rtol,
        atol = W_1_atol,
    )
    rho_t = interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    W_bath_t = interpolate_with_tspan(W_0_bath_t_itp, tspan, t)
    for a=2:elLen, b=1:elLen
        a1 = indicesMap[a][1]
        a2 = indicesMap[a][end]
        b1 = indicesMap[b][1]
        b2 = indicesMap[b][end]
        if abs(rho_t[a, b]) < _SAFE_DIV_TOL
            tmp1[a1:a2, b1:b2] = W_1_diff[a1:a2, b1:b2]
        else
            tmp1[a1:a2, b1:b2] = W_1_diff[a1:a2, b1:b2] / rho_t[a, b]
        end
    end
    W_1_bath = deepcopy(W_bath_t) - tmp1
    if normalize
        W_1_bath = normalize_bath(W_1_bath, aggCore, aggTools, aggOperators)
    end
    return W_1_bath
end

function W_abcd_1_markov1_bath_core(t, t1, t2, p, tmp1, tmp2)
    (; aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp) = p

    # W_1_bath = deepcopy(W0_bath.data)
    # elLen = aggCore.molCount+1
    # indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I

    Ham_II_t1 = get_interaction_ham_i_picture(Ham_0, Ham_I, t1)
    Ham_II_t2 = get_interaction_ham_i_picture(Ham_0, Ham_I, t2)

    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    rho_t2 = interpolate_with_tspan(rho_0_int_t_itp, tspan, t2)
    # W_bath_t1 instead of W_bath_t2
    W_bath_t1 = interpolate_with_tspan(W_0_bath_t_itp, tspan, t1)

    tmp1[:, :] = ad(rho_t2, W_bath_t1, aggCore, aggTools)
    # tmp1[:, :] = W_bath_t2[:, :]
    tmp2[:, :] = Ham_II_t2.data * tmp1 - tmp1 * Ham_II_t2.data
    W_1_bath[:, :] = Ham_II_t1.data * tmp2 - tmp2 * Ham_II_t1.data

    return W_1_bath
end

function W_1_markov1_bath(
        t, p, tmp1, tmp2;
        normalize=false,
        W_1_rtol=1e-12, W_1_atol=1e-12
    )
    (; aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp) = p
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap

    W_1_diff, err = QuadGK.quadgk(
        t1 -> begin
            W_1_diff_t1, err = QuadGK.quadgk(
                t2 -> W_abcd_1_markov1_bath_core(t, t1, t2, p, tmp1, tmp2),
                0,
                t1,
                rtol = W_1_rtol,
                atol = W_1_atol,
            );
            return W_1_diff_t1
        end,
        0,
        t,
        rtol = W_1_rtol,
        atol = W_1_atol,
    )
    rho_t = interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    W_bath_t = interpolate_with_tspan(W_0_bath_t_itp, tspan, t)
    for a=2:elLen, b=1:elLen
        a1 = indicesMap[a][1]
        a2 = indicesMap[a][end]
        b1 = indicesMap[b][1]
        b2 = indicesMap[b][end]
        if abs(rho_t[a, b]) < _SAFE_DIV_TOL
            tmp1[a1:a2, b1:b2] = W_1_diff[a1:a2, b1:b2]
        else
            tmp1[a1:a2, b1:b2] = W_1_diff[a1:a2, b1:b2] / rho_t[a, b]
        end
    end
    W_1_bath = deepcopy(W_bath_t) - tmp1
    if normalize
        W_1_bath = normalize_bath(W_1_bath, aggCore, aggTools, aggOperators)
    end
    return W_1_bath
end
