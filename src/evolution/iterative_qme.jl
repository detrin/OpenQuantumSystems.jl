# Master equations for an aggregate coupled to an arbitrary harmonic bath (finite or
# infinite) expressed through bath correlation functions in exponential-sum form, see
# [`ExponentialCorrelation`](@ref). Everything is formulated in the exciton basis of the
# single-excitation manifold and in the interaction picture with respect to H_S + H_B.
#
# Every memory-kernel contribution has the form
#     dρ_ab/dt ∋ coeff · e^{iΩt} ∫_0^t dτ e^{-Ẑ(t-τ)} B(τ) ρ_cd(τ),
# with Re Ẑ ≥ 0, so the history integrals are propagated recursively instead of being
# re-integrated at every time step.

import LinearAlgebra

# ---------------------------------------------------------------------------
# Exponential trapezoid rule
# ---------------------------------------------------------------------------

# ∫_0^h e^{-Z(h-s)} [(1 - s/h) f0 + (s/h) f1] ds = w0 f0 + w1 f1, and decay = e^{-Zh};
# exact for f linear on the step.
function _exp_trapezoid(Z::Number, h::Real)
    u = Z * h
    if abs(u) < 1e-3
        w0 = h * (1 / 2 - u / 3 + u^2 / 8 - u^3 / 30)
        w1 = h * (1 / 2 - u / 6 + u^2 / 24 - u^3 / 120)
    else
        e = exp(-u)
        em = -expm1(-u)
        w1 = (em * (u - 1) + u * e) / (Z * u)
        w0 = em / Z - w1
    end
    return ComplexF64(exp(-u)), ComplexF64(w0), ComplexF64(w1)
end

# y(t_n) = ∫_0^{t_n} e^{-Z(t_n - s)} f(s) ds on a uniform grid
function _exp_convolution(Z::Number, f::AbstractVector, h::Real)
    d, w0, w1 = _exp_trapezoid(Z, h)
    y = zeros(ComplexF64, length(f))
    for n in 2:length(f)
        y[n] = d * y[n-1] + w0 * f[n-1] + w1 * f[n]
    end
    return y
end

# ∫_0^t e^{-Z(t - s)} ds
_exp_integral(Z::Number, t::Real) = abs(Z * t) < 1e-12 ? complex(float(t)) : -expm1(-Z * t) / Z

# ---------------------------------------------------------------------------
# Operator chains and Wick pairings
# ---------------------------------------------------------------------------

# H_I(t_slot)_{αβ} = e^{i(E_α - E_β) t_slot} Σ_n c_αn c_βn ΔV_n(t_slot)
struct _ChainOp
    slot::Int
    α::Int
    β::Int
end

_manifold_paths(α, β, nops, N) =
    nops == 0 ? (α == β ? [(α,)] : NTuple{1,Int}[]) :
    nops == 1 ? [(α, β)] : [(α, g, β) for g in 1:N]

function _append_ops!(chain, slots, path)
    for (j, slot) in enumerate(slots)
        push!(chain, _ChainOp(slot, path[j], path[j+1]))
    end
    return chain
end

# <ΔV(t_p) ΔV(t_q)> for chain positions p < q, site n and exponential component k.
# Returns the coefficient, the exponent and the slot span (lo < hi) it decays over;
# slots are ordered t_1 ≥ t_2 ≥ t_3 ≥ t_4.
function _pair_factor(model::ExcitonModel, op_p::_ChainOp, op_q::_ChainOp, n::Int, k::Int)
    c = model.coefficients
    C = model.baths[n]
    K = c[op_p.α, n] * c[op_p.β, n] * c[op_q.α, n] * c[op_q.β, n]
    if op_p.slot < op_q.slot
        return K * C.a[k], C.z[k], op_p.slot, op_q.slot
    else
        return K * conj(C.a[k]), conj(C.z[k]), op_q.slot, op_p.slot
    end
end

_round_key(x::Real) = round(x; digits = 9)
_round_key(x::Complex) = (round(real(x); digits = 9), round(imag(x); digits = 9))

# Contribution coeff · e^{iΩt} ∫_0^t e^{-Ẑ(t-τ)} B(τ) ρ_cd(τ) dτ to dρ_ab/dt;
# B === nothing means B ≡ 1.
struct _ConvTerm
    out::Int
    inp::Int
    Ω::Float64
    Ẑ::ComplexF64
    coeff::ComplexF64
    B::Union{Nothing,Vector{ComplexF64}}
end

_lin(a, b, N) = a + (b - 1) * N

"""
Second-order memory kernel ``M_{abcd}(t, \\tau; \\hat w^{0})`` (thesis, eq. MF_w_0_2) as
convolution terms; it enters as ``\\dot\\rho_{ab} = -\\sum_{cd}\\int_0^t M_{abcd}\\rho_{cd}``.
"""
function _second_order_terms(model::ExcitonModel)
    N = length(model.energies)
    E = model.energies
    acc = Dict{Tuple{Int,Int,Float64,Tuple{Float64,Float64}},Tuple{Float64,ComplexF64,ComplexF64}}()
    for a in 1:N, b in 1:N, c in 1:N, d in 1:N
        for (s, L, R) in _double_commutator_terms(1, 2)
            for pL in _manifold_paths(a, c, length(L), N), pR in _manifold_paths(d, b, length(R), N)
                # tr_B{L O R} = <R L>
                chain = _ChainOp[]
                _append_ops!(chain, R, pR)
                _append_ops!(chain, L, pL)
                Ω = zeros(2)
                for op in chain
                    Ω[op.slot] += E[op.α] - E[op.β]
                end
                for n in 1:N, k in eachindex(model.baths[n].a)
                    coef, z, _, _ = _pair_factor(model, chain[1], chain[2], n, k)
                    coef == 0 && continue
                    ΩA = Ω[1] + Ω[2]
                    Ẑ = z + im * Ω[2]
                    key = (_lin(a, b, N), _lin(c, d, N), _round_key(ΩA), _round_key(Ẑ))
                    _, _, prev = get(acc, key, (ΩA, Ẑ, zero(ComplexF64)))
                    acc[key] = (ΩA, Ẑ, prev + s * coef)
                end
            end
        end
    end
    return [_ConvTerm(k[1], k[2], v[1], v[2], -v[3], nothing) for (k, v) in acc if v[3] != 0]
end

"""
First-order correction of the memory kernel for a known reduced density matrix ρ¹(t)
on the grid (time × exciton × exciton). The correction integrand is

``X_{abcd,ef}(t_1,t_2,t_3,t_4) = \\operatorname{tr}_B\\langle a|[H(t_1),[H(t_2),
\\langle c|[H(t_3),[H(t_4), |e\\rangle\\langle f| w_{eq}]]|d\\rangle\\, |c\\rangle\\langle d|]]|b\\rangle``,

integrated over ``t_1 ≥ t_2 ≥ t_3 ≥ t_4 ≥ 0``. Variants:
- `:thesis` — thesis eq. for ``M(t_1, t_2; \\hat w^{1})``: all Wick pairings, weight
  ``\\rho^1_{ef}(t_4)/\\rho^1_{cd}(t_3)``; enters as a kernel acting on ``\\rho_{cd}(t_2)``.
- `:projected` — as `:thesis`, plus the term ``-\\hat w\\,\\dot\\rho/\\rho`` that the
  quotient rule adds to the RBP equation; it removes the disconnected pairing.
- `:nz` — `:projected` with ``\\rho_{cd}(t_2)/\\rho^1_{cd}(t_3) \\to 1``: a source term
  ``\\sum \\int X\\, \\rho^1_{ef}(t_4)`` (fourth-order Nakajima-Zwanzig term).
"""
function _fourth_order_terms(model::ExcitonModel, ρ1::Array{ComplexF64,3}, h::Real, variant::Symbol)
    variant ∈ (:thesis, :projected, :nz) ||
        throw(ArgumentError("correction must be :thesis, :projected or :nz, got :$variant"))
    N = length(model.energies)
    E = model.energies
    pairings = (((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3)))
    gap_exponent(lo, hi, z, g) = (lo <= g && hi >= g + 1) ? z : zero(z)

    # plan[Kkey][Fkey] = (Ω, Ẑ12', coeff); Hkey/Kkey carry the inner exponents
    HZ = Dict{Any,ComplexF64}()
    KZ = Dict{Any,Tuple{Any,ComplexF64,Int,Int}}()
    plan = Dict{Any,Dict{Any,Tuple{Float64,ComplexF64,ComplexF64}}}()
    outer = _double_commutator_terms(1, 2)
    inner = _double_commutator_terms(3, 4)
    for a in 1:N, b in 1:N, c in 1:N, d in 1:N, e in 1:N, f in 1:N
        for (s_o, L_o, R_o) in outer, (s_i, L_i, R_i) in inner
            for p_Lo in _manifold_paths(a, c, length(L_o), N),
                p_Li in _manifold_paths(c, e, length(L_i), N),
                p_Ri in _manifold_paths(f, d, length(R_i), N),
                p_Ro in _manifold_paths(d, b, length(R_o), N)

                # tr_B{L_o L_i w R_i R_o} = <R_i R_o L_o L_i>
                chain = _ChainOp[]
                _append_ops!(chain, R_i, p_Ri)
                _append_ops!(chain, R_o, p_Ro)
                _append_ops!(chain, L_o, p_Lo)
                _append_ops!(chain, L_i, p_Li)
                Ω = zeros(4)
                for op in chain
                    Ω[op.slot] += E[op.α] - E[op.β]
                end
                for ((p, q), (r, s)) in pairings
                    disconnected = Set((chain[p].slot, chain[q].slot)) in (Set((1, 2)), Set((3, 4)))
                    disconnected && variant != :thesis && continue
                    for n in 1:N, m in 1:N,
                        k in eachindex(model.baths[n].a), l in eachindex(model.baths[m].a)

                        c1, z1, lo1, hi1 = _pair_factor(model, chain[p], chain[q], n, k)
                        c1 == 0 && continue
                        c2, z2, lo2, hi2 = _pair_factor(model, chain[r], chain[s], m, l)
                        c2 == 0 && continue
                        Zg = [gap_exponent(lo1, hi1, z1, g) + gap_exponent(lo2, hi2, z2, g) for g in 1:3]
                        Ẑ34 = Zg[3] + im * Ω[4]
                        Ẑ23 = Zg[2] + im * (Ω[3] + Ω[4])
                        Ẑ12 = Zg[1] + im * (Ω[2] + Ω[3] + Ω[4])
                        Ωtot = sum(Ω)
                        hkey = (e, f, _round_key(Ẑ34))
                        kkey = (c, d, _round_key(Ẑ23), hkey)
                        fkey = (_lin(a, b, N), _lin(c, d, N), _round_key(Ωtot), _round_key(Ẑ12))
                        HZ[hkey] = Ẑ34
                        KZ[kkey] = (hkey, Ẑ23, c, d)
                        inner_plan = get!(plan, kkey, Dict{Any,Tuple{Float64,ComplexF64,ComplexF64}}())
                        _, _, prev = get(inner_plan, fkey, (Ωtot, Ẑ12, zero(ComplexF64)))
                        inner_plan[fkey] = (Ωtot, Ẑ12, prev + s_o * s_i * c1 * c2)
                    end
                end
            end
        end
    end

    nt = size(ρ1, 1)
    Hhat = Dict{Any,Vector{ComplexF64}}()
    for (hkey, Ẑ34) in HZ
        e, f, _ = hkey
        Hhat[hkey] = _exp_convolution(Ẑ34, ρ1[:, e, f], h)
    end
    Bacc = Dict{Any,Tuple{Float64,ComplexF64,Vector{ComplexF64}}}()
    for (kkey, (hkey, Ẑ23, c, d)) in KZ
        r = variant == :nz ? ones(ComplexF64, nt) : map(_safe_inv, ρ1[:, c, d])
        Khat = _exp_convolution(Ẑ23, r .* Hhat[hkey], h)
        for (fkey, (Ωtot, Ẑ12, coeff)) in plan[kkey]
            coeff == 0 && continue
            _, _, B = get!(Bacc, fkey, (Ωtot, Ẑ12, zeros(ComplexF64, nt)))
            B .+= coeff .* Khat
        end
    end
    return [_ConvTerm(k[1], k[2], v[1], v[2], one(ComplexF64), v[3]) for (k, v) in Bacc]
end

# ---------------------------------------------------------------------------
# Time stepping
# ---------------------------------------------------------------------------

# dρ_ab/dt = Σ_g coeff e^{iΩt} ∫_0^t e^{-Ẑ(t-τ)} B(τ) ρ_cd(τ) dτ + S_ab(t),
# trapezoid rule for ρ, exponential trapezoid for the history integrals.
function _solve_convolution(N::Int, ρ0::AbstractMatrix, h::Real, nsteps::Int,
                            terms::Vector{_ConvTerm}, source::Union{Nothing,Array{ComplexF64,3}})
    G = length(terms)
    dec = Vector{ComplexF64}(undef, G); w0 = similar(dec); w1 = similar(dec)
    for (g, term) in enumerate(terms)
        dec[g], w0[g], w1[g] = _exp_trapezoid(term.Ẑ, h)
    end
    Bval(g, i) = terms[g].B === nothing ? terms[g].coeff : terms[g].coeff * terms[g].B[i]
    N2 = N * N
    ρ = zeros(ComplexF64, nsteps + 1, N, N)
    ρ[1, :, :] .= ρ0
    x = vec(Matrix{ComplexF64}(ρ0))
    Y = zeros(ComplexF64, G)
    S(i) = source === nothing ? zeros(ComplexF64, N2) : vec(source[i, :, :])
    F = S(1)  # all history integrals vanish at t = 0
    Yhat = similar(Y)
    for i in 2:nsteps+1
        t = (i - 1) * h
        Fhat = S(i)
        Gmat = zeros(ComplexF64, N2, N2)
        for (g, term) in enumerate(terms)
            Yhat[g] = dec[g] * Y[g] + w0[g] * Bval(g, i - 1) * x[term.inp]
            phase = cis(term.Ω * t)
            Fhat[term.out] += phase * Yhat[g]
            Gmat[term.out, term.inp] -= phase * w1[g] * Bval(g, i)
        end
        A = LinearAlgebra.I(N2) + (h / 2) * Gmat
        x_new = A \ (x + (h / 2) * (F + Fhat))
        for (g, term) in enumerate(terms)
            Y[g] = Yhat[g] + w1[g] * Bval(g, i) * x_new[term.inp]
        end
        F = Fhat - Gmat * x_new
        x = x_new
        ρ[i, :, :] .= reshape(x, N, N)
    end
    return ρ
end

# Source term S_ab(t) = Σ_g coeff e^{iΩt} ∫_0^t e^{-Ẑ(t-τ)} B(τ) dτ
function _source_from_terms(N::Int, terms::Vector{_ConvTerm}, h::Real, nt::Int)
    S = zeros(ComplexF64, nt, N, N)
    for term in terms
        y = _exp_convolution(term.Ẑ, term.coeff .* term.B, h)
        a = mod1(term.out, N); b = div(term.out - 1, N) + 1
        for i in 1:nt
            S[i, a, b] += cis(term.Ω * (i - 1) * h) * y[i]
        end
    end
    return S
end

_time_grid(t_max, dt) = (n = round(Int, t_max / dt); (n, t_max / n, collect(range(0.0, t_max; length = n + 1))))

"""
    QME_cf_born(model, ρ0, t_max, dt)

Second-order time-nonlocal quantum master equation with the bath kept in equilibrium
(constant ansatz, thesis kernel ``M(t,\\tau;\\hat w^0)``) for an [`ExcitonModel`](@ref).
`ρ0` is the initial reduced density matrix in the exciton basis. Returns
`(t, ρ)` with `ρ[i, :, :]` the exciton-basis density matrix in the interaction picture.
"""
function QME_cf_born(model::ExcitonModel, ρ0::AbstractMatrix, t_max::Real, dt::Real)
    nsteps, h, t = _time_grid(t_max, dt)
    N = length(model.energies)
    ρ = _solve_convolution(N, ρ0, h, nsteps, _second_order_terms(model), nothing)
    return (t = t, ρ = ρ)
end

"""
    QME_cf_redfield(model, ρ0, t_max, dt)

Time-local second-order master equation (Redfield without secular or Markov
approximation, TCL2): ``\\dot\\rho(t) = -\\int_0^t d\\tau\\, M(t,\\tau;\\hat w^0)\\,\\rho(t)``.
"""
function QME_cf_redfield(model::ExcitonModel, ρ0::AbstractMatrix, t_max::Real, dt::Real)
    nsteps, h, t = _time_grid(t_max, dt)
    N = length(model.energies)
    terms = _second_order_terms(model)
    N2 = N * N
    R(τ) = begin
        M = zeros(ComplexF64, N2, N2)
        for term in terms
            M[term.out, term.inp] -= term.coeff * cis(term.Ω * τ) * _exp_integral(term.Ẑ, τ)
        end
        M
    end
    ρ = zeros(ComplexF64, nsteps + 1, N, N)
    ρ[1, :, :] .= ρ0
    x = vec(Matrix{ComplexF64}(ρ0))
    Rprev = R(0.0)
    for i in 2:nsteps+1
        Rnew = R((i - 1) * h)
        x = (LinearAlgebra.I(N2) + (h / 2) * Rnew) \ (x - (h / 2) * (Rprev * x))
        Rprev = Rnew
        ρ[i, :, :] .= reshape(x, N, N)
    end
    return (t = t, ρ = ρ)
end

"""
    QME_cf_iterative(model, ρ0, t_max, dt; correction=:nz, iterations=1)

Iterative quantum master equation with the first correction of the relative bath part
(thesis, section "Correction of memory kernel"), for finite or infinite harmonic baths.

1. ``\\rho^1`` solves the second-order QME with the equilibrium bath ([`QME_cf_born`](@ref)).
2. ``\\rho^{k+1}`` solves the QME whose memory kernel contains the first-order bath
   correction evaluated with ``\\rho^k``.

`correction` selects the form of the correction:
- `:thesis` — the kernel ``M(t_1,t_2;\\hat w^1)`` exactly as derived in the thesis.
- `:projected` — the thesis kernel with the quotient-rule term of the RBP equation,
  which removes the disconnected (secular) Wick pairing.
- `:nz` — `:projected` with ``\\rho_{cd}(t_2)/\\rho^k_{cd}(t_3) \\to 1``, i.e. the
  fourth-order Nakajima-Zwanzig term evaluated with ``\\rho^k`` (no division by ρ).

Returns `(t, ρ, iterates)` where `iterates[k]` is ``\\rho^k``.
"""
function QME_cf_iterative(model::ExcitonModel, ρ0::AbstractMatrix, t_max::Real, dt::Real;
                          correction::Symbol = :nz, iterations::Int = 1)
    nsteps, h, t = _time_grid(t_max, dt)
    N = length(model.energies)
    second = _second_order_terms(model)
    ρk = _solve_convolution(N, ρ0, h, nsteps, second, nothing)
    iterates = [ρk]
    for _ in 1:iterations
        fourth = _fourth_order_terms(model, ρk, h, correction)
        if correction == :nz
            S = _source_from_terms(N, fourth, h, nsteps + 1)
            ρk = _solve_convolution(N, ρ0, h, nsteps, second, S)
        else
            # M¹ = M⁰ - ΔM and dρ/dt = -∫ M¹ ρ, so the correction enters with +1
            ρk = _solve_convolution(N, ρ0, h, nsteps, vcat(second, fourth), nothing)
        end
        push!(iterates, ρk)
    end
    return (t = t, ρ = ρk, iterates = iterates)
end

"""
    exciton_to_schroedinger(model, t, ρ)

Convert interaction-picture exciton-basis density matrices `ρ[i, :, :]` at times `t`
to the Schrödinger picture.
"""
function exciton_to_schroedinger(model::ExcitonModel, t::AbstractVector, ρ::Array{ComplexF64,3})
    E = model.energies
    out = similar(ρ)
    for i in eachindex(t), a in axes(ρ, 2), b in axes(ρ, 3)
        out[i, a, b] = ρ[i, a, b] * cis(-(E[a] - E[b]) * t[i])
    end
    return out
end
