using OpenQuantumSystems, LinearAlgebra, Printf
const OQS = OpenQuantumSystems

# --- disconnected term identity: ΔM_thesis - ΔM_projected = -M⁰ ln(ρ¹_cd(τ)/ρ¹_cd(0))
λ, γ = 20.0, 100.0
Hsite = [12500.0+λ 100.0; 100.0 12700.0+λ]
bath = bath_correlation(DrudeLorentz(λ, γ), 300.0; n_matsubara = 4)
model = ExcitonModel(Hsite, [bath, bath])
ρ0 = site_to_exciton(model, ComplexF64[1 0; 0 0])
dt = 2.5e-4
born = QME_cf_born(model, ρ0, 0.2, dt)
Tth = OQS._fourth_order_terms(model, born.ρ, dt, :thesis)
Tpr = OQS._fourth_order_terms(model, born.ρ, dt, :projected)
T2 = OQS._second_order_terms(model)
kern(terms, i1, i2, withB) = begin
    M = zeros(ComplexF64, 4, 4); t1 = (i1 - 1) * dt; t2 = (i2 - 1) * dt
    for g in terms
        M[g.out, g.inp] += g.coeff * cis(g.Ω * t1) * exp(-g.Ẑ * (t1 - t2)) * (withB ? g.B[i2] : 1)
    end
    M
end
println("== disconnected part of the thesis kernel vs -M⁰ ln(ρ¹(τ)/ρ¹(0))")
for τ in (0.02, 0.1, 0.19)
    i2 = round(Int, τ / dt) + 1; i1 = i2 + 8
    D = kern(Tth, i1, i2, true) - kern(Tpr, i1, i2, true)
    M0 = -kern(T2, i1, i2, false)          # second-order terms store -M⁰
    pops = [1, 4]  # (1,1), (2,2): real ρ, no branch of the logarithm
    pred = zeros(ComplexF64, 4, 2)
    for (j, cd) in enumerate(pops)
        c = mod1(cd, 2)
        pred[:, j] = -M0[:, cd] * log(real(born.ρ[i2, c, c]) / real(born.ρ[1, c, c]))
    end
    @printf("  τ=%.2f  populations: |D - pred|/|D| = %.2e\n", τ, norm(D[:, pops] - pred) / norm(D[:, pops]))
end

# --- order of accuracy against exact finite-bath dynamics
function finite_case(S; nvib = 6, J = 50.0, T = 10.0, t_max = 0.1, dt = 2e-4)
    mode = Mode(; omega = 200.0, hr_factor = S)
    aggCore = AggregateCore([Molecule([mode], nvib, [0.0, 12500.0]), Molecule([mode], nvib, [0.0, 12700.0])])
    aggCore.coupling[2, 3] = J; aggCore.coupling[3, 2] = J
    agg = setup_aggregate(aggCore; vib_basis = GroundGround())
    ops, tools = agg.operators, agg.tools
    W_eq = equilibrium_bath_state(agg.core, tools, ops, T)
    ρsite = zeros(ComplexF64, 3, 3); ρsite[2, 2] = 0.7; ρsite[3, 3] = 0.3
    W0 = ad(ρsite, W_eq, agg.core, tools)
    model = ExcitonModel(aggCore; T = T)
    ρ0 = site_to_exciton(model, ρsite[2:3, 2:3])
    F = eigen(Hermitian(Matrix(ops.Ham.data))); F0 = eigen(Hermitian(Matrix(ops.Ham_0.data)))
    nsteps = round(Int, t_max / dt)
    exact = zeros(ComplexF64, nsteps + 1, 2, 2)
    for i in 1:nsteps+1
        t = (i - 1) * dt
        U = F.vectors * Diagonal(cis.(-F.values * t)) * F.vectors'
        U0 = F0.vectors * Diagonal(cis.(-F0.values * t)) * F0.vectors'
        W = U0' * U * W0 * U' * U0
        exact[i, :, :] = site_to_exciton(model, trace_bath(W, agg.core, ops, tools)[2:3, 2:3])
    end
    change = maximum(norm(exact[i, :, :] - ρ0) for i in 1:nsteps+1)
    err(ρ) = maximum(norm(ρ[i, :, :] - exact[i, :, :]) for i in 1:nsteps+1) / change
    res = Dict{String,Float64}()
    res["Redfield"] = err(QME_cf_redfield(model, ρ0, t_max, dt).ρ)
    res["Born (M⁰)"] = err(QME_cf_born(model, ρ0, t_max, dt).ρ)
    for v in (:thesis, :projected, :nz)
        res["IQME $v"] = err(QME_cf_iterative(model, ρ0, t_max, dt; correction = v).ρ)
    end
    return res, change
end

println("== relative error vs exact (finite bath, nvib = 6, t ≤ 0.1)")
Ss = (0.00625, 0.0125, 0.025, 0.05, 0.1)
results = Dict(S => finite_case(S) for S in Ss)
methods = ["Redfield", "Born (M⁰)", "IQME thesis", "IQME projected", "IQME nz"]
@printf("%-16s", "S"); foreach(S -> @printf("%12.4f", S), Ss); println("   slope(two weakest)")
for m in methods
    @printf("%-16s", m)
    e = [results[S][1][m] for S in Ss]
    foreach(x -> @printf("%12.2e", x), e)
    slope = (log(e[2]) - log(e[1])) / (log(Ss[2]) - log(Ss[1]))
    @printf("   %5.2f\n", slope)
end
