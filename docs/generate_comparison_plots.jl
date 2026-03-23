using OpenQuantumSystems
using LinearAlgebra
using Plots
using LaTeXStrings

outdir = joinpath(@__DIR__, "src", "assets")
mkpath(outdir)

"""
    exciton_initial_W0(agg, coeffs, exciton_idx, T)

Prepare initial density matrix W0 for a pure exciton eigenstate.
The electronic part is |ψ⟩⟨ψ| where |ψ⟩ = Σ_k c_{m,k} |k⟩ (exciton m),
tensored with the thermal bath at temperature T.
"""
function exciton_initial_W0(agg, coeffs, exciton_idx::Int, T::Real)
    aggCore, aggTools, aggOps = agg.core, agg.tools, agg.operators
    N = aggCore.molCount
    elLen = N + 1
    indicesMap = aggTools.indicesMap
    bSize = aggTools.bSize

    # Get thermal bath density matrix (from ground state block)
    a1, a2 = indicesMap[1][1], indicesMap[1][end]
    H_0_g = real.(aggOps.Ham_0.data[a1:a2, a1:a2])
    rho_bath = exp(-H_0_g / (BOLTZMANN_CM * T))
    rho_bath ./= tr(rho_bath)
    vib_dim = a2 - a1 + 1

    # Build full W0 = |ψ⟩⟨ψ| ⊗ rho_bath
    # where |ψ⟩ = Σ_k c_{m,k} |k+1⟩ (k+1 because index 1 is ground state)
    W0_data = zeros(ComplexF64, bSize, bSize)
    m = exciton_idx
    for k in 1:N, l in 1:N
        # Electronic indices k+1 and l+1 (1-indexed, ground=1)
        ak1, ak2 = indicesMap[k+1][1], indicesMap[k+1][end]
        al1, al2 = indicesMap[l+1][1], indicesMap[l+1][end]
        W0_data[ak1:ak2, al1:al2] .= coeffs[m, k] * coeffs[m, l] .* rho_bath
    end

    basis = agg.operators.Ham.basis_l
    return DenseOperator(basis, basis, W0_data)
end

function run_comparison(S; J=100.0, T=300.0, nvib=5, t_max=0.5, dt=0.001,
                        include_corrected=true)
    mode = Mode(; omega=200.0, hr_factor=S)
    mol1 = Molecule([mode], nvib, [0.0, 12500.0])
    mol2 = Molecule([mode], nvib, [0.0, 12700.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = J
    aggCore.coupling[3, 2] = J

    agg = setup_aggregate(aggCore)
    Ham = agg.operators.Ham
    N = aggCore.molCount
    elLen = N + 1
    tspan = collect(0.0:dt:t_max)

    aggTools = AggregateTools(aggCore)
    energies, coeffs = exciton_basis(aggCore, aggTools)

    results = Dict{String, Vector{Float64}}()
    results["tspan"] = tspan

    # --- Initial condition: upper exciton eigenstate (m=2) ---
    println("  Preparing upper exciton eigenstate...")
    W0 = exciton_initial_W0(agg, coeffs, 2, T)

    # Verify initial exciton population
    rho_site_0 = trace_bath(W0.data, agg)
    p2_0 = sum(coeffs[2, k] * coeffs[2, l] * rho_site_0[k+1, l+1] for k in 1:2, l in 1:2)
    println("  Initial upper exciton pop: ", real(p2_0))

    # --- Exact evolution ---
    println("  Exact evolution...")
    _, W_t = evolution_approximate(W0, tspan, Ham)
    rho_exact_upper = zeros(Float64, length(tspan))
    for t_i in 1:length(tspan)
        rho_site = trace_bath(W_t[t_i].data, agg)
        p2 = sum(coeffs[2, k] * coeffs[2, l] * rho_site[k+1, l+1] for k in 1:2, l in 1:2)
        rho_exact_upper[t_i] = real(p2)
    end
    results["Exact"] = rho_exact_upper

    # --- Modified Redfield (exciton basis natively) ---
    println("  Modified Redfield...")
    _, pop_mr = modified_redfield_dynamics(aggCore, [0.0, 1.0], tspan; T=T)
    results["Modified Redfield"] = pop_mr[:, 2]

    # --- Förster (site basis rates, transform to exciton-basis rate matrix) ---
    println("  Förster...")
    K_site = forster_rate_matrix(aggCore; T=T, sigma=30.0)
    # Transform site rates to exciton basis: K_exc = C * K_site * C^{-1}
    # where C[m,k] = |c_{m,k}|^2 (population transform)
    C_pop = [coeffs[m, k]^2 for m in 1:N, k in 1:N]
    K_exc = C_pop * K_site * inv(C_pop)
    p = [0.0, 1.0]  # start in upper exciton
    pop_forster_exc = zeros(Float64, length(tspan), 2)
    pop_forster_exc[1, :] .= p
    for i in 2:length(tspan)
        ddt = tspan[i] - tspan[i-1]
        dp = K_exc * p
        p .= p .+ ddt .* dp
        clamp!(p, 0.0, 1.0)
        s = sum(p); s > 0 && (p ./= s)
        pop_forster_exc[i, :] .= p
    end
    results["Förster"] = pop_forster_exc[:, 2]

    # --- Corrected kernel (exciton basis natively) ---
    if include_corrected
        println("  Corrected kernel (iterative)...")
        t_corr, pop_corr, all_iters = corrected_qme_rdm(aggCore, [0.0, 1.0], tspan;
            T=T, t_ref=0.1, max_iter=2, rtol=1e-2, atol=1e-4)
        results["Corrected (0th)"] = all_iters[1][:, 2]
        results["Corrected (1st)"] = all_iters[2][:, 2]
    end

    return results
end

function make_plot(results, S; include_corrected=true)
    tspan = results["tspan"]
    plt = plot(size=(700, 450), dpi=150)

    plot!(plt, tspan, results["Exact"], label="Exact", lw=2.5, color=:black)
    plot!(plt, tspan, results["Modified Redfield"],
        label="Mod. Redfield", ls=:dot, lw=2, color=:red)
    plot!(plt, tspan, results["Förster"],
        label="Förster", ls=:dot, lw=2, color=:purple)

    if include_corrected && haskey(results, "Corrected (0th)")
        plot!(plt, tspan, results["Corrected (0th)"],
            label="Corrected (0th)", ls=:dashdot, lw=1.5, color=:orange)
        plot!(plt, tspan, results["Corrected (1st)"],
            label="Corrected (1st)", ls=:dashdot, lw=1.5, color=:green)
    end

    xlabel!(plt, L"t\;\;(\mathrm{cm}^{-1}\;\mathrm{units})")
    ylabel!(plt, L"\rho^\mathrm{exc}_{22}(t)")
    title!(plt, "Upper exciton population — S = $S")
    ylims!(plt, 0.0, 1.05)
    return plt
end

for (S, include_corr) in [(0.01, true), (0.05, true), (0.1, true)]
    println("=== S = $S ===")
    @time results = run_comparison(S; include_corrected=include_corr)
    plt = make_plot(results, S; include_corrected=include_corr)
    fname = joinpath(outdir, "comparison_S$(replace(string(S), "." => "p")).png")
    savefig(plt, fname)
    println("  Saved: $fname")
end

println("Done!")
