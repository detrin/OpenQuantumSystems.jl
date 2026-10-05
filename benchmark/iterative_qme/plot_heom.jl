using Serialization, Plots, LaTeXStrings, LinearAlgebra

store = deserialize(joinpath(@__DIR__, "heom_benchmark.jls"))
outdir = ARGS[1]
methods = ["Redfield", "Born (M⁰)", "IQME thesis", "IQME projected", "IQME nz"]
labels = Dict("Redfield" => "Redfield (TCL2)", "Born (M⁰)" => "Born QME (M⁰)",
              "IQME thesis" => "IQME, thesis kernel", "IQME projected" => "IQME, projected",
              "IQME nz" => "IQME, NZ form")
styles = Dict("Redfield" => (:dash, :blue), "Born (M⁰)" => (:dash, :cyan),
              "IQME thesis" => (:dashdot, :orange), "IQME projected" => (:dot, :green),
              "IQME nz" => (:solid, :red))

for λ in (10.0, 40.0)
    t, ref, res = store[λ]
    plt = plot(size = (700, 450), dpi = 150, legend = :topright)
    plot!(plt, t, [real(r[1, 1]) for r in ref], label = "HEOM (exact)", lw = 3, color = :black)
    for m in methods
        ls, col = styles[m]
        plot!(plt, t, [real(r[1, 1]) for r in res[m]], label = labels[m], ls = ls, lw = 1.8, color = col)
    end
    xlabel!(plt, L"t\;\;(\mathrm{cm}^{-1}\;\mathrm{units})")
    ylabel!(plt, L"\rho_{11}(t)")
    title!(plt, "Drude-Lorentz bath, λ = $(Int(λ)) cm⁻¹, γ = 100 cm⁻¹, 300 K")
    savefig(plt, joinpath(outdir, "iqme_heom_lambda$(Int(λ)).png"))
end

λs = sort(collect(keys(store)))
relerr(xs, refs) = maximum(norm(x - r) for (x, r) in zip(xs, refs)) / maximum(norm(r - refs[1]) for r in refs)
plt = plot(size = (650, 450), dpi = 150, xscale = :log10, yscale = :log10, legend = :bottomright)
for m in methods
    ls, col = styles[m]
    errs = [relerr(store[λ][3][m], store[λ][2]) for λ in λs]
    plot!(plt, λs, errs, label = labels[m], ls = ls, lw = 1.8, color = col, marker = :circle)
end
xlabel!(plt, "reorganization energy λ (cm⁻¹)")
ylabel!(plt, "relative error vs HEOM")
title!(plt, "Accuracy against HEOM, t ≤ 0.3")
savefig(plt, joinpath(outdir, "iqme_heom_error.png"))
println("plots written")
