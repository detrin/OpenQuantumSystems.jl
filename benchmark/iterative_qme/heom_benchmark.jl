using HierarchicalEOM, OpenQuantumSystems, LinearAlgebra, Printf, Serialization
const OQS = OpenQuantumSystems

const E1, E2, J = 12500.0, 12700.0, 100.0
const γ, T, NMATS = 100.0, 300.0, 1
const kT = OQS.BOLTZMANN_CM * T

# --- 1. bath convention: HierarchicalEOM vs ours ---------------------------------
let λ = 20.0
    P1 = Qobj([1.0 0.0; 0.0 0.0])
    bath = Boson_DrudeLorentz_Matsubara(P1, λ, γ, kT, NMATS)
    ts = [0.001, 0.005, 0.02]
    C_heom = HierarchicalEOM.correlation_function(bath, ts)
    C_ours = bath_correlation(DrudeLorentz(λ, γ), T; n_matsubara = NMATS).(ts)
    println("bath convention check, max |ΔC|/|C| = ", maximum(abs.(C_heom .- C_ours) ./ abs.(C_ours))); flush(stdout)
end

function heom_reference(λ; tier = 8, t_max = 0.3, dt = 5e-4)
    Eshift = (E1 + E2) / 2
    Hs = Qobj([E1+λ-Eshift J; J E2+λ-Eshift])
    P1 = Qobj([1.0 0.0; 0.0 0.0]); P2 = Qobj([0.0 0.0; 0.0 1.0])
    baths = [Boson_DrudeLorentz_Matsubara(P1, λ, γ, kT, NMATS),
             Boson_DrudeLorentz_Matsubara(P2, λ, γ, kT, NMATS)]
    M = M_Boson(Hs, tier, baths; verbose = false)
    println("  HEOM λ=$λ tier=$tier dim=$(size(M, 1))"); flush(stdout)
    ρ0 = Qobj(ComplexF64[1 0; 0 0])
    tlist = collect(0:dt:t_max)
    sol = HEOMsolve(M, ρ0, tlist; reltol = 1e-9, abstol = 1e-11, progress_bar = Val(false))
    ρs = [Matrix(getRho(ados).data) for ados in sol.ados]
    return tlist, ρs
end

function ours(λ; t_max = 0.3, dt = 5e-4)
    bath = bath_correlation(DrudeLorentz(λ, γ), T; n_matsubara = NMATS)
    model = ExcitonModel([E1+λ J; J E2+λ], [bath, bath])
    ρ0 = site_to_exciton(model, ComplexF64[1 0; 0 0])
    to_site(res) = begin
        ρS = exciton_to_schroedinger(model, res.t, res.ρ)
        [exciton_to_site(model, ρS[i, :, :]) for i in eachindex(res.t)]
    end
    out = Dict{String,Vector{Matrix{ComplexF64}}}()
    out["Redfield"] = to_site(QME_cf_redfield(model, ρ0, t_max, dt))
    out["Born (M⁰)"] = to_site(QME_cf_born(model, ρ0, t_max, dt))
    for v in (:thesis, :projected, :nz)
        out["IQME $v"] = to_site(QME_cf_iterative(model, ρ0, t_max, dt; correction = v))
    end
    return out
end

# HEOM phases are relative to the mean site energy; compare in a frame where that
# global phase is irrelevant (it cancels in ρ for a fixed electronic manifold).
relerr(xs, refs) = maximum(norm(x - r) for (x, r) in zip(xs, refs)) / maximum(norm(r - refs[1]) for r in refs)

λs = (2.5, 5.0, 10.0, 20.0, 40.0)
# tier convergence for the strongest coupling
_, r8 = heom_reference(40.0; tier = 8)
_, r12 = heom_reference(40.0; tier = 12)
println("HEOM tier 8 vs 12 at λ = 40: rel diff = ", relerr(r8, r12)); flush(stdout)

methods = ["Redfield", "Born (M⁰)", "IQME thesis", "IQME projected", "IQME nz"]
table = Dict{String,Vector{Float64}}(m => Float64[] for m in methods)
store = Dict{Float64,Any}()
for λ in λs
    t, ref = heom_reference(λ; tier = 10)
    res = ours(λ)
    store[λ] = (t, ref, res)
    for m in methods
        push!(table[m], relerr(res[m], ref))
    end
    println("λ = $λ done"); flush(stdout)
end
serialize(joinpath(@__DIR__, "heom_benchmark.jls"), store)

@printf("%-16s", "λ (cm⁻¹)"); foreach(l -> @printf("%10.1f", l), λs); println("   slope(two weakest)")
for m in methods
    e = table[m]
    @printf("%-16s", m); foreach(x -> @printf("%10.2e", x), e)
    @printf("   %5.2f\n", (log(e[2]) - log(e[1])) / (log(λs[2]) - log(λs[1])))
end
