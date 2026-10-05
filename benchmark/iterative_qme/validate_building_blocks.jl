using OpenQuantumSystems, LinearAlgebra, QuadGK
const OQS = OpenQuantumSystems
import DelayDiffEq, OrdinaryDiffEq

println("== V1 correlation functions")
sd = SpectralDensity([200.0, 350.0], [0.05, 0.02])
C = bath_correlation(sd, 300.0)
println("discrete vs analytic_correlation_nn: ",
    maximum(abs(C(τ) - analytic_correlation_nn(sd, τ, 300.0)) for τ in (-0.02, 0.0, 0.003, 0.05)))
dl = DrudeLorentz(20.0, 100.0)
T = 300.0
β = 1 / (OQS.BOLTZMANN_CM * T)
Cnum(τ) = quadgk(ω -> OQS.spectral_density(dl, ω) / π * (coth(β * ω / 2) * cos(ω * τ) - im * sin(ω * τ)), 0, Inf; rtol = 1e-10)[1]
Cdl = bath_correlation(dl, T; n_matsubara = 2000)
println("Drude-Lorentz vs numerical integral (rel): ",
    maximum(abs(Cdl(τ) - Cnum(τ)) / abs(Cnum(τ)) for τ in (0.002, 0.01, 0.05)))
println("C(-τ) = conj C(τ): ", Cdl(-0.01) ≈ conj(Cdl(0.01)))

# finite dimer
function dimer(nvib; S = 0.05, J = 50.0)
    mode = Mode(; omega = 200.0, hr_factor = S)
    aggCore = AggregateCore([Molecule([mode], nvib, [0.0, 12500.0]), Molecule([mode], nvib, [0.0, 12700.0])])
    aggCore.coupling[2, 3] = J; aggCore.coupling[3, 2] = J
    aggCore
end
aggCore = dimer(3)
Tlow = 10.0
model = ExcitonModel(aggCore; T = Tlow)

evalM(terms, t1, t2, N) = begin
    M = zeros(ComplexF64, N, N, N, N)
    for g in terms
        a = mod1(g.out, N); b = div(g.out - 1, N) + 1; c = mod1(g.inp, N); d = div(g.inp - 1, N) + 1
        M[a, b, c, d] -= g.coeff * cis(g.Ω * t1) * exp(-g.Ẑ * (t1 - t2))
    end
    M
end

println("== V2 second-order kernel vs zeroth_order_memory_kernel_cf")
terms2 = OQS._second_order_terms(model)
for (t1, t2) in ((0.003, 0.001), (0.02, 0.005))
    Mref = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = Tlow)[2:3, 2:3, 2:3, 2:3]
    M = evalM(terms2, t1, t2, 2)
    println("  ($t1,$t2) rel err = ", maximum(abs.(M .- Mref)) / maximum(abs.(Mref)))
end

println("== V3 fourth-order (thesis, constant ρ) vs first_order_memory_kernel_cf")
R = ComplexF64[0.6 0.2+0.1im; 0.2-0.1im 0.4]
h = 2e-6
t2max = 0.004
nt = round(Int, t2max / h) + 1
ρc = zeros(ComplexF64, nt, 2, 2); for i in 1:nt; ρc[i, :, :] .= R; end
terms4 = OQS._fourth_order_terms(model, ρc, h, :thesis)
println("  number of convolution terms: ", length(terms4))
Rfull = zeros(ComplexF64, 3, 3); Rfull[1, 1] = 1; Rfull[2:3, 2:3] .= R
for (t1, i2) in ((0.005, nt), (0.006, div(nt, 2)))
    t2 = (i2 - 1) * h
    ΔM = zeros(ComplexF64, 2, 2, 2, 2)
    for g in terms4
        a = mod1(g.out, 2); b = div(g.out - 1, 2) + 1; c = mod1(g.inp, 2); d = div(g.inp - 1, 2) + 1
        ΔM[a, b, c, d] += cis(g.Ω * t1) * exp(-g.Ẑ * (t1 - t2)) * g.B[i2]
    end
    M0 = zeroth_order_memory_kernel_cf(t1, t2, aggCore; T = Tlow)[2:3, 2:3, 2:3, 2:3]
    M1 = first_order_memory_kernel_cf(t1, t2, aggCore, Rfull; T = Tlow, rtol = 1e-10, atol = 1e-12)[2:3, 2:3, 2:3, 2:3]
    ref = M0 .- M1
    println("  t1=$t1 t2=$t2 rel err = ", maximum(abs.(ΔM .- ref)) / maximum(abs.(ref)))
end

println("== V4 solvers vs finite-basis QME (same equations)")
agg = setup_aggregate(aggCore; vib_basis = GroundGround())
ops, tools = agg.operators, agg.tools
W_eq = equilibrium_bath_state(agg.core, tools, ops, Tlow)
ρsite = zeros(ComplexF64, 3, 3); ρsite[2, 2] = 1.0
W0 = DenseOperator(tools.basis, tools.basis, ad(ρsite, W_eq, agg.core, tools))
ρ0_exc = site_to_exciton(model, ρsite[2:3, 2:3])
t_max, dt = 0.03, 1e-5
kw = (reltol = 1e-8, abstol = 1e-10, int_reltol = 1e-8, int_abstol = 1e-10,
      alg = DelayDiffEq.MethodOfSteps(OrdinaryDiffEq.Tsit5()))
tspan = collect(range(0, t_max; length = 7))
to_exc(r) = site_to_exciton(model, r.data[2:3, 2:3])
grid_at(res, t) = res.ρ[round(Int, t / (res.t[2] - res.t[1])) + 1, :, :]

born = QME_cf_born(model, ρ0_exc, t_max, dt)
_, ref = QME_sI_ansatz(W0, tspan, agg; ansatz = :const_int, kw...)
change = maximum(norm(to_exc(r) - ρ0_exc) for r in ref)
println("  Born vs QME_sI_ansatz(:const_int): rel err = ",
    maximum(norm(grid_at(born, t) - to_exc(r)) for (t, r) in zip(tspan, ref)) / change)
red = QME_cf_redfield(model, ρ0_exc, t_max, dt)
_, refR = QME_sI_Redfield(W0, tspan, agg; kw...)
println("  Redfield vs QME_sI_Redfield: rel err = ",
    maximum(norm(grid_at(red, t) - to_exc(r)) for (t, r) in zip(tspan, refR)) / change)
