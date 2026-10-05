using OpenQuantumSystems, LinearAlgebra
const OQS = OpenQuantumSystems

mode = Mode(; omega = 200.0, hr_factor = 0.05)
aggCore = AggregateCore([Molecule([mode], 4, [0.0, 12500.0]), Molecule([mode], 4, [0.0, 12700.0])])
aggCore.coupling[2, 3] = 50.0; aggCore.coupling[3, 2] = 50.0
agg = setup_aggregate(aggCore; vib_basis = GroundGround())
ops, tools = agg.operators, agg.tools
im_ = tools.indicesMap
W_eq = equilibrium_bath_state(agg.core, tools, ops, 10.0)
w_eq = W_eq[im_[1], im_[1]]
ρsite = zeros(ComplexF64, 3, 3); ρsite[2, 2] = 0.7; ρsite[3, 3] = 0.3; ρsite[2, 3] = ρsite[3, 2] = 0.2
W0 = ad(ρsite, W_eq, agg.core, tools)
H, H0 = ops.Ham.data, ops.Ham_0.data
WI(t) = (U = exp(-im * H * t); U0 = exp(-im * H0 * t); U0' * U * W0 * U' * U0)
block(W, c, d) = W[im_[c], im_[d]]

println("== C1: thesis RBP equation dŵ/dt = Ẇ/ρ vs exact dŵ/dt")
for t in (0.005, 0.02, 0.05)
    δ = 1e-6
    Wm, W, Wp = WI(t - δ), WI(t), WI(t + δ)
    for (c, d) in ((2, 2), (2, 3))
        ρm, ρ, ρp = tr(block(Wm, c, d)), tr(block(W, c, d)), tr(block(Wp, c, d))
        ŵdot = (block(Wp, c, d) / ρp - block(Wm, c, d) / ρm) / (2δ)
        Wdot = (block(Wp, c, d) - block(Wm, c, d)) / (2δ)
        ρdot = (ρp - ρm) / (2δ)
        thesis = Wdot / ρ
        quotient = thesis - (block(W, c, d) / ρ) * ρdot / ρ
        println("  t=$t (c,d)=($c,$d): |ŵ̇ - Ẇ/ρ| / |ŵ̇| = ", round(norm(ŵdot - thesis) / norm(ŵdot); sigdigits = 3),
                "   with -ŵρ̇/ρ term: ", round(norm(ŵdot - quotient) / norm(ŵdot); sigdigits = 3),
                "   tr_B of thesis ŵ̇ = ", round(abs(tr(thesis)); sigdigits = 3))
    end
end

println("== C2: odd number of H_I under tr_B vanishes (thesis: K_w term dropped)")
HI(t) = (U0 = exp(-im * H0 * t); U0' * ops.Ham_I.data * U0)
X = zeros(ComplexF64, size(H)); X[im_[2], im_[2]] .= w_eq
Y = HI(0.004) * HI(0.002) * HI(0.001) * X
println("  max |tr_B <a|H H H (|2><2| ⊗ w)|b>| = ",
    maximum(abs(tr(block(Y, a, b))) for a in 1:3, b in 1:3), "   (scale of two H: ",
    maximum(abs(tr(block(HI(0.004) * HI(0.002) * X, a, b))) for a in 1:3, b in 1:3), ")")

println("== C3: secular growth of the thesis correction (Drude-Lorentz bath)")
λ, γ = 20.0, 100.0
Hsite = [12500.0+λ 100.0; 100.0 12700.0+λ]
bath = bath_correlation(DrudeLorentz(λ, γ), 300.0; n_matsubara = 4)
model = ExcitonModel(Hsite, [bath, bath])
ρ0 = site_to_exciton(model, ComplexF64[1 0; 0 0])
dt = 5e-4
born = QME_cf_born(model, ρ0, 0.4, dt)
lag = round(Int, 0.002 / dt)
for variant in (:thesis, :projected)
    terms = OQS._fourth_order_terms(model, born.ρ, dt, variant)
    vals = Float64[]
    for τ in (0.05, 0.1, 0.2, 0.39)
        i2 = round(Int, τ / dt) + 1; i1 = i2 + lag; t1 = (i1 - 1) * dt
        ΔM = zeros(ComplexF64, 4, 4)
        for g in terms
            ΔM[g.out, g.inp] += cis(g.Ω * t1) * exp(-g.Ẑ * (t1 - (i2 - 1) * dt)) * g.B[i2]
        end
        push!(vals, maximum(abs.(ΔM[[1, 4], [1, 4]])))  # population block
    end
    println("  $variant |ΔM(τ + 0.002, τ)| for τ = 0.05, 0.1, 0.2, 0.39: ", round.(vals; sigdigits = 3))
end
