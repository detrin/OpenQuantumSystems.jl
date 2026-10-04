# Shared fixtures for master-equation tests: a small dimer and its exact dynamics.
using OpenQuantumSystems
using LinearAlgebra
import DelayDiffEq
import OrdinaryDiffEq

# Dimer with one 200 cm⁻¹ mode per molecule, 2 vibrational levels per mode,
# laser-excited from a 10 K thermal ground state.
function qme_test_dimer(S; nvib = 2, J = 50.0)
    mode = Mode(; omega = 200.0, hr_factor = S)
    aggCore = AggregateCore([
        Molecule([mode], nvib, [0.0, 12500.0]),
        Molecule([mode], nvib, [0.0, 12700.0]),
    ])
    aggCore.coupling[2, 3] = J
    aggCore.coupling[3, 2] = J
    agg = setup_aggregate(aggCore; vib_basis = GroundGround())
    W0, _, _ = ultrafast_laser_excitation(10.0, [0.0, 0.7, 0.3], agg)
    return agg, W0
end

# A third of a vibrational period: long enough for visible bath-induced dynamics.
const QME_TEST_TSPAN = collect(range(0.0, 0.01; length = 6))

qme_test_kwargs(tol = 1e-4) = (
    reltol = tol, abstol = tol / 100, int_reltol = tol, int_abstol = tol / 100,
    alg = DelayDiffEq.MethodOfSteps(OrdinaryDiffEq.Tsit5()),
)

# Exact W(t) by diagonalisation; `picture` is :schrodinger or :interaction,
# `reduced` traces over the bath.
function qme_exact(agg, W0, tspan; picture = :interaction, reduced = true)
    ops = agg.operators
    map(tspan) do t
        U = evolution_operator(ops.Ham, t)
        W = U * W0 * U'
        if picture == :interaction
            U0 = evolution_operator(ops.Ham_0, t)
            W = U0' * W * U0
        end
        reduced ? trace_bath(W, agg.core, ops, agg.tools; vib_basis = ops.vib_basis).data : W.data
    end
end

_qme_data(x) = x isa AbstractMatrix ? x : x.data

# Largest deviation from the reference, relative to how much the reference changes.
function qme_relative_error(xs, refs)
    dev = maximum(norm(_qme_data(x) - r) for (x, r) in zip(xs, refs))
    change = maximum(norm(r - refs[1]) for r in refs)
    return dev / change
end

qme_trace_drift(xs) = maximum(abs(tr(_qme_data(x)) - tr(_qme_data(xs[1]))) for x in xs)
qme_hermiticity(xs) = maximum(norm(_qme_data(x) - _qme_data(x)') for x in xs)
