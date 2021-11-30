
import OrdinaryDiffEq
import QuantumOpticsBase

"""
    liouvilleVonNeumann(rho0, tspan, H; 
    \treltol=1.0e-12, abstol=1.0e-12, fout=nothing, alg=OrdinaryDiffEq.Tsit5())

Integrate Liouville-von Neumann equation to evolve states or compute propagators

``\\frac{d}{d t} \\rho(t) = - \\frac{i}{\\hbar} [ \\hat H, \\rho(t) ], \\quad \\hbar = 1``.

# Arguments
* `rho0`: Initial state vector (can be a bra or a ket) or initial propagator.
* `tspan`: Vector specifying the points of time for which output should be displayed.
* `H`: Arbitrary operator specifying the Hamiltonian.
* `reltol`: Relative tolerance for OrdinaryDiffEq solver and its inner states.
* `abstol`: Absolute tolerance for OrdinaryDiffEq solver and its inner states.
* `fout=nothing`: If given, this function `fout(t, rho)` is called every time
        an output should be displayed. ATTENTION: The state `rho` is neither
        normalized nor permanent! It is still in use by the ode solver and
        therefore must not be changed.
* `alg`: Algorithm with which OrdinaryDiffEq will solve LvN equation.
"""
function LvN_SS(
    W0::T,
    tspan::Array,
    agg::Aggregate;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.Tsit5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    p = (agg.operators, eltype(W0))
    dLvN_(t, W::T, dW::T) = dLvN_SS(t, W, dW, p)

    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = W0.data
    state = T(W0.basis_l, W0.basis_r, W0.data)
    dstate = T(W0.basis_l, W0.basis_r, W0.data)

    OpenQuantumSystems.integrate(
        tspan_,
        dLvN_,
        x0,
        state,
        dstate,
        fout;
        reltol = reltol,
        abstol = abstol,
        alg = alg,
        kwargs...,
    )
end

function dLvN_SS(
    t::AbstractFloat,
    rho::T,
    drho::T,
    p
) where {B<:Basis,T<:Operator{B,B}}
    aggOperators, elementtype = p
    Ham = aggOperators.Ham.data
    drho.data[:, :] = -elementtype(im) * (Ham * rho.data - rho.data * Ham)
    return drho
end
