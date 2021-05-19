
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
function liouvilleVonNeumann(
    rho0::T,
    tspan::Array,
    H::AbstractOperator{B,B};
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.DP5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    # tmp = copy(rho0)
    dliouvilleVonNeumann_(t, rho::T, drho::T) = dliouvilleVonNeumann(rho, H, drho)
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate(
        tspan_,
        dliouvilleVonNeumann_,
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

function dliouvilleVonNeumann(
    rho::T,
    H::AbstractOperator{B,B},
    drho::T,
) where {B<:Basis,T<:Operator{B,B}}
    QuantumOpticsBase.mul!(drho, H, rho, -eltype(rho)(im), zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho, rho, H, eltype(rho)(im), one(eltype(rho)))
    return drho
end
