import OrdinaryDiffEq

"""
    schroedinger(psi0, tspan, H; 
    \treltol=1.0e-12, abstol=1.0e-12, fout=nothing, alg=OrdinaryDiffEq.Tsit5())

Integrate Schroedinger equation to evolve states or compute propagators

`` \\frac{d}{d t}\\vert\\psi(t)\\rangle = - \\frac{i}{\\hbar} \\hat H\\vert\\psi(t)\\rangle, \\quad \\hbar = 1``.

# Arguments
* `psi0`: Initial state vector (can be a bra or a ket) or initial propagator.
* `tspan`: Vector specifying the points of time for which output should be displayed.
* `H`: Arbitrary operator specifying the Hamiltonian.
* `reltol`: Relative tolerance for OrdinaryDiffEq solver and its inner states.
* `abstol`: Absolute tolerance for OrdinaryDiffEq solver and its inner states.
* `fout=nothing`: If given, this function `fout(t, psi)` is called every time
        an output should be displayed. ATTENTION: The state `psi` is neither
        normalized nor permanent! It is still in use by the ode solver and
        therefore must not be changed.
* `alg`: Algorithm with which OrdinaryDiffEq will solve Schroedinger equation.
"""
function schroedinger(
    psi0::T,
    tspan::Array,
    H::AbstractOperator{B,B};
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.Tsit5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Union{AbstractOperator{B,B},StateVector{B}}}
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    dschroedinger_(t, psi::T, dpsi::T, p) = dschroedinger(psi, H, dpsi)
    x0 = psi0.data
    state = copy(psi0)
    dstate = copy(psi0)
    OpenQuantumSystems.integrate(
        tspan_,
        dschroedinger_,
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


"""
    schroedinger_dynamic(psi0, tspan, f; 
    \treltol=1.0e-12, abstol=1.0e-12, fout=nothing, alg=OrdinaryDiffEq.Tsit5())

Integrate time-dependent Schroedinger equation to evolve states or compute propagators

``\\frac{d}{d t}\\vert\\psi(t)\\rangle = - \\frac{i}{\\hbar} \\hat H(t)\\vert\\rho(t)\\rangle, \\quad \\hbar = 1``.

# Arguments
* `psi0`: Initial state vector (can be a bra or a ket) or initial propagator.
* `tspan`: Vector specifying the points of time for which output should be displayed.
* `f`: Function `f(t, psi) -> H` returning the time and or state dependent Hamiltonian.
* `reltol`: Relative tolerance for OrdinaryDiffEq solver and its inner states.
* `abstol`: Absolute tolerance for OrdinaryDiffEq solver and its inner states.
* `fout=nothing`: If given, this function `fout(t, psi)` is called every time
        an output should be displayed. ATTENTION: The state `psi` is neither
        normalized nor permanent! It is still in use by the ode solver and
        therefore must not be changed.
* `alg`: Algorithm with which OrdinaryDiffEq will solve Schroedinger equation.
"""
function schroedinger_dynamic(
    psi0::T,
    tspan::Array,
    f::Function;
    reltol::Float64 = 1.0e-6,
    abstol::Float64 = 1.0e-8,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.DP5(),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {T<:Union{StateVector,AbstractOperator}}
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    dschroedinger_(t, psi::T, dpsi::T, p) = dschroedinger_dynamic(t, psi, f, dpsi)
    x0 = psi0.data
    state = copy(psi0)
    dstate = copy(psi0)
    OpenQuantumSystems.integrate(
        tspan_,
        dschroedinger_,
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


recast!(x::D, psi::StateVector{B,D}) where {B<:Basis,D} = (psi.data = x);
recast!(psi::StateVector{B,D}, x::D) where {B<:Basis,D} = nothing


function dschroedinger(
    psi::Union{Ket{B},AbstractOperator{B,B}},
    H::AbstractOperator{B,B},
    dpsi::Union{Ket{B},AbstractOperator{B,B}},
) where {B<:Basis}
    OpenQuantumSystems.mul!(dpsi, H, psi, eltype(psi)(-im), zero(eltype(psi)))
    return dpsi
end

function dschroedinger(psi::Bra{B}, H::AbstractOperator{B,B}, dpsi::Bra{B}) where {B<:Basis}
    OpenQuantumSystems.mul!(dpsi, psi, H, eltype(psi)(im), zero(eltype(psi)))
    return dpsi
end


function dschroedinger_dynamic(
    t,
    psi0::T,
    f::Function,
    dpsi::T,
) where {T<:Union{StateVector,AbstractOperator}}
    H = f(t, psi0)
    dschroedinger(psi0, H, dpsi)
end


function check_schroedinger(psi::Ket, H::AbstractOperator)
    check_multiplicable(H, psi)
    check_samebases(H)
end

function check_schroedinger(psi::Bra, H::AbstractOperator)
    check_multiplicable(psi, H)
    check_samebases(H)
end
