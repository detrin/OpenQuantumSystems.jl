import OrdinaryDiffEq, DiffEqCallbacks, DelayDiffEq

function recast! end

"""
    integrate(tspan, df::Function, x0::Vector{ComplexF64},
            state::T, dstate::T, fout::Function; kwargs...)

Integrate using OrdinaryDiffEq
"""
function integrate(
    tspan,
    df::Function,
    x0::X,
    state::T,
    dstate::T,
    fout::Function;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.DP5(),
    steady_state = false,
    tol = 1e-3,
    save_everystep = false,
    saveat = tspan,
    callback = nothing,
    kwargs...,
) where {T,X}

    function df_(dx::T, x::T, p, t) where {T}
        recast!(x, state)
        recast!(dx, dstate)
        df(t, state, dstate)
        recast!(dstate, dx)
    end
    function fout_(x, t, integrator)
        recast!(x, state)
        fout(t, state)
    end

    out_type = pure_inference(fout, Tuple{eltype(tspan),typeof(state)})

    out = DiffEqCallbacks.SavedValues(eltype(tspan), out_type)

    scb = DiffEqCallbacks.SavingCallback(
        fout_,
        out,
        saveat = saveat,
        save_everystep = save_everystep,
        save_start = false,
    )

    prob = OrdinaryDiffEq.ODEProblem{true}(df_, x0, (tspan[1], tspan[end]))

    if steady_state
        affect! = function (integrator)
            !save_everystep && scb.affect!(integrator, true)
            OrdinaryDiffEq.terminate!(integrator)
        end
        _cb = OrdinaryDiffEq.DiscreteCallback(
            SteadyStateCondtion(copy(state), tol, state),
            affect!;
            save_positions = (false, false),
        )
        cb = OrdinaryDiffEq.CallbackSet(_cb, scb)
    else
        cb = scb
    end

    full_cb = OrdinaryDiffEq.CallbackSet(callback, cb)

    sol = OrdinaryDiffEq.solve(
        prob,
        alg;
        reltol = reltol,
        abstol = abstol,
        save_everystep = false,
        save_start = false,
        save_end = false,
        callback = full_cb,
        kwargs...,
    )
    out.t, out.saveval
end

function integrate(
    tspan,
    df::Function,
    x0::X,
    state::T,
    dstate::T,
    ::Nothing;
    kwargs...,
) where {T,X}
    function fout(t, state::T)
        copy(state)
    end
    integrate(tspan, df, x0, state, dstate, fout; kwargs...)
end


function integrate_delayed(
    tspan,
    df::Function,
    h::Function,
    x0::X,
    state::T,
    dstate::T,
    fout::Function;
    p = nothing,
    reltol::Float64 = 1.0e-6,
    abstol::Float64 = 1.0e-6,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    steady_state = false,
    tol = 1e-3,
    save_everystep = false,
    saveat = tspan,
    callback = nothing,
    kwargs...,
) where {T,X}

    function df_(dx::T, x::T, h, p, t) where {T}
        recast!(x, state)
        recast!(dx, dstate)
        df(t, state, dstate, h, p)
        recast!(dstate, dx)
    end
    function fout_(x, t, integrator)
        recast!(x, state)
        fout(t, state)
    end

    out_type = pure_inference(fout, Tuple{eltype(tspan),typeof(state)})

    out = DiffEqCallbacks.SavedValues(eltype(tspan), out_type)

    scb = DiffEqCallbacks.SavingCallback(
        fout_,
        out,
        saveat = saveat,
        save_everystep = save_everystep,
        save_start = false,
    )

    prob = DelayDiffEq.DDEProblem{true}(df_, x0, h, (tspan[1], tspan[end]))

    if steady_state
        affect! = function (integrator)
            !save_everystep && scb.affect!(integrator, true)
            DelayDiffEq.terminate!(integrator)
        end
        _cb = DelayDiffEq.DiscreteCallback(
            SteadyStateCondtion(copy(state), tol, state),
            affect!;
            save_positions = (false, false),
        )
        cb = DelayDiffEq.CallbackSet(_cb, scb)
    else
        cb = scb
    end

    full_cb = DelayDiffEq.CallbackSet(callback, cb)

    sol = DelayDiffEq.solve(
        prob,
        alg;
        reltol = reltol,
        abstol = abstol,
        save_everystep = false,
        save_start = false,
        save_end = false,
        callback = full_cb,
        kwargs...,
    )
    out.t, out.saveval
end

function integrate_delayed(
    tspan,
    df::Function,
    h::Function,
    x0::X,
    state::T,
    dstate::T,
    ::Nothing;
    kwargs...,
) where {T,X}
    function fout(t, state::T)
        copy(state)
    end
    integrate_delayed(tspan, df, h::Function, x0, state, dstate, fout; kwargs...)
end

struct SteadyStateCondtion{T,T2,T3}
    rho0::T
    tol::T2
    state::T3
end
function (c::SteadyStateCondtion)(rho, t, integrator)
    recast!(rho, c.state)
    dt = integrator.dt
    drho = tracedistance(c.rho0, c.state)
    c.rho0.data[:] = c.state.data
    drho / dt < c.tol
end


const QO_CHECKS = Ref(true)
"""
    @skiptimechecks

Macro to skip checks during time-dependent problems.
Useful for `master_dynamic` and similar functions.
"""
macro skiptimechecks(ex)
    return quote
        QO_CHECKS.x = false
        local val = $(esc(ex))
        QO_CHECKS.x = true
        val
    end
end

Base.@pure pure_inference(fout, T) = Core.Compiler.return_type(fout, T)

function recast!(x::Union{Vector,SubArray}, rho::Operator{B,B,T}) where {B<:Basis,T}
    rho.data = reshape(x, size(rho.data))
end
recast!(state::Operator{B,B}, x::SubArray) where {B<:Basis} = (x[:] = state.data)
recast!(state::Operator{B,B}, x::Vector) where {B<:Basis} = nothing

# Recasting needed for the ODE solver is just providing the underlying data
function recast!(x::T, rho::Operator{B,B,T}) where {B<:Basis,T}
    rho.data = x
end
recast!(rho::Operator{B,B,T}, x::T) where {B<:Basis,T} = nothing
