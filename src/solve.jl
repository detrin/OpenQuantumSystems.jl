
const _SOLVE_METHODS = (
    :lvn_si, :lvn_ss, :lvn_SI, :lvn_SS,
    :qme_exact_si, :qme_exact_ss, :qme_exact_SI, :qme_exact_SS,
    :qme_redfield, :qme_ansatz, :qme_iterative,
)

"""
    solve(W0, tspan, agg; method=:qme_exact_si, kwargs...)

Unified entry point for all solvers in OpenQuantumSystems.

Dispatches to the appropriate solver function based on `method`:

| `method`           | Solver function       | Notes                                     |
|:-------------------|:----------------------|:------------------------------------------|
| `:lvn_si`          | `LvN_sI`             | Liouville–von Neumann, small, interaction |
| `:lvn_ss`          | `LvN_sS`             | Liouville–von Neumann, small, Schrödinger |
| `:lvn_SI`          | `LvN_SI`             | Liouville–von Neumann, big, interaction   |
| `:lvn_SS`          | `LvN_SS`             | Liouville–von Neumann, big, Schrödinger   |
| `:qme_exact_si`    | `QME_sI_exact`        | QME exact, small, interaction             |
| `:qme_exact_ss`    | `QME_sS_exact`        | QME exact, small, Schrödinger             |
| `:qme_exact_SI`    | `QME_SI_exact`        | QME exact, big, interaction               |
| `:qme_exact_SS`    | `QME_SS_exact`        | QME exact, big, Schrödinger               |
| `:qme_redfield`    | `QME_sI_Redfield`     | Redfield QME                              |
| `:qme_ansatz`      | `QME_sI_ansatz`       | Ansatz QME (pass `ansatz=:test` etc.)     |
| `:qme_iterative`   | `QME_sI_iterative`    | Iterative QME (requires extra positional args via `iterative_args`) |

All keyword arguments except `method` and `iterative_args` are forwarded to the
underlying solver.

For `:qme_iterative`, pass the required positional arguments `(rho_0_int_t, W_0_bath_t)`
as a tuple via the `iterative_args` keyword:

```julia
solve(W0, tspan, agg; method=:qme_iterative,
      iterative_args=(rho_0_int_t, W_0_bath_t))
```
"""
function solve(
    W0,
    tspan,
    agg::Aggregate;
    method::Symbol = :qme_exact_si,
    iterative_args::Union{Tuple,Nothing} = nothing,
    kwargs...,
)
    method == :lvn_si       && return LvN_sI(W0, tspan, agg; kwargs...)
    method == :lvn_ss       && return LvN_sS(W0, tspan, agg; kwargs...)
    method == :lvn_SI       && return LvN_SI(W0, tspan, agg; kwargs...)
    method == :lvn_SS       && return LvN_SS(W0, tspan, agg; kwargs...)
    method == :qme_exact_si && return QME_sI_exact(W0, tspan, agg; kwargs...)
    method == :qme_exact_ss && return QME_sS_exact(W0, tspan, agg; kwargs...)
    method == :qme_exact_SI && return QME_SI_exact(W0, tspan, agg; kwargs...)
    method == :qme_exact_SS && return QME_SS_exact(W0, tspan, agg; kwargs...)
    method == :qme_redfield && return QME_sI_Redfield(W0, tspan, agg; kwargs...)
    method == :qme_ansatz   && return QME_sI_ansatz(W0, tspan, agg; kwargs...)
    if method == :qme_iterative
        iterative_args === nothing &&
            throw(ArgumentError("method=:qme_iterative requires iterative_args=(rho_0_int_t, W_0_bath_t)"))
        return QME_sI_iterative(W0, iterative_args[1], iterative_args[2], tspan, agg; kwargs...)
    end
    throw(ArgumentError(
        "Unknown method :$method. Supported methods: $(join(_SOLVE_METHODS, ", :"))"
    ))
end
