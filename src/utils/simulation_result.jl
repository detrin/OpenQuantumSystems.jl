"""
    SimulationResult

Wrapper around solver output that provides a uniform interface regardless of
which solver was used.

# Fields
- `tspan::Vector{Float64}` — time points at which states were recorded.
- `states::Vector{<:Operator}` — density matrices at each time point.
- `method::Symbol` — solver method that produced this result.
- `extra::Dict{Symbol,Any}` — optional extra data (e.g. `:W_1_bath_t`).
"""
struct SimulationResult{T<:Operator}
    tspan::Vector{Float64}
    states::Vector{T}
    method::Symbol
    extra::Dict{Symbol,Any}
end

function SimulationResult(tspan, states, method::Symbol)
    SimulationResult(tspan, states, method, Dict{Symbol,Any}())
end

"""
    populations(r::SimulationResult) -> Matrix{Float64}

Return a `length(r) x dim` real matrix whose `(t, k)` entry is the `k`-th
diagonal element of the density matrix at time step `t`.
"""
function populations(r::SimulationResult)
    n = length(r.states)
    dim = size(r.states[1].data, 1)
    pops = Matrix{Float64}(undef, n, dim)
    for i in 1:n
        for k in 1:dim
            pops[i, k] = real(r.states[i].data[k, k])
        end
    end
    pops
end

Base.getindex(r::SimulationResult, i) = r.states[i]
Base.length(r::SimulationResult) = length(r.states)
Base.iterate(r::SimulationResult) = iterate(r.states)
Base.iterate(r::SimulationResult, state) = iterate(r.states, state)
Base.eltype(::Type{SimulationResult{T}}) where {T} = T
Base.firstindex(r::SimulationResult) = firstindex(r.states)
Base.lastindex(r::SimulationResult) = lastindex(r.states)
