
# include("core.jl")

abstract type AbstractAggregate end

mutable struct Aggregate <: AbstractAggregate
    core::AggregateCore
    tools::AggregateTools
    operators::AggregateOperators
end

AggregateTools(agg::Aggregate)::AggregateTools = AggregateTools(agg.core)

AggregateOperators(agg::Aggregate; groundEnergy::Bool = true, vib_basis::Symbol = :ground_ground)::AggregateOperators =
    AggregateOperators(agg.core, agg.tools; groundEnergy=groundEnergy, vib_basis=vib_basis)

"""
    setupAggregate(agg; groundEnergy=true, verbose=false)

Generate all basic data from the [`Aggregate`](@ref). Returns
`aggInds, vibindices, aggIndLen, basis, FCFact, FCProd, Ham, Ham_0, Ham_I`.

"""
function setupAggregate(aggCore::AggregateCore; groundEnergy::Bool = true, vib_basis::Symbol = :ground_ground)::Aggregate
    aggTools = AggregateTools(aggCore)
    aggOperators = AggregateOperators(aggCore, aggTools; groundEnergy = groundEnergy, vib_basis = vib_basis)
    return Aggregate(aggCore, aggTools, aggOperators)
end

function setupAggregate!(agg::Aggregate; groundEnergy::Bool = true, vib_basis::Symbol = :ground_ground)::Aggregate
    agg.tools = AggregateTools(agg.core)
    agg.operators = AggregateOperators(agg.core, agg.tools; groundEnergy = groundEnergy, vib_basis = vib_basis)
    return agg
end

Base.:(==)(x::Aggregate, y::Aggregate) =
    x.core == y.core && x.tools == y.tools && x.operators == y.operators
