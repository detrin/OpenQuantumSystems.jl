
using Revise
using Pkg; Pkg.activate("../"); Pkg.build()
# using Pkg; Pkg.activate(); Pkg.build()
using OpenQuantumSystems
Pkg.status("OpenQuantumSystems")
using LinearAlgebra
using SparseArrays

# using Plots

mode1 = Mode(0.2, 1.0)
mode2 = Mode(0.2, 1.0)
Energy = [0.0, 200.0]
mol1 = Molecule([mode1, mode2], 3, [2.0, 200.0])
mol2 = Molecule([mode1], 3, [3.0, 300.0])
aggCore = AggregateCore([mol1, mol2])
agg = Aggregate(aggCore)
agg.core.coupling[2, 3] = 50
agg.core.coupling[3, 2] = 50
agg.tools = AggregateTools(agg)


# aggregateOperators = AggregateOperators(agg)
