using OpenQuantumSystems
using LinearAlgebra
using SparseArrays
# using Plots

mode1 = Mode(0.2, 1.0)
Energy = [0.0, 200.0]
mol1 = Molecule([mode1], 3, [2.0, 200.0])
mol2 = Molecule([mode1], 3, [3.0, 300.0])
mol3 = Molecule([mode1], 3, [4.0, 400.0])
agg = Aggregate([mol1, mol2, mol3])
aggInds = getIndices(agg)
FCFact = getFranckCondonFactors(agg, aggInds)

agg.coupling[2, 3] = 50
agg.coupling[3, 2] = 50

println(aggInds)

Ham_I = getAggHamInteraction(agg)
print(size(Ham_I.data))
print(Ham_I)
# heatmap(Ham_I.data)
# show()

Ham = getAggHamiltonian(agg)
print(size(Ham.data))
print(Ham)