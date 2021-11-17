using OpenQuantumSystems
using LinearAlgebra
using SparseArrays


mode1 = Mode(0.2, 1.0)
Energy = [0.0, 200.0]
mol1 = Molecule([mode1], 2, [2.0, 200.0])
mol2 = Molecule([mode1], 2, [3.0, 300.0])
agg = Aggregate([mol1, mol2])
aggInds = getIndices(agg)
FCFact = getFranckCondonFactors(agg, aggInds)

agg.coupling[2, 3] = 50
agg.coupling[3, 2] = 50

Ham_S = getAggHamiltonianSystemSmall(agg; groundEnergy=true)
print(size(Ham_S.data))
print(Ham_S)

Ham_S = getAggHamiltonianSystemBig(agg, aggInds, FCFact; groundEnergy=true)
print(size(Ham_S.data))
print(Ham_S)

Ham_B = getAggHamiltonianBathSmall(agg; groundEnergy=true)
print(size(Ham_B.data))
print(Ham_B)

Ham_B = getAggHamiltonianBathBig(agg; groundEnergy=true)
print(size(Ham_B.data))
print(Ham_B)

Ham_0 = getAggHamSystemBath(agg, aggInds, FCFact; groundEnergy=true)
print(size(Ham_0.data))
print(Ham_0)