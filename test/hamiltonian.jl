
using OpenQuantumSystems
using DelimitedFiles

Nvib = 3

HR = 0.01
shift = (2.0 * HR)^0.5
omega = 300.
mols = [
    Molecule([Mode(omega, shift)], Nvib, [12500., 12750.]),
    Molecule([Mode(omega, shift)], Nvib, [12500., 12800.]),
]

agg = Aggregate(mols)
for mol_i in 2:length(agg.molecules)
    agg.coupling[mol_i, mol_i+1] = 100
    agg.coupling[mol_i+1, mol_i] = 100
end

aggInds = getIndices(agg; groundState=true)
aggIndLen = length(aggInds)
base = GenericBasis([aggIndLen])
FCFact = getFranckCondonFactors(agg, aggInds; groundState=true)
Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState=true, groundEnergy=true)
println(size(Ham.data))
writedlm("test/hamiltonian.csv",  Ham.data, ',')