
using BenchmarkTools
using Plots
using ProgressBars
using OpenQuantumSystems
# using OpenQuantumSystemsPrivate
using LinearAlgebra
using Random
using QuadGK
using ProgressMeter
using Expokit
# using DifferentialEquations
using DataFrames
using JSON
using Git

Random.seed!(0)

import OrdinaryDiffEq, DiffEqCallbacks, DelayDiffEq
import SparseArrays: sparse
import QuantumOpticsBase

# dic_commits = JSON.parsefile("test/benchmark_commits.json")
open("benchmark/benchmark_commits.json","r") do f 
    global dic_commits = JSON.parse(f)
end

commit_hash = readlines(`git log -n 1 --no-merges --format=%h`)[1]
commit_date = readlines(`git log -n 1 --no-merges --format=%cs`)[1]
dic_commit = Dict()
dic_commit["commit_date"] = commit_date

### agg_dimer_small
t1_1 = @benchmark begin
    HR = 0.01
    shift = (2.0 * HR)
    modes = [Mode(180., shift)]
    mols = [
        Molecule([Mode(300., shift)], 3, [12500., 12710.]),
        Molecule([Mode(300., shift)], 3, [12500., 12690.])
    ]

    agg = Aggregate(mols)
    for mol_i in 2:length(agg.molecules)
        agg.coupling[mol_i, mol_i+1] = 200
        agg.coupling[mol_i+1, mol_i] = 200
    end

    aggIndices = getIndices(agg; groundState=true)
    vibindices = getVibIndices(agg, aggIndices)
    aggIndLen = length(aggIndices)
    base = GenericBasis([aggIndLen])
    FCFact = getFranckCondonFactors(agg, aggIndices; groundState=true)
    FCProd = getFCProd(agg, FCFact, aggIndices, vibindices; groundState = true)
    Ham = getAggHamiltonian(agg, aggIndices, FCFact; groundState=true)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundState=true, groundEnergy=false)
    b = GenericBasis([aggIndLen])
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact; groundState=true)
    Ham_S = Ham - Ham_int
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])
    Ham_B = tensor(Ham_bath, OneDenseOperator(b_sys))
    E0 = Ham_B.data[1, 1]
    for I = 1:aggIndLen
        Ham_B.data[I, I] -= E0
    end
    Ham_B = DenseOperator(b, b, Ham_B.data)
end

### agg_dimer_medium
t1_2 = @benchmark begin
    HR = 0.01
    shift = (2.0 * HR)
    modes = [Mode(180., shift)]
    mols = [
        Molecule([Mode(300., shift)], 7, [12500., 12710.]),
        Molecule([Mode(300., shift)], 7, [12500., 12690.])
    ]

    agg = Aggregate(mols)
    for mol_i in 2:length(agg.molecules)
        agg.coupling[mol_i, mol_i+1] = 200
        agg.coupling[mol_i+1, mol_i] = 200
    end

    aggIndices = getIndices(agg; groundState=true)
    vibindices = getVibIndices(agg, aggIndices)
    aggIndLen = length(aggIndices)
    base = GenericBasis([aggIndLen])
    FCFact = getFranckCondonFactors(agg, aggIndices; groundState=true)
    FCProd = getFCProd(agg, FCFact, aggIndices, vibindices; groundState = true)
    Ham = getAggHamiltonian(agg, aggIndices, FCFact; groundState=true)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundState=true, groundEnergy=false)
    b = GenericBasis([aggIndLen])
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact; groundState=true)
    Ham_S = Ham - Ham_int
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])
    Ham_B = tensor(Ham_bath, OneDenseOperator(b_sys))
    E0 = Ham_B.data[1, 1]
    for I = 1:aggIndLen
        Ham_B.data[I, I] -= E0
    end
    Ham_B = DenseOperator(b, b, Ham_B.data)
end

### agg_dimer_big
t1_3 = @benchmark begin
    HR = 0.01
    shift = (2.0 * HR)
    modes = [Mode(180., shift)]
    mols = [
        Molecule([Mode(300., shift)], 6, [12500., 12710.]),
        Molecule([Mode(300., shift)], 6, [12500., 12690.]),
        Molecule([Mode(300., shift)], 6, [12500., 12690.])
    ]

    agg = Aggregate(mols)
    for mol_i in 2:length(agg.molecules)
        agg.coupling[mol_i, mol_i+1] = 200
        agg.coupling[mol_i+1, mol_i] = 200
    end

    aggIndices = getIndices(agg; groundState=true)
    vibindices = getVibIndices(agg, aggIndices)
    aggIndLen = length(aggIndices)
    base = GenericBasis([aggIndLen])
    FCFact = getFranckCondonFactors(agg, aggIndices; groundState=true)
    FCProd = getFCProd(agg, FCFact, aggIndices, vibindices; groundState = true)
    Ham = getAggHamiltonian(agg, aggIndices, FCFact; groundState=true)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundState=true, groundEnergy=false)
    b = GenericBasis([aggIndLen])
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact; groundState=true)
    Ham_S = Ham - Ham_int
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])
    Ham_B = tensor(Ham_bath, OneDenseOperator(b_sys))
    E0 = Ham_B.data[1, 1]
    for I = 1:aggIndLen
        Ham_B.data[I, I] -= E0
    end
    Ham_B = DenseOperator(b, b, Ham_B.data)
end

dic_commit["agg_dimer_small"] = t1_1.times # 27
dic_commit["agg_dimer_medium"] = t1_2.times # 147
dic_commit["agg_dimer_big"] = t1_3.times # 864

dic_commits[commit_hash] = dic_commit
json_string = JSON.json(dic_commits)

open("benchmark/benchmark_commits.json","w") do f 
    write(f, json_string) 
end