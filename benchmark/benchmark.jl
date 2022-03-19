
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

import OrdinaryDiffEq, DiffEqCallbacks, DelayDiffEq
import SparseArrays: sparse
import QuantumOpticsBase

BenchmarkTools.DEFAULT_PARAMETERS.samples = 100

# dic_commits = JSON.parsefile("test/benchmark_commits.json")
open("benchmark/benchmark_commits.json", "r") do f
    global dic_commits = JSON.parse(f)
end

commit_hash = readlines(`git log -n 1 --no-merges --format=%h`)[1]
commit_date = readlines(`git log -n 1 --no-merges --format=%cs`)[1]
dic_commit = Dict()
dic_commit["commit_date"] = commit_date

println

## Aggregate building
### agg_dimer_small
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5.0
t1_1 = @benchmark begin
    HR = 0.01
    shift = (2.0 * HR)
    modes = [Mode(180.0, shift)]
    mols = [
        Molecule([Mode(300.0, shift)], 3, [12500.0, 12710.0]),
        Molecule([Mode(300.0, shift)], 3, [12500.0, 12690.0]),
    ]

    agg = Aggregate(mols)
    for mol_i = 2:length(agg.molecules)
        agg.coupling[mol_i, mol_i+1] = 200
        agg.coupling[mol_i+1, mol_i] = 200
    end

    aggIndices = getIndices(agg)
    vibindices = getVibIndices(agg, aggIndices)
    aggIndLen = length(aggIndices)
    base = GenericBasis([aggIndLen])
    FCFact = getFranckCondonFactors(agg, aggIndices)
    FCProd = getFCProd(agg, FCFact, aggIndices, vibindices)
    Ham = getAggHamiltonian(agg, aggIndices, FCFact)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundEnergy = false)
    b = GenericBasis([aggIndLen])
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact)
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
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 20.0
t1_2 = @benchmark begin
    HR = 0.01
    shift = (2.0 * HR)
    modes = [Mode(180.0, shift)]
    mols = [
        Molecule([Mode(300.0, shift)], 7, [12500.0, 12710.0]),
        Molecule([Mode(300.0, shift)], 7, [12500.0, 12690.0]),
    ]

    agg = Aggregate(mols)
    for mol_i = 2:length(agg.molecules)
        agg.coupling[mol_i, mol_i+1] = 200
        agg.coupling[mol_i+1, mol_i] = 200
    end

    aggIndices = getIndices(agg)
    vibindices = getVibIndices(agg, aggIndices)
    aggIndLen = length(aggIndices)
    base = GenericBasis([aggIndLen])
    FCFact = getFranckCondonFactors(agg, aggIndices)
    FCProd = getFCProd(agg, FCFact, aggIndices, vibindices)
    Ham = getAggHamiltonian(agg, aggIndices, FCFact)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg; groundEnergy = false)
    b = GenericBasis([aggIndLen])
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact)
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

### agg_building
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
t1_4 = @benchmark begin
    HR = 0.01
    shift = (2.0 * HR)
    modes = [Mode(180.0, shift)]
    mols = [
        Molecule([Mode(300.0, shift)], 7, [12500.0, 12710.0]),
        Molecule([Mode(300.0, shift)], 7, [12500.0, 12690.0]),
    ]

    agg = Aggregate(mols)
end
HR = 0.01
shift = (2.0 * HR)
modes = [Mode(180.0, shift)]
mols = [
    Molecule([Mode(300.0, shift)], 7, [12500.0, 12710.0]),
    Molecule([Mode(300.0, shift)], 7, [12500.0, 12690.0]),
]

agg = Aggregate(mols)
for mol_i = 2:length(agg.molecules)
    agg.coupling[mol_i, mol_i+1] = 200
    agg.coupling[mol_i+1, mol_i] = 200
end

aggIndices = getIndices(agg)
vibindices = getVibIndices(agg, aggIndices)
aggIndLen = length(aggIndices)
base = GenericBasis([aggIndLen])

### getFranckCondonFactors
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
t1_5 = @benchmark begin
    FCFact = getFranckCondonFactors(agg, aggIndices)
end
FCFact = getFranckCondonFactors(agg, aggIndices)

### getFCProd
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100.0
t1_6 = @benchmark begin
    FCProd = getFCProd(agg, FCFact, aggIndices, vibindices)
end
FCProd = getFCProd(agg, FCFact, aggIndices, vibindices)

### getFCProd
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
t1_7 = @benchmark begin
    Ham = getAggHamiltonian(agg, aggIndices, FCFact)
end
Ham = getAggHamiltonian(agg, aggIndices, FCFact)
Ham_bath = getAggHamiltonianBath(agg)
Ham_sys = getAggHamiltonianSystem(agg; groundEnergy = false)
b = GenericBasis([aggIndLen])
b_sys = GenericBasis([size(Ham_sys, 1)])
b_bath = GenericBasis([size(Ham_bath, 1)])

### getAggHamiltonianInteraction
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
t1_8 = @benchmark begin
    Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact)
end

dic_commit["agg_dimer_small"] = t1_1.times # 27
dic_commit["agg_dimer_medium"] = t1_2.times # 147
dic_commit["agg_building"] = t1_4.times
dic_commit["getFranckCondonFactors"] = t1_5.times
dic_commit["getFCProd"] = t1_6.times
dic_commit["getAggHamiltonian"] = t1_7.times
dic_commit["getAggHamiltonianInteraction"] = t1_8.times


## Simulations
HR = 0.01
shift = (2.0 * HR)
modes = [Mode(180.0, shift)]
mols = [
    Molecule([Mode(300.0, shift)], 3, [12500.0, 12710.0]),
    Molecule([Mode(300.0, shift)], 3, [12500.0, 12690.0]),
]

agg = Aggregate(mols)
for mol_i = 2:length(agg.molecules)
    agg.coupling[mol_i, mol_i+1] = 200
    agg.coupling[mol_i+1, mol_i] = 200
end

aggIndices = getIndices(agg)
vibindices = getVibIndices(agg, aggIndices)
aggIndLen = length(aggIndices)
base = GenericBasis([aggIndLen])
FCFact = getFranckCondonFactors(agg, aggIndices)
FCProd = getFCProd(agg, FCFact, aggIndices, vibindices)
Ham = getAggHamiltonian(agg, aggIndices, FCFact)
basis = GenericBasis([length(aggIndices)])

Ham_bath = getAggHamiltonianBath(agg)
Ham_sys = getAggHamiltonianSystem(agg; groundEnergy = false)
b = GenericBasis([aggIndLen])
b_sys = GenericBasis([size(Ham_sys, 1)])
b_bath = GenericBasis([size(Ham_bath, 1)])

Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact)
Ham_S = Ham - Ham_int
Ham_B = tensor(Ham_bath, OneDenseOperator(b_sys))
E0 = Ham_B.data[1, 1]
for I = 1:aggIndLen
    Ham_B.data[I, I] -= E0
end
Ham_B = DenseOperator(b, b, Ham_B.data)

H_lambda, H_S = eigen(Ham_S.data)
H_Sinv = inv(H_S)
H_lambda = diagm(H_lambda)

t_max = 0.01
t_count = 200
t0 = 0.0
t_step = (t_max - t0) / (t_count)
tspan = [t0:t_step:t_max;]
T = 300
mu_array = [[2, 1]]
W0_1 = thermal_state(T, mu_array, Ham, vibindices, aggIndices; diagonalize = false)
mu_array = [[1, 2]]
W0_2 = thermal_state(T, mu_array, Ham, vibindices, aggIndices; diagonalize = false)
W0 = 0.8 * W0_1 + 0.2 * W0_2
W0 = DenseOperator(W0.basis_l, W0.basis_r, complex(W0.data))
W0_bath = get_rho_bath(W0, agg, FCProd, aggIndices, vibindices)
W0_bath = DenseOperator(W0_bath.basis_l, W0_bath.basis_r, complex(W0_bath.data))
rho0 = trace_bath(W0, agg, FCProd, aggIndices, vibindices)
rho0 = DenseOperator(rho0.basis_l, rho0.basis_r, complex(rho0.data))

### trace_bath
t2_0 = @benchmark begin
    rho0 = trace_bath(W0, agg, FCProd, aggIndices, vibindices)
end

### evolutionExact
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
t2_1 = @benchmark begin
    op_array = evolutionExact(W0, tspan, Ham)
end

### evolutionApproximate
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2.0
t2_2 = @benchmark begin
    op_array = evolutionApproximate(W0, tspan, Ham)
end

### schroedinger
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
psi0 = randstate(basis)
t2_3 = @benchmark begin
    Tspan, psi_t = schroedinger(
        psi0,
        tspan,
        Ham;
        reltol = 1e-4,
        abstol = 1e-4,
        alg = OrdinaryDiffEq.Tsit5(),
    )
end

### schroedinger
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
t2_4 = @benchmark begin
    Tspan, rho_t = liouvilleVonNeumann(
        W0,
        tspan,
        Ham;
        reltol = 1e-4,
        abstol = 1e-4,
        alg = OrdinaryDiffEq.Tsit5(),
    )
end

### evolutionOperatorIterator
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10.0
elLen = length(agg.molecules)
rho_t_exact = zeros(ComplexF64, length(tspan), elLen + 1, elLen + 1)
t2_5 = @benchmark begin
    t_i = 0
    foreach(
        evolutionOperatorIterator(Ham, tspan; diagonalize = false, approximate = true),
    ) do U_op
        t_i = t_i + 1
        # println(t_i / 150. * 100)
        W = U_op * W0 * U_op'
        rho_traced = trace_bath(W.data, agg, FCProd, aggIndices, vibindices)
        rho_t_exact[t_i, :, :] = rho_traced
    end
end

### master_int
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 20.0
p = (
    Ham_S,
    Ham_int,
    H_lambda,
    H_S,
    H_Sinv,
    W0,
    W0_bath,
    agg,
    FCProd,
    aggIndices,
    vibindices,
    ComplexF64,
)
t2_6 = @benchmark begin
    Tspan, W_t_master_int = master_int(
        W0,
        tspan,
        Ham_S,
        Ham_int;
        reltol = 1.0e-4,
        abstol = 1.0e-4,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5()),
    )

    elLen = length(agg.molecules)
    rho_t_master = zeros(ComplexF64, length(tspan), elLen + 1, elLen + 1)
    for t_i = 1:length(tspan)
        t = tspan[t_i]
        U_op_S = evolutionOperator(Ham_sys, t)
        rho = trace_bath(W_t_master_int[t_i], agg, FCProd, aggIndices, vibindices)
        rho = U_op_S * rho * U_op_S'
        rho_t_master[t_i, :, :] = rho.data
    end
end

### evolutionOperatorIterator
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 20.0
p = (
    Ham_S,
    Ham_int,
    H_lambda,
    H_S,
    H_Sinv,
    Ham_B,
    W0,
    W0_bath,
    agg,
    FCProd,
    aggIndices,
    vibindices,
    ComplexF64,
)
t2_7 = @benchmark begin
    Tspan, rho_t_ansatz_int = master_ansatz(
        rho0,
        tspan,
        p;
        reltol = 1.0e-4,
        abstol = 1.0e-4,
        alg = DelayDiffEq.MethodOfSteps(DelayDiffEq.Tsit5()),
    )

    elLen = length(agg.molecules)
    rho_t_ansatz = zeros(ComplexF64, length(tspan), elLen + 1, elLen + 1)
    for t_i = 1:length(tspan)
        t = tspan[t_i]
        U_op_S = evolutionOperator(Ham_sys, t)
        rho = U_op_S * rho_t_ansatz_int[t_i] * U_op_S'
        rho_t_ansatz[t_i, :, :] = rho.data
    end
end

dic_commit["trace_bath"] = t2_0.times
dic_commit["evolutionExact"] = t2_1.times
dic_commit["evolutionApproximate"] = t2_2.times
dic_commit["schroedinger"] = t2_3.times
dic_commit["liouvilleVonNeumann"] = t2_4.times
dic_commit["evolutionOperatorIterator"] = t2_5.times
dic_commit["master_int"] = t2_6.times
dic_commit["master_ansatz"] = t2_7.times

dic_commits[commit_hash] = dic_commit
json_string = JSON.json(dic_commits)

open("benchmark/benchmark_commits.json", "w") do f
    write(f, json_string)
end
