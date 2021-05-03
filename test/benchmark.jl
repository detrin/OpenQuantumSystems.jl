# OpenQuantumSystems.jl has to be installed locally
# julia
# ]
# add <path-to-OpenQuantumSystems.jl>

using OpenQuantumSystems
using BenchmarkTools
using Random

Random.seed!(0)

println("# Benchmarks")
println("\nThis benchmark was automatically generated by test/benchmark.jl. In order to run this file you have to install BenchmarkTools.jl first.\n")
println("\nBenchmarked on: ", Sys.cpu_info()[1].model, "\n")
println("## Building")
println("### System 1\n")
modes = [Mode(2., 2.)]
mols = [
    Molecule(modes, 5, [0., 200.]),
    Molecule(modes, 5, [0., 300.])
    ]

agg1 = Aggregate(mols)
agg1.coupling[2, 3] = 100
agg1.coupling[3, 2] = 100
print("""```julia 
modes = [Mode(2., 2.)]
mols = [
    Molecule(modes, 5, [0., 200.]),
    Molecule(modes, 5, [0., 300.])
    ]

agg1 = Aggregate(mols)
agg1.coupling[2, 3] = 100
agg1.coupling[3, 2] = 100
```""")

aggIndices1 = getIndices(agg1; groundState=false)
basisLen = length(aggIndices1)
println("\nBasis size: ", basisLen,"\n")
basis1 = GenericBasis([length(aggIndices1)])
println("\nFranck Condon factors\n")
@time FCFact1 = getFranckCondonFactors(agg1, aggIndices1; groundState=false)
println("\nHamiltonian\n")
@time Ham1 = getAggHamiltonian(agg1, aggIndices1, FCFact1; groundState=false)

println("### System 2\n")
modes = [Mode(2., 2.)]
mols = [
    Molecule(modes, 5, [0., 200.]),
    Molecule(modes, 5, [0., 300.]),
    Molecule(modes, 5, [0., 300.])
    ]

agg2 = Aggregate(mols)
agg2.coupling[2, 3] = 100
agg2.coupling[3, 2] = 100
agg2.coupling[3, 4] = 100
agg2.coupling[4, 3] = 100
print("""```julia 
modes = [Mode(2., 2.)]
mols = [
    Molecule(modes, 5, [0., 200.]),
    Molecule(modes, 5, [0., 300.]),
    Molecule(modes, 5, [0., 300.])
    ]

agg2 = Aggregate(mols)
agg2.coupling[2, 3] = 100
agg2.coupling[3, 2] = 100
agg2.coupling[3, 4] = 100
agg2.coupling[4, 3] = 100
```""")

aggIndices2 = getIndices(agg2; groundState=false)
basisLen = length(aggIndices2)
println("\nBasis size: ", basisLen,"\n")
basis2 = GenericBasis([length(aggIndices2)])
println("\nFranck Condon factors\n")
@time FCFact2 = getFranckCondonFactors(agg2, aggIndices2; groundState=false)
println("\nHamiltonian\n")
@time Ham2 = getAggHamiltonian(agg2, aggIndices2, FCFact2; groundState=false)

println("### System 3\n")
modes = [Mode(2., 2.), Mode(3., 2.)]
mols = [
    Molecule(modes, 5, [0., 200.]),
    Molecule(modes, 5, [0., 300.])
    ]

agg3 = Aggregate(mols)
agg3.coupling[2, 3] = 100
agg3.coupling[3, 2] = 100
print("""```julia 
modes = [Mode(2., 2.), Mode(3., 2.)] 
mols = [ 
    Molecule(modes, 5, [0., 200.]), 
    Molecule(modes, 5, [0., 300.]) 
    ] 
    
agg3 = Aggregate(mols) 
agg3.coupling[2, 3] = 100 
agg3.coupling[3, 2] = 100
```""")

aggIndices3 = getIndices(agg3; groundState=false)
basisLen = length(aggIndices3)
println("\nBasis size: ", basisLen,"\n")
basis3 = GenericBasis([length(aggIndices3)])
println("\nFranck Condon factors\n")
@time FCFact3 = getFranckCondonFactors(agg3, aggIndices3; groundState=false)
println("\nHamiltonian\n")
@time Ham3 = getAggHamiltonian(agg3, aggIndices3, FCFact3; groundState=false)

if true
    println("## Evolution")
    println("\n `tspan = [0.0:0.02:1.0;]`\n")
    tspan = [0.0:0.05:1.0;]

    println("### System 1\n")
    psi0 = randstate(basis1)
    rho0 = dm(psi0)
    println("\nevolutionExact(rho0, Ham1, tspan)\n")
    @btime let
        op_array = evolutionExact(rho0, Ham1, tspan)
    end
    println("\nevolutionApproximate(rho0, Ham1, tspan)\n")
    @btime let
        op_array = evolutionApproximate(rho0, Ham1, tspan)
    end

    println("### System 2\n")
    psi0 = randstate(basis2)
    rho0 = dm(psi0)
    println("\nevolutionExact(rho0, Ham2, tspan)\n")
    @btime let
        op_array = evolutionExact(rho0, Ham2, tspan)
    end
    println("\nevolutionApproximate(rho0, Ham2, tspan)\n")
    @btime let
        op_array = evolutionApproximate(rho0, Ham2, tspan)
    end

    println("### System 3\n")
    psi0 = randstate(basis3)
    rho0 = dm(psi0)
    println("\nevolutionExact(rho0, Ham3, tspan)\n")
    @btime let
        op_array = evolutionExact(rho0, Ham3, tspan)
    end
    println("\nevolutionApproximate(rho0, Ham3, tspan)\n")
    @btime let
        op_array = evolutionApproximate(rho0, Ham3, tspan)
    end
end
if true
    println("## Schrödinger")
    tspan = [0.0:0.02:1.0;]
    println("\n `tspan = [0.0:0.02:1.0;]`\n")
    println("### System 1\n")
    psi0 = randstate(basis1)
    println("\nschroedinger(psi0, Ham1, tspan)\n")
    @btime let
        t, psi_t = schroedinger(psi0, Ham1, tspan)
    end

    println("### System 2\n")
    psi0 = randstate(basis2)
    println("\nschroedinger(psi0, Ham2, tspan)\n")
    @btime let
        t, psi_t = schroedinger(psi0, Ham2, tspan)
    end

    println("### System 3\n")
    psi0 = randstate(basis3)
    println("\nschroedinger(psi0, Ham3, tspan)\n")
    @btime let
        t, psi_t = schroedinger(psi0, Ham3, tspan)
    end
end
# tout, psi_t = schroedinger(psi0, Ham1, tspan)