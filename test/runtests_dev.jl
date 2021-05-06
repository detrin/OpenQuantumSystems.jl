
using Revise
using Test

wait_for_key(prompt) = (println(stdout, prompt); read(stdin, 1); nothing)

names = [
    "test_operators_dense.jl",
    "test_superoperators.jl",
    "test_metrics.jl",
    "test_molecules.jl",
    "test_aggregate.jl",
    "test_evolution.jl",
    "test_schroedinger.jl",
    "test_liouville.jl",
    "test_master.jl",
]

names = ["test_master.jl"]

while true 
    for name in names
        if startswith(name, "test_") && endswith(name, ".jl")
            try
                include(name)
            catch
                nothing
            end
        end
    end

    wait_for_key("press any key to run tests")
end