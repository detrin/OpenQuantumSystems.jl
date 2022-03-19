if length(ARGS) != 0 && ARGS[1] == "local"
    include("local.jl"); exit()
end

names = [
    "test_operators_dense.jl",
    "test_superoperators.jl",
    "test_metrics.jl",
    "test_molecules.jl",
    "test_aggregateCore.jl",
    "test_aggregate.jl",
    "test_evolution.jl",
    "test_schroedinger.jl",
    "test_liouville.jl",
    "test_interaction_picture.jl",
    "test_master_exact.jl",
    "test_trace.jl",
    "test_initial_state.jl",
    "test_memory_kernel.jl",
    "test_master_ansatz.jl",
    "test_postprocessing.jl",
    "test_scoring.jl",
]

#=
detected_tests =
    filter(name -> startswith(name, "test_") && endswith(name, ".jl"), readdir("."))

unused_tests = setdiff(detected_tests, names)
if length(unused_tests) != 0
    error("The following tests are not used:\n", join(unused_tests, "\n"))
end

unavailable_tests = setdiff(names, detected_tests)
if length(unavailable_tests) != 0
    error("The following tests could not be found:\n", join(unavailable_tests, "\n"))
end
=#

if length(ARGS) != 0
    names = map(x -> "test_" * x * ".jl", ARGS)
end
# names = ["test_trace.jl"]
for name in names
    if startswith(name, "test_") && endswith(name, ".jl")
        include(name)
    end
end
