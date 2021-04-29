
names = [
    "test_sortedindices.jl",
    "test_bases.jl",
    "test_states.jl",
    "test_operators.jl",
    "test_operators_dense.jl",
    "test_sparsematrix.jl",
    "test_operators_sparse.jl",
    "test_fock.jl",
    "test_spin.jl",
    "test_nlevel.jl",
    "test_state_definitions.jl",
    "test_metrics.jl",
    "test_embed.jl",
    "test_superoperators.jl",
    "test_pauli.jl",
]

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

for name in names
    if startswith(name, "test_") && endswith(name, ".jl")
        include(name)
    end
end
