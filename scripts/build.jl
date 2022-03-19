import SimpleContainerGenerator

pkgs = [
    "DelayDiffEq", 
    "DiffEqCallbacks",
    "LinearAlgebra",
    "OrdinaryDiffEq",
    "QuadGK",
    "QuantumOpticsBase",
    "Random",
    "Random",
    "Reexport",
    "ResumableFunctions",
    "SparseArrays",
    "StableRNGs",
]

julia_version = v"1.7.0"

SimpleContainerGenerator.create_dockerfile(
    pkgs;
    julia_version = julia_version,
    output_directory = pwd()
)