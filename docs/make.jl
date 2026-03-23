# push!(LOAD_PATH,"../src/")

using Documenter, OpenQuantumSystems

makedocs(
    modules = [OpenQuantumSystems],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://detrin.github.io/OpenQuantumSystems.jl/latest/",
        analytics = "G-F12F2Y92VY",
    ),
    clean = false,
    sitename = "OpenQuantumSystems.jl",
    authors = "Daniel Herman",
    linkcheck = !("skiplinks" in ARGS),
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "Hamiltonian & Bath Model" => "theory/hamiltonian.md",
            "Quantum Master Equation" => "theory/qme.md",
            "Redfield Equations" => "theory/redfield.md",
            "Bath State Ansatze" => "theory/ansatz.md",
            "Iterative QME" => "theory/iterative.md",
            "Corrected Memory Kernel" => "theory/memory_kernel.md",
        ],
        "Solver Guide" => "solver_guide.md",
        "API Reference" => "documentation.md",
        "Naming Conventions & Glossary" => "glossary.md",
        "Tutorials" => [
            "Dimer" => "tutorials/dimer.md",
            "Förster Theory & Dipole Couplings" => "tutorials/forster.md",
            "Comparing Simulation Methods" => "tutorials/comparison.md",
        ],
    ],
)

deploydocs(repo = "github.com/detrin/OpenQuantumSystems.jl.git")
