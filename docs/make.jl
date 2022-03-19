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
        "Documentation" => "documentation.md",
        "Tutorials" => ["Dimer" => "tutorials/dimer.md"],
    ],
)

deploydocs(repo = "github.com/detrin/OpenQuantumSystems.jl.git")
