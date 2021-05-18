# push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterTools
using OpenQuantumSystems

makedocs(
    format=     Documenter.HTML(
                    prettyurls= !("local" in ARGS),
                    canonical=  "https://detrin.github.io/OpenQuantumSystems.jl//latest/",
                    analytics=  "G-F12F2Y92VY",
                    ),
    clean=      false,
    sitename=   "OpemQuantumSystems.jl",
    pages=      [
                "Home"      =>  "index.md",
                ]
)

deploydocs(

    repo=       "github.com/detrin/OpenQuantumSystems.jl.git",
    target=     "build",
    deps=       nothing,
    make=       nothing,
)
#     modules = [OpemQuantumSystems],