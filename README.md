# OpenQuantumSystems.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://detrin.github.io/OpenQuantumSystems.jl/dev/)
[![CI](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/detrin/OpenQuantumSystems.jl)

**OpenQuantumSystems.jl** is a numerical framework written in Julia that makes it easy to simulate various kinds of open quantum systems with main focus on quantum biology in finite basis. It is inspired by the
[quantarhei](https://github.com/tmancal74/quantarhei) and
[QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

## Installation

You can obtain OpenQuantumSystems using Julia's Pkg REPL-mode (hitting `]` as the first character of the command prompt):

```julia
(v1.6) pkg> add OpenQuantumSystems
```

or with `using Pkg; Pkg.add("OpenQuantumSystems")`.

### Questions & Contributions

The package is still under development. If you have any questions or need help,
see [gitter channel](https://gitter.im/OpenQuantumSystems-jl/community) and ask
away. Also, contributions of any kind are always welcome! If you have any ideas
for bug fixes, new features, optimisation or unit tests suggest it straight away
or create a pull request.

### Motivation
This package is a backbone for calculations for my master thesis supervised by [doc. Tomáš Mančal](http://www.mancal.cz/).

### Roadmap

What is OpenQuantumSystems.jl capable of right now (v0.5.0)?
- Supports systems with a finite basis.
- Construct Hamiltonians for aggregates of molecules. Bath degrees of freedom are modelled as LHOs with shifted potentials for excited local electronic states.
- Trace over bath degrees of freedom using Franck-Condon factors and reduced density matrix.
- Dynamics solvers:
  - Exact dynamics (U(t) operator)
  - Schrödinger equation
  - Liouville-von Neumann equation
  - Quantum Master Equation — ansatz, iterative, and Redfield variants
- Unified `solve()` entry point with `SimulationResult` return type.
- Prepare initial density matrix from a laser-excited pulse.
- Evaluate the bath part of the density matrix using various ansatzes (closed-form QME with RDM).
- Iterative correction of the bath part of the reduced density matrix.
- Memory kernel as superoperator.
- Rate constants from dynamics and memory kernel.
- Physical validation of simulation results (trace, positivity).
- Scoring and comparison of reduced density matrices over time.
- Convenience constructors for common systems (dimer, trimer, linear chain).
- Solver selection guide and API glossary in documentation.

What would be nice to have?
- Loading Hamiltonian and data storage for aggregate.
- Anharmonic oscillators.
- Double excited states.
- Interface with quantarhei, or creating an interface in quantarhei package.

### Branches status

Master

[![CI](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml)
[![CI-nightly-julia](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-nightly-julia.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-nightly-julia.yml)
[![CI-short](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml)
[![Deploy Nightly](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/deploy-nightly.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/deploy-nightly.yml)

Develop

[![CI-short](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml/badge.svg?branch=devel)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml)
[![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/devel)

### Testing
You can use `make.bat` or `make` file for testing, but it will install all the dependencies before executing tests. I advise using jupyter notebook `src/local_test.ipynb` with julia kernel for tests. Package `Revise.jl` will load all the functions and by using `include("src/test_file.jl")` you can quickly debug and add new tests.

### Developing
Please use formatting in `make.bat` or `make` file. I suggest using `Revise.jl` with jupyter notebook in `src/local.ipynb` for adding new features. It will spare you of the painful precompilation time.

### Credit

[quantarhei](https://github.com/tmancal74/quantarhei) - Concepts of aggregate construction, Hamiltonian of aggregate construction and trace over bath degrees of freedom were implemented in `quantarhei` in `python`.

[QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl) - Provides building elements for this package, such as the implementation of operators and superoperators (and many more).