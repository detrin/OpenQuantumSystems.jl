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

What is OpenQuantumSystems.jl capable of right now?

- It only supports systems with a finite basis.
- Calculate Hamiltonian for Aggregates of molecules. Molecules have system Hamiltonian that consists only of relevant degrees of freedom (electronic basis) and bath degrees of freedom. Bath consists of LHOs with shifted potentials for excited local electronic states.
- Trace over bath degrees of freedom using Franck-Condon Factors and reduced density matrix.
- It can calculate dynamics for aggregate named as follows
  - Exact dynamics (calculate U(t) operator)
  - Schrodinger dynamics
  - Liuville-von Neumann dynamics
  - Quantum Master Equation dynamics
- OQS can also prepare the initial condition for density matrix as a laser-excited pulse.
- It can evaluate the remaining part of the bath. The majority of my master thesis proposes various kinds of ansatzes that aim to find a closed-form of QME with RDM.
- Scoring of obtained RDM in time in comparison to another RDM.
- Calculate memory kernel as superoperator.

What can we expect in the near future?

- Iterative approach for bath part of the reduced density matrix.
- Calculation of rate constants from dynamics and memory kernel.
- Better documentation with examples of how to use this package.
- QME with iterative correction of bath part.
- Calculation of dipoles for molecules given coordinates.
- Implementation of Foerster and modified Redfield theories of excitation energy transfer for finite systems (and perhaps special infinite cases).

What would be nice to have?

- Decomposition of mixed states into a linear combination of pure states using linear programming.
- Loading Hamiltonian and data storage for aggregate.
- GPU support for Schrodinger equation and possibly state decomposition.
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
