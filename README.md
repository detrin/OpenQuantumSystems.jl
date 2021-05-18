# OpenQuantumSystems.jl

[![Join the chat at https://gitter.im/OpenQuantumSystems-jl/community](https://badges.gitter.im/OpenQuantumSystems-jl/community.svg)](https://gitter.im/OpenQuantumSystems-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://detrin.github.io/OpenQuantumSystems.jl/dev/)
[![CI](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/detrin/OpenQuantumSystems.jl)

**OpenQuantumSystems.jl** is a numerical framework written in [Julia] that makes it easy to simulate various kinds of open quantum systems with main focus on quantum biology. It is inspired by the
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

### Roadmap

1. - Efficient representation of operators given Hamiltonians in different
     Hilbert bases. ✓
   - Vibrational basis, spin basis, exciton basis. ✓
   - Calculation of Franck-Condon factors for multidimensional linear harmonic
     oscilators. ✓
2. - Schrodinger equation for solving dynamics. ✓
   - Solving dynamics for mixed states using exponentials, Liouville equation,
     QME. ✓
   - Aggregate with yml preferences.
   - Decomposition of mixed states into a linear combination of pure states
     using linear programming.
   - Loading Hamiltonian and data storage.
3. Implementing Foerster and modified Redfield theories of excitation energy
   transfer.
4. - GPU support for Schrodinger equation and possibly state decomposition.
   - Anharmonic oscillators.
   - Double excited states.

### Branches status

Master

[![CI](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml)
[![CI-nightly-julia](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-nightly-julia.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-nightly-julia.yml)
[![CI-short](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml)
[![Deploy Nightly](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/deploy-nightly.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/deploy-nightly.yml)

Develop

[![CI-short](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml/badge.svg?branch=devel)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci-short.yml)
[![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/devel)

### Credit

[quantarhei](https://github.com/tmancal74/quantarhei) - Concepts of aggregate construction, Hamiltonian of aggregate construction and trace over bath degrees of freedom were implemented in `quantarhei` in `python`.

[QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl) - Provides building elements for this package, such as the implementation of operators and superoperators (and many more).