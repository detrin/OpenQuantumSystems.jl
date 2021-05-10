# OpenQuantumSystems.jl

[![Join the chat at https://gitter.im/OpenQuantumSystems-jl/community](https://badges.gitter.im/OpenQuantumSystems-jl/community.svg)](https://gitter.im/OpenQuantumSystems-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CI](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/detrin/OpenQuantumSystems.jl)

**OpenQuantumSystems.jl** is a numerical framework written in [Julia] that makes
it easy to simulate various kinds of open quantum systems. It is inspired by the
[quantarhei](https://github.com/tmancal74/quantarhei) and
[QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

### Questions & Contributions

The package is still under development. If you have any questions or need help,
see [gitter channel](https://gitter.im/OpenQuantumSystems-jl/community) and ask
away. Also, contributions of any kind are always welcome! If you have any ideas
for bug fixes, new features, optimisation or unit tests suggest it straight away
or create a pull request.

### Roadmap

1. - Efficient representation of operators given Hamiltonians in different
     Hilbert basis.
   - Vibrational basis, spin basis, exciton basis.
   - Calculation of Franck-Condon factors for multidimensional linear harmonic
     oscilators.
   - Base transformation.
2. - Schrodinger equation for solving dynamics.
   - Solving dynamics for mixed states using exponentials, Liouville equation,
     QME and their modifications.
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

