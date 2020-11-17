# OpenQuantumSystems.jl

| Service  | Master  | Develop  |
| :------- | :------ | :------- |
| Stable | [![Actions Status](https://github.com/detrin/OpenQuantumSystems.jl/workflows/CI/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions) | [![Actions Status](https://github.com/detrin/OpenQuantumSystems.jl/workflows/CI/badge.svg?branch=devel)](https://github.com/detrin/OpenQuantumSystems.jl/actions) |
| Nightly Julia | [![Actions Status](https://github.com/detrin/OpenQuantumSystems.jl/workflows/CI-nightly-julia/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions) | [![Actions Status](https://github.com/detrin/OpenQuantumSystems.jl/workflows/CI-nightly-julia/badge.svg?branch=devel)](https://github.com/detrin/OpenQuantumSystems.jl/actions) |
| Nigtly deploy | [![Actions Status](https://github.com/detrin/OpenQuantumSystems.jl/workflows/Deploy%20Nightly/badge.svg?branch=master)](https://github.com/detrin/OpenQuantumSystems.jl/actions) | [![Actions Status](https://github.com/detrin/OpenQuantumSystems.jl/workflows/Deploy%20Nightly/badge.svg?branch=devel)](https://github.com/detrin/OpenQuantumSystems.jl/actions) |
| coveralls | [![Coverage Status](https://coveralls.io/repos/JuliaLang/Example.jl/badge.svg?branch=master)](https://coveralls.io/r/detrin/OpenQuantumSystems.jl?branch=master) | [![Coverage Status](https://coveralls.io/repos/github/detrin/OpenQuantumSystems.jl/badge.svg?branch=devel)](https://coveralls.io/github/detrin/OpenQuantumSystems.jl?branch=devel) |
| codecov | [![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/detrin/OpenQuantumSystems.jl) | [![codecov](https://codecov.io/gh/detrin/OpenQuantumSystems.jl/branch/devel/graph/badge.svg)](https://codecov.io/gh/detrin/OpenQuantumSystems.jl) |


**OpenQuantumSystems.jl** is a numerical framework written in [Julia] that makes it easy to simulate various kinds of quantum systems. It is inspired by the [quantarhei](https://github.com/tmancal74/quantarhei) and [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

### Questions & Contributions

The package is still under development. If you have any questions or need help, see [gitter channel](https://gitter.im/OpenQuantumSystems-jl/community) and ask away. Also, contributions of any kind are always welcome! If you have any ideas for bug fixes, new features, optimisation or unit tests suggest it straight away or create a pull request.

### Roadmap
1. * Efficient representation of operators given Hamiltonians in different Hilbert basis.
   * Vibrational basis, spin basis, exciton basis.
   * Calculation of Franck-Condon factors for multidimensional linear harmonic oscilators.
   * Base transformation.
2. * Schrodinger equation for solving dynamics.
   * Solving dynamics for mixed states using exponentials, Liouville equation, QME and their modifications.
   * Decomposition of mixed states into a linear combination of pure states using linear programming.
   * Loading Hamiltonian and data storage.
3. Implementing Foerster and modified Redfield theories of excitation energy transfer.
4. * GPU support for Schrodinger equation and possibly state decomposition.
   * Anharmonic oscillators.
   * Double excited states.

  