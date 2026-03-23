```@meta
CurrentModule=OpenQuantumSystems
```

# Hamiltonian and Bath Model

This page summarises the Hamiltonian structure, the linear harmonic oscillator
(LHO) bath, correlation functions, and the initial-state preparation used
throughout the package. The derivations follow Chapter 1 of [Her23].

## Hamiltonian decomposition

The total Hamiltonian of a molecular aggregate coupled to a vibrational bath is
split into three parts:

```math
\hat{H} = \hat{H}_S + \hat{H}_B + \hat{H}_I,
```

where ``\hat{H}_S`` is the electronic (system) Hamiltonian, ``\hat{H}_B`` is
the bath of intramolecular vibrations, and ``\hat{H}_I`` is the
system--bath interaction.

For an aggregate of ``N`` molecules, each with a ground and excited electronic
state, the electronic Hamiltonian in the site basis reads

```math
\hat{H}_S = \sum_n \epsilon_n |n\rangle\langle n|
           + \sum_{n \neq m} J_{nm} |n\rangle\langle m|,
```

where ``\epsilon_n`` is the excitation energy of molecule ``n`` and ``J_{nm}``
is the electronic coupling between molecules ``n`` and ``m``.

## Linear harmonic oscillator bath

Each molecule ``n`` has a set of vibrational modes modelled as displaced harmonic
oscillators. The bath Hamiltonian for molecule ``n`` is

```math
\hat{H}_{B,n} = \sum_k \hbar\omega_{nk}
  \left(\hat{a}_{nk}^\dagger \hat{a}_{nk} + \tfrac{1}{2}\right),
```

and the system--bath interaction arises from the displacement of the equilibrium
position upon electronic excitation:

```math
\hat{H}_I = \sum_n \Delta\hat{V}_n |n\rangle\langle n|, \qquad
\Delta\hat{V}_n = \sum_k \hbar\omega_{nk} d_{nk}
  (\hat{a}_{nk}^\dagger + \hat{a}_{nk}),
```

where ``d_{nk}`` is the dimensionless displacement (shift) of mode ``k`` on
molecule ``n``. The Huang--Rhys factor is ``S_{nk} = d_{nk}^2 / 2``.

### Shifted vs non-shifted basis

In the **shifted (displaced) basis**, the excited-state potential energy surface
is centred at the new equilibrium. Franck--Condon factors connect ground and
excited vibrational states through the shift operator

```math
\hat{D}_n = \exp\!\Bigl(\sum_k d_{nk}
  (\hat{a}_{nk}^\dagger - \hat{a}_{nk})\Bigr).
```

In the **non-shifted basis**, the vibrational states are the same in both
electronic states, but the interaction Hamiltonian ``\hat{H}_I`` is non-zero.
OpenQuantumSystems.jl uses the non-shifted basis for numerical propagation, as
the bath trace is simpler to evaluate and delayed differential equations converge
faster.

### Reorganisation energy

Each mode contributes a reorganisation energy

```math
\lambda_{nk} = \hbar\omega_{nk} \frac{d_{nk}^2}{2}
             = \hbar\omega_{nk} S_{nk},
```

which represents the energy dissipated into the bath upon electronic excitation.
The total reorganisation energy for molecule ``n`` is
``\lambda_n = \sum_k \lambda_{nk}``.

## Correlation functions

The first-order bath correlation function for a single mode ``n`` is

```math
C_n(\tau) = \frac{\hbar^2}{2}\omega_n^2 d_n^2
  \bigl(n(\omega_n) e^{-i\omega_n\tau}
       + (n(\omega_n)+1) e^{+i\omega_n\tau}\bigr),
```

where ``n(\omega) = (e^{\hbar\omega/k_BT} - 1)^{-1}`` is the Bose--Einstein
distribution. For multiple modes on the same molecule, the correlation functions
are additive: ``C_n(\tau) = \sum_k C_{nk}(\tau)``.

The correlation function is diagonal in molecule index:
``C_{nm}(t_1, t_2) = \delta_{nm} C_n(t_1 - t_2)``.

### Higher-order correlation functions

For a Gaussian (LHO) bath, higher-order correlation functions decompose into
products of first-order ones. The second-order correlation function
(four-point function) decomposes as

```math
C_{nmkl}(t_1,t_2,t_3,t_4) = \delta_{nm}\delta_{kl} C_{nn}(t_1,t_2)C_{mm}(t_3,t_4)
  + \delta_{nk}\delta_{ml} C_{nn}(t_1,t_3)C_{mm}(t_2,t_4)
  + \delta_{nl}\delta_{mk} C_{nn}(t_1,t_4)C_{mm}(t_2,t_3).
```

This Gaussian property is essential for expressing the corrected memory kernel
entirely in terms of first-order correlation functions (see
[Corrected Memory Kernel](@ref)).

### High-temperature limit

At high temperatures, the correlation functions become real:

```math
C_{nm}(t_1,t_2) = \hbar^2\delta_{nm}\sum_u \omega_{nu}^2 d_{nu}^2
  n(\omega_{nu})\cos(\omega_{nu}(t_1-t_2)),
```

which simplifies many derivations and makes the second-order correlation
functions invariant under permutations of their time arguments.

## Lineshape function

The lineshape function ``g(t)`` encodes the bath-induced dephasing and energy
shift for a single molecule:

```math
g(t) = \sum_k S_k \bigl[(2\bar{n}_k+1)(1-\cos\omega_k t)
       + i(\sin\omega_k t - \omega_k t)\bigr],
```

where ``S_k`` is the Huang--Rhys factor and ``\bar{n}_k`` the thermal occupation.
The lineshape function and its derivatives are used by the modified Redfield
and Forster rate calculations.

## Initial state

The initial condition models ultrafast laser excitation. The bath is in thermal
equilibrium

```math
\rho_{eq} = \frac{e^{-H_S/k_BT}}{\operatorname{Tr} e^{-H_S/k_BT}}, \qquad
w_{eq} = \frac{e^{-H_B/k_BT}}{\operatorname{Tr} e^{-H_B/k_BT}},
```

and the initial density matrix of the whole system factorises as

```math
\hat{W}(t=0) = \rho_0 \otimes w_{eq},
```

where ``\rho_0 = \sum_n w_n |e_n\rangle\langle e_n|`` distributes the
population among singly excited states with weights ``w_n``.
This is justified by the Franck--Condon principle: the excitation is so fast
that the vibrational wavepacket does not change its shape.

## Reduced density matrix and relative bath part

The reduced density matrix (RDM) is obtained by tracing over the bath:

```math
\hat{\rho}(t) = \operatorname{tr}_B\{\hat{W}(t)\}.
```

The **relative bath part** (RBP) captures the remainder:

```math
\hat{w}_{nm}(t) = \frac{\langle n|\hat{W}(t)|m\rangle}
                       {\operatorname{tr}_B\{\langle n|\hat{W}(t)|m\rangle\}},
\qquad \operatorname{tr}_B\{\hat{w}_{nm}(t)\} = 1.
```

The full density matrix can be reconstructed as

```math
\hat{W}(t) = \sum_{nm} \rho_{nm}(t)\,\hat{w}_{nm}(t)\,|n\rangle\langle m|
           = \hat{\rho}(t) \oplus \hat{w}(t),
```

where ``\oplus`` denotes the formal product (not to be confused with the tensor
product ``\otimes``). This factorisation is the starting point for the
iterative treatment of the QME.
