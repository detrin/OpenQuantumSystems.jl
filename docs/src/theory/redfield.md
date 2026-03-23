```@meta
CurrentModule=OpenQuantumSystems
```

# Redfield Equations

The Redfield equations are obtained from the second-order QME by making two
key approximations: the bath remains in thermal equilibrium, and the system
density matrix varies slowly compared to the bath correlation time. The
derivation follows Section 2.2 of [Her23].

## Assumptions

1. **Weak system--bath coupling.** The excitonic basis (eigenstates of
   ``\hat{H}_S``) is the natural choice, since these are nearly stationary:

```math
\frac{\partial}{\partial t}\rho(t) = -\frac{i}{\hbar}[H_S, \rho(t)], \qquad
\rho_{\alpha\alpha}(t) \approx \text{const}, \qquad
\rho_{\alpha\beta}(t) = \rho_{\alpha\beta}(0)e^{-i\omega_{\alpha\beta}t}.
```

2. **Bath in equilibrium.** The RBP is approximated as constant:
   ``\hat{w}_{cd}^{(I)}(t) \approx \hat{w}_{eq}``, so the bath state factorises
   from the system at all times.

3. **Markov approximation.** The density matrix inside the memory integral is
   evaluated at the current time rather than at the integration variable:
   ``\hat{\rho}^{(I)}(t - \tau) \approx \hat{\rho}^{(I)}(t)``.

## Derivation

Starting from the QME for the RDM (see [Quantum Master Equation](@ref)) with
a bath in equilibrium and switching to the Schroedinger picture for the RDM:

```math
\frac{\partial}{\partial t}\hat{\rho}^{(I)}(t) =
  -\frac{i}{\hbar}[\hat{H}_S(t), \hat{\rho}(t)]
  -\frac{1}{\hbar^2}\int_0^{t}\mathrm{d}\tau\,
  U_S(t)\operatorname{tr}_B\bigl\{
  [\hat{H}_I^{(I)}(t), [\hat{H}_I^{(I)}(t-\tau),
  \hat{\rho}^{(I)}(t) \oplus \hat{w}_{eq}\,|c\rangle\langle d|]]
  \bigr\}U_S^\dagger(t).
```

Applying the Markov approximation and expanding the double commutator using
the interaction Hamiltonian
``\hat{H}_I = \sum_n \Delta\hat{V}_n \hat{K}_n``
(where ``\hat{K}_n = |n\rangle\langle n|``), the bath traces reduce to
correlation functions:

```math
\operatorname{tr}_B\{\Delta V_m(-\tau)\Delta V_n \hat{w}_{eq}\}
  = C_n^*(\tau)\,\delta_{nm}.
```

## Final form

Defining the ``\Lambda`` operator:

```math
\Lambda_n(t) = \int_0^{t}\mathrm{d}\tau\,
  C_n(\tau) U_S(\tau) \hat{K}_n U_S^\dagger(\tau),
```

the Redfield equation in the Schroedinger picture reads

```math
\frac{\partial}{\partial t}\rho(t) = -\frac{i}{\hbar}[H_S, \rho(t)]
  + \frac{1}{\hbar^2}\sum_n \Bigl[
    K_n \rho(t) \Lambda_n^\dagger(t) + \Lambda_n(t) \rho(t) K_n
    - K_n \Lambda_n(t) \rho(t) - \rho(t) \Lambda_n^\dagger(t) K_n
  \Bigr].
```

In the excitonic basis (``|\alpha\rangle, |\beta\rangle, \ldots`` being
eigenstates of ``H_S``), the matrix elements of ``\Lambda`` are

```math
\langle\alpha|\Lambda_n(t)|\beta\rangle =
  \int_0^{t}\mathrm{d}\tau\, C_n(\tau) e^{-i\omega_{\alpha\beta}\tau}
  \langle\alpha|\hat{K}_n|\beta\rangle,
```

and the Redfield equations become

```math
\frac{\partial}{\partial t}\rho(t) = -i\omega_{\alpha\beta}\rho_{\alpha\beta}
  + \frac{1}{\hbar^2}\sum_n \sum_{\gamma\delta}\bigl[
    K_{\alpha\gamma}^n \rho_{\gamma\delta}(t)(\Lambda_{\delta\beta}^n)^*(t)
    + \Lambda_{\alpha\gamma}^n(t)\rho_{\gamma\delta}(t)K_{\delta\beta}^n
    - K_{\alpha\gamma}^n \Lambda_{\gamma\delta}^n(t)\rho_{\delta\beta}(t)
    - \rho_{\alpha\gamma}(t)\Lambda_{\gamma\delta}^{n*}(t)K_{\delta\beta}^n
  \bigr],
```

where ``K_{\alpha\beta}^n = \langle\alpha|\hat{K}_n|\beta\rangle``.

## Relationship to other methods

The Redfield equations are the simplest approximate QME and serve as a
baseline:

- **Constant ansatz** (``\hat{w}^{(I)}_{ab}(t) = \hat{w}_{eq}``) gives the
  same bath approximation but keeps ``\hat{\rho}^{(I)}(\tau)`` inside the
  integral (no Markov approximation). This is the `QME_sI_ansatz` solver with
  `ansatz=:const_sch`.
- **Iterative QME** systematically improves beyond Redfield by refining the
  bath state.
- **Modified Redfield** treats diagonal fluctuations non-perturbatively.

## Implementation

In the package, `QME_sI_Redfield` solves the Redfield equations in the
interaction picture, with the memory integral evaluated numerically at each
time step. The correlation function ``C_n(\tau)`` depends on
``t - \tau`` and is constructed from the vibrational modes of each molecule.
