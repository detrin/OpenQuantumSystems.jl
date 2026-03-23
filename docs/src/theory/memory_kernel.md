```@meta
CurrentModule=OpenQuantumSystems
```

# Corrected Memory Kernel

The standard Redfield equations use a memory kernel that assumes the bath
remains in thermal equilibrium. By deriving the first correction to the bath
state (RBP) and substituting it into the memory kernel, we obtain a
**corrected Redfield equation** that accounts for the bath's response to the
system evolution. The derivation follows Section 3.7 of [Her23].

## Perturbation of the RBP in the memory kernel

The memory kernel in the QME-RDM is

```math
M_{abcd}(t,\tau) = \langle a|\bigl[
  \hat{H}_I^{(I)}(t), [\hat{H}_I^{(I)}(\tau),
  \hat{\rho}^{(I)}(\tau) \oplus \hat{w}^{(I)}(\tau)]
\bigr]|b\rangle,
```

which depends on the bath state ``\hat{w}^{(I)}(\tau)``. In the standard
(zeroth-order) Redfield equations, ``\hat{w}^{(I)}(\tau) = \hat{w}_{eq}``
(constant). The corrected memory kernel uses the first-order bath correction
``\hat{w}^{1,(I)}`` from the iterative scheme.

## Zeroth-order memory kernel

With ``\hat{w}^{0,(I)} = \hat{w}_{eq}``:

```math
M_{abcd}(t_1,t_2; \hat{w}^{0,(I)}) =
  \sum_e \bigl[
    e^{i\omega_{ae}t_1 + i\omega_{ec}t_2} C_{aeec}(t_1,t_2)\delta_{db}
    + e^{i\omega_{de}t_1 + i\omega_{eb}t_2} C_{deeb}(t_1,t_2)\delta_{ac}
    - e^{i\omega_{ae}t_2 + i\omega_{db}t_2} C_{dbac}(t_2,t_1)
    - e^{i\omega_{ae}t_2 + i\omega_{db}t_1} C_{dbac}(t_1,t_2)
  \bigr],
```

which involves only first-order correlation functions.

## First-order correction

The corrected memory kernel ``M_{abcd}(t_1,t_2; \hat{w}^{1,(I)})`` contains
the zeroth-order term plus a correction involving a double integral over
``t_3`` and ``t_4``:

```math
M_{abcd}(t_1,t_2; \hat{w}^{1,(I)}) = M_{abcd}(t_1,t_2; \hat{w}^{0,(I)})
  - \frac{1}{\hbar^2}\int_0^{t_2}\mathrm{d}t_3 \int_0^{t_3}\mathrm{d}t_4\,
  \bigl[M_{1a} - M_{1b} + M_{1c} - M_{1d}
       + M_{2a} - M_{2b} + M_{2c} - M_{2d}
       + M_{3a} - M_{3b} + M_{3c} - M_{3d}
       + M_{4a} - M_{4b} + M_{4c} - M_{4d}\bigr],
```

where each ``M_{ia}`` etc. involves four interaction Hamiltonians under the
bath trace. The key result is that for an LHO bath, all these terms can be
expressed as **products of first-order correlation functions**, thanks to the
Gaussian decomposition property.

## Gaussian decomposition

For a linear harmonic oscillator bath, the second-order correlation function
decomposes as:

```math
C_{nmkl}(t_1,t_2,t_3,t_4) = \delta_{nm}\delta_{kl}C_{nn}(t_1,t_2)C_{mm}(t_3,t_4)
  + \delta_{nk}\delta_{ml}C_{nn}(t_1,t_3)C_{mm}(t_2,t_4)
  + \delta_{nl}\delta_{mk}C_{nn}(t_1,t_4)C_{mm}(t_2,t_3).
```

This allows the full corrected memory kernel (Eq. 3.64 in the thesis) to be
written entirely in terms of first-order correlation functions, multi-frequency
phase factors ``\Omega_{ab,cd,ef,gh}``, and basis transformation coefficients
``K_{ab,cd}^{nm}``.

## High-temperature limit

At high temperatures, the correlation functions become real, and the five
distinct groups of second-order correlation functions in the memory kernel
simplify. The corrected memory kernel in the high-temperature limit (Eq. 3.76)
is more compact but retains the same structure of corrections.

## Corrected Redfield iteration

The corrected memory kernel enables a new iterative scheme that avoids
explicitly computing the RBP:

```math
\hat{w}^{0,(I)}(t) \xrightarrow{\text{QME-RDM}} \hat{\rho}^{1,(I)}(t)
\xrightarrow{\text{QME-RDM}} \hat{\rho}^{2,(I)}(t)
\xrightarrow{\text{QME-RDM}} \hat{\rho}^{3,(I)}(t) \to \ldots
```

At each step, the corrected memory kernel
``M_{abcd}(t_1,t_2; \hat{\rho}^{k,(I)})`` from the previous iteration is used.
The RBP never needs to be explicitly stored or propagated -- it is implicit in
the corrected kernel. This makes the approach feasible for infinite baths
(via the spectral density function), where the RBP cannot be represented
explicitly.

## Limitations

While the corrected memory kernel provides a principled improvement over
standard Redfield, numerical experiments show:

- The improvement is modest for finite systems -- typically less than one
  order of magnitude in the scoring function.
- The corrected kernel is significantly more complex (involving double
  integrals over ``t_3, t_4``), increasing computational cost.
- For the population-only case (diagonal elements in the excitonic basis),
  two of the sixteen terms in the correction vanish, simplifying the
  expression.

The approach is most promising for weak interactions where the spectral
density can be modelled by many weakly-coupled modes.
