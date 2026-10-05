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

## Iterative QME with the corrected kernel

The corrected kernel gives an iterative scheme that never stores the RBP:

```math
\hat{w}^{0} \xrightarrow{\text{QME-RDM}} \hat{\rho}^{1}
\xrightarrow{M(\,\cdot\,;\hat w^{1}[\rho^{1}])} \hat{\rho}^{2} \to \ldots
```

``\hat\rho^1`` solves the second-order QME with the equilibrium bath; ``\hat\rho^2``
solves the QME whose kernel contains the first bath correction evaluated with
``\hat\rho^1``. Only bath correlation functions enter, so the scheme works for
infinite baths. It is implemented in [`QME_cf_iterative`](@ref) for an
[`ExcitonModel`](@ref) with one [`ExponentialCorrelation`](@ref) per site: discrete
modes ([`bath_correlation`](@ref) of a `SpectralDensity`) or a Drude-Lorentz bath
([`DrudeLorentz`](@ref)). Because every correlation function is a sum of
exponentials, all history integrals, including the nested ``t_3, t_4`` integrals
of the correction, are propagated recursively in time; the sixteen operator
orderings and the Wick pairings are generated mechanically.

Three forms of the correction are available (`correction` keyword):

| `correction` | kernel correction | weight |
|---|---|---|
| `:thesis` | as derived in the thesis, all Wick pairings | ``\rho^1_{ef}(t_4)/\rho^1_{cd}(t_3)``, acts on ``\rho_{cd}(t_2)`` |
| `:projected` | thesis + quotient-rule term: disconnected pairing removed | same |
| `:nz` | connected pairings, ``\rho_{cd}(t_2)/\rho^1_{cd}(t_3) \to 1`` | ``\rho^1_{ef}(t_4)`` |

The `:nz` form is exactly the fourth-order Nakajima-Zwanzig (time-convolution)
term ``\int \operatorname{tr}_B L(t)L(t_2)\,\mathcal{Q}\,L(t_3)L(t_4)\,
\rho^1(t_4)\otimes w_{eq}`` with ``\rho^1`` in place of ``\rho``, which makes the
first iteration a fourth-order theory.

## Checking the derivation

All statements below are verified numerically on a dimer, either against the exact
finite-bath dynamics obtained by diagonalisation or against HEOM.

- **Second-order kernel** (Eq. MF_w_0_2) and the **Wick decomposition** of the
  sixteen fourth-order terms are correct: the correlation-function kernels agree
  with the exact finite-basis kernels to ``10^{-12}`` (second order) and
  ``10^{-6}`` of the correction (first order, after the fixes of issues #102–#106).
- **Odd numbers of** ``H_I`` vanish under ``\operatorname{tr}_B`` for the
  equilibrium bath, so dropping the ``\mathcal{K}_w`` term is correct.
- **The QME-RBP misses the quotient-rule term** (see
  [Iterative Quantum Master Equation](@ref)). Consequence for the kernel: the
  thesis correction contains the disconnected pairing
  ``\langle H(t_1)H(t_2)\rangle\langle H(t_3)H(t_4)\rangle``, which evaluates to
  ``-M^0_{abcd}(t,\tau)\,\ln[\rho^1_{cd}(\tau)/\rho^1_{cd}(0)]``. It is not a
  small fourth-order correction but the second-order kernel times an
  ``\mathcal{O}(1)`` factor, and it diverges for density-matrix elements that
  start at zero.
- **The ratio** ``\rho^1_{ef}(t_4)/\rho^1_{cd}(t_3)`` is singular whenever
  ``\rho^1_{cd}`` passes through zero, e.g. for an initially unpopulated exciton.
  The `:nz` form avoids the division entirely.

## Accuracy

Relative error ``\max_t\|\rho-\rho_\text{exact}\|/\max_t\|\rho_\text{exact}-\rho(0)\|``
against the exact dynamics of a dimer with one 200 cm⁻¹ mode per molecule
(6 levels, ``J = 50`` cm⁻¹, ``\Delta E = 200`` cm⁻¹, 10 K, ``t \le 0.1``):

| Huang-Rhys ``S`` | 0.00625 | 0.0125 | 0.025 | 0.05 | 0.1 | order |
|---|---|---|---|---|---|---|
| Redfield (TCL2) | 0.053 | 0.105 | 0.21 | 0.39 | 0.67 | 1.00 |
| Born QME (``M^0``, constant ansatz) | 0.119 | 0.250 | 0.55 | 1.13 | 2.12 | 1.07 |
| IQME `:thesis` | 0.025 | 0.064 | 0.18 | 0.47 | 1.24 | 1.34 |
| IQME `:projected` | 0.0088 | 0.035 | 0.13 | 0.39 | 0.94 | 1.98 |
| IQME `:nz` | 0.0068 | 0.029 | 0.13 | 0.54 | 2.02 | 2.08 |

"Order" is the slope of the error against ``S`` between the two weakest couplings:
second-order theories have relative error ``\propto S``, a fourth-order theory
``\propto S^2``. The corrected iteration (`:projected`, `:nz`) is fourth order and
an order of magnitude more accurate than the Born QME at weak coupling; the
thesis form improves on Born but loses an order to the disconnected term. With a
single undamped mode per molecule the memory never decays, so every
perturbative scheme fails at ``S \gtrsim 0.05``.
