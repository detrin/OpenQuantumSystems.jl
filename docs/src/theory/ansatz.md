```@meta
CurrentModule=OpenQuantumSystems
```

# Bath State Ansatze

The second-order QME contains the full density matrix
``\hat{W}^{(I)}(\tau) = \hat{\rho}^{(I)}(\tau) \oplus \hat{w}^{(I)}(\tau)``
inside the memory integral. To avoid propagating the full system+bath state,
we can replace the relative bath part (RBP) ``\hat{w}^{(I)}(\tau)`` with an
**ansatz** -- an approximate form that is cheaper to evaluate. Different
ansatze offer different accuracy/cost trade-offs. The derivations follow
Section 3.3 of [Her23].

## Constant ansatz

The simplest approximation freezes the bath at its initial equilibrium state:

```math
\hat{w}_{ab}^{(I)}(t) = \hat{w}_{eq}
```

for all electronic states ``a, b`` and all times ``t``. This is justified when
the bath is infinite and initially in equilibrium -- any effect of the system
on the bath is negligible. For finite systems, this is accurate only for short
times (typically ``t < 400`` fs for weak coupling).

In the package: `QME_sI_ansatz(W0, tspan, agg; ansatz=:const_sch)` or
`:const_int`.

## Linear ansatz (L1)

For short times, a Taylor expansion of the evolution operators yields a
first-order correction:

```math
\hat{w}(t) = \hat{w}_{eq} - \frac{i}{\hbar}[\hat{H}, \hat{w}_{eq}]\,t.
```

In the interaction picture:

```math
\hat{w}_{ab}^{(I)}(t) = \langle a|U_0^\dagger(t)
  \bigl[\hat{w}_{eq} - \tfrac{i}{\hbar}[\hat{H}, \hat{w}_{eq}]\,t\bigr]
  U_0(t)|b\rangle.
```

This is more accurate than the constant ansatz for short times but the
polynomial correction diverges for long simulations.

In the package: `ansatz=:linear_sch`.

## Piecewise-linear ansatz (L2)

To extend the validity of the linear correction, the evolution is broken into
``k`` steps of length ``t/k``:

```math
\hat{w}^{(I)}(t) = \mathcal{U}_0^\dagger(t)\,
  \tilde{\mathcal{U}}^k\!\left(\frac{t}{k}\right)\,\hat{w}(0),
```

where ``\tilde{\mathcal{U}}(t)\hat{O} = \hat{O} - \frac{i}{\hbar}[\hat{H}, \hat{O}]\,t``
is the linearised evolution superoperator. With more corrections (larger ``k``),
the L2 ansatz approaches the exact evolution, but divergence is still possible
for ``k`` too small.

In the package: `ansatz=:linear2_sch`.

## U1 ansatz (block-diagonal unitary evolution)

Instead of a polynomial approximation, we evolve only the **diagonal**
(population) blocks of the RBP with the corresponding block Hamiltonian:

```math
\hat{w}^{(I)}(t) = \mathcal{U}_0^\dagger(t)
  \left[\sum_a \hat{U}_a(t)\,\hat{w}_{eq}\,\hat{U}_a^\dagger(t)\,|a\rangle\langle a|
  + \sum_{a \neq b} \hat{w}_{eq}\,|a\rangle\langle b|\right],
```

where ``\hat{U}_a(t) = \exp\bigl[-\frac{i}{\hbar}\langle a|\hat{H}|a\rangle\,t\bigr]``
is the evolution operator for the ``a``-th population block of the Hamiltonian.
Off-diagonal blocks of the RBP are left at equilibrium.

This ansatz correctly captures the separate evolution of each population's bath
state without the divergence issues of polynomial approximations.

In the package: `ansatz=:upart1_sch` or `:upart1_int`.

## U2 ansatz (full-block unitary evolution)

The U2 ansatz extends U1 by evolving **all** blocks of the RBP, including
coherences:

```math
\hat{w}^{(I)}(t) = \mathcal{U}_0^\dagger(t)
  \sum_{ab} \hat{U}_{ab}(t)\,\hat{w}_{eq}\,\hat{U}_{ab}^\dagger(t),
```

where ``\hat{U}_{ab}(t) = \exp\bigl[-\frac{i}{\hbar}\langle a|\hat{H}|b\rangle\,t\bigr]``
uses the full ``(a,b)`` block of the Hamiltonian. Numerically, U1 and U2 give
nearly identical results for the systems tested, suggesting that the
off-diagonal bath evolution is less important than the diagonal part.

In the package: `ansatz=:upart2_sch` or `:upart2_int`.

## Comparison of ansatze

| Ansatz | Stability | Accuracy | Cost |
|--------|-----------|----------|------|
| Constant | Excellent | Low (short times only) | Lowest |
| L1 (linear) | Limited (diverges) | Good for short times | Low |
| L2 (piecewise-linear) | Better with large ``k`` | Improves with ``k`` | ``O(k)`` |
| U1 (block-diagonal) | Good | Moderate | Moderate |
| U2 (full-block) | Good | ``\approx`` U1 | Moderate |
| QME test (exact bath) | Excellent | Exact (reference) | Highest |

The constant ansatz performs similarly to Redfield equations. The L1, U1, and
U2 ansatze do not systematically outperform Redfield for finite systems,
though L2 with many corrections can. For best results in the weak-coupling
regime, use the **iterative ansatz** (see [Iterative Quantum Master Equation](@ref)) which
systematically improves the bath state.
