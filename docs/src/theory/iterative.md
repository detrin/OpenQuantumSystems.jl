```@meta
CurrentModule=OpenQuantumSystems
```

# Iterative Quantum Master Equation

The iterative quantum master equation (IQME) systematically improves the
reduced density matrix by alternating between solving for the system (RDM) and
the bath (RBP). The derivation follows Section 3.5 of [Her23].

## Motivation

The direct ansatze (constant, linear, U1, U2) approximate the bath state
without accounting for the system's back-action on the bath. The iterative
approach solves the RDM and RBP **sequentially**: first compute the RDM with
an approximate bath, then use the RDM to improve the bath, and repeat.

## Separating QME into RDM and RBP equations

Using the factorisation ``\hat{W}(t) = \hat{\rho}(t) \oplus \hat{w}(t)``,
the QME splits into coupled equations:

**QME-RDM** (for the reduced density matrix):

```math
\frac{\partial}{\partial t}\rho_{ab}^{(I)}(t) =
  -\frac{i}{\hbar}\operatorname{tr}_B\!\left\{
    \bigl[\hat{H}_I^{(I)}(t), \hat{W}^{(I)}(0)\bigr]_{ab}
  \right\}
  -\frac{1}{\hbar^2}\int_0^{t}\mathrm{d}\tau\,
  \operatorname{tr}_B\!\left\{
    \bigl[\hat{H}_I^{(I)}(t),
    [\hat{H}_I^{(I)}(\tau), \hat{\rho}^{(I)}(\tau) \oplus \hat{w}^{(I)}(\tau)]
    \bigr]_{ab}
  \right\}.
```

**QME-RBP** (for the relative bath part):

```math
\frac{\partial}{\partial t}\hat{w}_{ab}^{(I)}(t) =
  -\frac{i}{\hbar\rho_{ab}^{(I)}(t)}
  \bigl[\hat{H}_I^{(I)}(t), \hat{W}^{(I)}(0)\bigr]_{ab}
  -\frac{1}{\hbar^2\rho_{ab}^{(I)}(t)}\int_0^{t}\mathrm{d}\tau\,
  \bigl[\hat{H}_I^{(I)}(t),
  [\hat{H}_I^{(I)}(\tau), \hat{\rho}^{(I)}(\tau) \oplus \hat{w}^{(I)}(\tau)]
  \bigr]_{ab}.
```

## Iterative scheme

Starting from the constant ansatz ``\hat{w}^{0,(I)}(t) = \hat{w}_{eq}``:

```math
\hat{w}^{0,(I)}(t) \xrightarrow{\text{QME-RDM}} \hat{\rho}^{1,(I)}(t)
\xrightarrow{\text{QME-RBP}} \hat{w}^{1,(I)}(t)
\xrightarrow{\text{QME-RDM}} \hat{\rho}^{2,(I)}(t)
\xrightarrow{\text{QME-RBP}} \hat{w}^{2,(I)}(t) \to \ldots
```

At each step, we solve a single IDE using the results from the previous
iteration. Convergence is not guaranteed, but numerical experiments show
improvement for weak interactions.

## I ansatz (exact iterative)

The ``(k+1)``-th corrected RBP uses the ``k``-th RDM and RBP:

```math
\hat{w}_{ab}^{k+1,(I)}(t) = \hat{w}_{ab}^{0,(I)}(0)
  + \int_0^{t}\mathrm{d}t_1\,\left[
  -\frac{i}{\hbar}\mathcal{K}_w(t_1)
  - \frac{1}{\hbar^2}\sum_{cd}\int_0^{t_1}\mathrm{d}t_2\,
  M_{w,abcd}(t_1,t_2)\frac{\rho_{cd}^{k,(I)}(t_2)}{\rho_{ab}^{k,(I)}(t_1)}
  \hat{w}_{cd}^{k,(I)}(t_2)
  \right].
```

The corresponding RDM equation uses the newly computed RBP:

```math
\frac{\partial}{\partial t}\rho_{ab}^{k+1,(I)}(t) =
  -\frac{i}{\hbar}\mathcal{K}_{ab}(t)
  - \frac{1}{\hbar^2}\sum_{cd}\int_0^{t}\mathrm{d}\tau\,
  M_{abcd}(t,\tau;\hat{w}^{k,(I)})\,\rho_{cd}^{k+1,(I)}(\tau).
```

In the package: `QME_sI_iterative(W0, rho_0, W_0_bath, tspan, agg)`.

## I.M1 ansatz (first Markov approximation)

The I.M1 ansatz simplifies the I ansatz by assuming that the RDM evolves
slowly in the interaction picture: ``\rho^{(I)}(t_1) \approx \rho^{(I)}(t)``.
This eliminates the need to store intermediate RDM values inside the RBP
integral.

The RBP equation becomes:

```math
\hat{w}_{ab}^{k+1,(I)}(t) = \hat{w}_{ab}^{k,(I)}(0)
  - \frac{i}{\hbar}\int_0^{t}\frac{1}{\rho_{ab}^{k,(I)}(t)}\mathrm{d}t_1\,
  \langle a|[\hat{H}_I^{(I)}(t_1), \hat{W}(0)]|b\rangle
  - \frac{1}{\hbar^2}\sum_{cd}\int_0^{t}\mathrm{d}t_1
  \int_0^{t_1}\mathrm{d}t_2\,
  \langle a|[\hat{H}_I^{(I)}(t_1),
  [\hat{H}_I^{(I)}(t_2), \hat{w}_{cd}^{k,(I)}(t_1)\,|c\rangle\langle d|]]|b\rangle
  \frac{\rho_{cd}^{k,(I)}(t_2)}{\rho_{ab}^{k,(I)}(t)}.
```

The I.M1 and I.M2 ansatze give the same results in the first iteration
(since ``\hat{w}^{0,(I)}`` is constant). In subsequent iterations, I.M1 is
preferred because it does not require intermediate RBP calculations.

In the package: `QME_sI_iterative(W0, rho_0, W_0_bath, tspan, agg; method=:markov0)`.

## I.M2 ansatz (second Markov approximation)

Additionally assumes that the RBP itself varies slowly:
``\hat{w}^{k,(I)}(t_2) \approx \hat{w}^{k,(I)}(t_1)``.

In the package: `QME_sI_iterative(W0, rho_0, W_0_bath, tspan, agg; method=:markov1)`.

## Numerical results

For a two-molecule aggregate with Huang--Rhys factors ``S_1 = S_2 = 0.05``:

- The **first iteration** of I.M1 already improves over Redfield and the
  constant ansatz, diverging from the exact solution later (around 300 fs
  vs 200 fs for the constant ansatz).
- The **second iteration** continues to improve, demonstrating that the
  iterative approach is genuinely converging.
- For stronger interactions (``S_1 = S_2 = 0.1``), I.M1 still outperforms
  Redfield in the first iteration, but the second iteration shows signs of
  divergence.

## Scoring function

To compare ansatze quantitatively, the package uses the scoring function:

```math
\text{score}(\rho_{nm}) = \log_{10}\!\left(
  \frac{1}{N}\sum_i \frac{|\rho_{nm}(t_i)| - |\rho_{ref,nm}(t_i)|}
                        {|\rho_{ref,nm}(t_i)| + c_{\text{smooth}}}
\right),
```

where ``c_{\text{smooth}} = 10^{-9}``. More negative scores indicate better
agreement with the exact solution. A score near zero corresponds to about
10% average relative error.
