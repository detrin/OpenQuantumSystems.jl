```@meta
CurrentModule=OpenQuantumSystems
```

# Quantum Master Equation

This page derives the second-order quantum master equation (QME) used as the
foundation for all solvers in the package. The derivation follows Chapter 2 of
[Her23].

## From Liouville--von Neumann to QME

The starting point is the Liouville--von Neumann (LvN) equation for the full
system+bath density matrix ``\hat{W}(t)``:

```math
\frac{\partial}{\partial t}\hat{W}(t) = -\frac{i}{\hbar}[\hat{H}, \hat{W}(t)]
  = -\frac{i}{\hbar}[\hat{H}_0, \hat{W}(t)]
    -\frac{i}{\hbar}[\hat{H}_I, \hat{W}(t)],
```

where ``\hat{H}_0 = \hat{H}_S + \hat{H}_B`` and ``\hat{H}_I`` is the
system--bath interaction.

Moving to the **interaction picture** with respect to ``\hat{H}_0``:

```math
\frac{\partial}{\partial t}\hat{W}^{(I)}(t)
  = -\frac{i}{\hbar}[\hat{H}_I^{(I)}(t), \hat{W}^{(I)}(t)],
```

where ``\hat{H}_I^{(I)}(t) = U_0^\dagger(t)\hat{H}_I U_0(t)`` and
``U_0(t) = e^{-i\hat{H}_0 t/\hbar}``.

## Formal integration

Integrating once gives

```math
\hat{W}^{(I)}(t) = \hat{W}^{(I)}(t_0) - \frac{i}{\hbar}\int_{t_0}^{t}
  \mathrm{d}\tau\,[\hat{H}_I^{(I)}(\tau), \hat{W}^{(I)}(\tau)].
```

Substituting back into the LvN equation yields the **second-order QME** for
system and bath:

```math
\frac{\partial}{\partial t}\hat{W}^{(I)}(t) =
  -\frac{i}{\hbar}[\hat{H}_I^{(I)}(t), \hat{W}^{(I)}(t_0)]
  -\frac{1}{\hbar^2}\int_{t_0}^{t}\mathrm{d}\tau\,
  [\hat{H}_I^{(I)}(t), [\hat{H}_I^{(I)}(\tau), \hat{W}^{(I)}(\tau)]].
```

## QME for the reduced density matrix

Applying the trace over bath degrees of freedom to both sides, and using the
factorised initial condition ``\hat{W}(0) = \rho_0 \otimes w_{eq}`` so that
``\operatorname{tr}_B\{[\hat{H}_I^{(I)}(t), \hat{W}^{(I)}(0)]\} = 0``,
we obtain

```math
\frac{\partial}{\partial t}\hat{\rho}^{(I)}(t) =
  -\frac{1}{\hbar^2}\int_0^{t}\mathrm{d}\tau\,
  \operatorname{tr}_B\bigl\{
  [\hat{H}_I^{(I)}(t), [\hat{H}_I^{(I)}(\tau), \hat{W}^{(I)}(\tau)]]
  \bigr\}.
```

This is exact to second order in ``\hat{H}_I``. The key difficulty is that
``\hat{W}^{(I)}(\tau)`` still appears inside the integral -- the system and
bath remain coupled. Different approximations for ``\hat{W}^{(I)}(\tau)``
lead to different solver families.

## Superoperator form

Using the Liouville superoperator
``\mathcal{L}(t)\hat{O} = [\hat{H}_I^{(I)}(t), \hat{O}]``, the QME can be
written compactly as

```math
\hat{W}^{(I)}(t) = \exp_\rightarrow\!\left(
  -\frac{i}{\hbar}\int_{t_0}^{t}\mathrm{d}\mathcal{L}_I^{(I)}(\tau)
\right)\hat{W}^{(I)}(t_0),
```

where ``\exp_\rightarrow`` denotes the time-ordered exponential. Expanding this
gives the Dyson series:

```math
\hat{W}^{(I)}(t) = \hat{W}^{(I)}(t_0)
  - \frac{i}{\hbar}\int_{t_0}^{t}\mathrm{d}t_1\,
    [\hat{H}_I^{(I)}(t_1), \hat{W}^{(I)}(t_1)]
  - \frac{1}{\hbar^2}\int_{t_0}^{t}\mathrm{d}t_1
    \int_{t_0}^{t_1}\mathrm{d}t_2\,
    [\hat{H}_I^{(I)}(t_1), [\hat{H}_I^{(I)}(t_2), \hat{W}^{(I)}(t_2)]]
  + \ldots
```

Truncation at second order and tracing over the bath yields the QME used
throughout the package.

## Solving QME as a delayed differential equation

The QME is an integro-differential equation (IDE): the right-hand side depends
on the history ``\hat{W}^{(I)}(\tau)`` for all ``0 \le \tau \le t``. Direct
numerical integration (e.g. fourth-order Runge--Kutta with stored history) is
possible but suffers from a "cold start" problem at early times.

A more robust approach rewrites the IDE as a **delayed differential equation**
(DDE). At each time step ``t``, the memory integral is approximated by a
function ``F`` that uses previously computed values
``\hat{W}^{(I)}(\tau)`` on ``[0, t)``:

```math
\frac{\partial}{\partial t}\hat{W}^{(I)}(t) =
  -\frac{i}{\hbar}[\hat{H}_I^{(I)}(t), \hat{W}^{(I)}(0)]
  -\frac{1}{\hbar^2}F\!\left(t; \hat{W}^{(I)}(\tau),\, 0 \le \tau < t\right),
```

where ``F`` is evaluated using Gauss--Kronrod quadrature (QuadGK.jl). The DDE
is then solved with a standard adaptive ODE solver (Tsit5 from
DifferentialEquations.jl). This approach handles the cold-start problem
gracefully and controls the integration error through tolerance parameters
`int_reltol` and `int_abstol`.

## Interaction picture vs Schroedinger picture

The QME is naturally formulated in the interaction picture, where the fast
free evolution under ``\hat{H}_0`` is factored out, leaving only the slow
bath-induced dynamics. Results can be converted to the Schroedinger picture via

```math
\hat{\rho}(t) = U_S(t)\hat{\rho}^{(I)}(t)U_S^\dagger(t).
```

In the package, the `_sI` suffix denotes reduced dynamics in the interaction
picture, `_sS` in the Schroedinger picture, and `_SI`/`_SS` for the full
density matrix.
