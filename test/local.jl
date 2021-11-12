using BenchmarkTools
using Plots
using ProgressBars
using OpenQuantumSystems
# using OpenQuantumSystemsPrivate
using LinearAlgebra
using Random
using QuadGK
using ProgressMeter
using Expokit
# using DifferentialEquations

Random.seed!(0)

import OrdinaryDiffEq, DiffEqCallbacks, DelayDiffEq
import SparseArrays: sparse
import QuantumOpticsBase

D(op1::Array, op2::Array) = abs(norm(op1 - op2))
D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
D(op1::AbstractOperator, op2::AbstractOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))
D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) = abs(tracedistance_nh(dense(op1), dense(op2)))

function master_ansatz2(
    rho0::T,
    tspan::Array,
    p;
    reltol::Float64 = 1.0e-12,
    abstol::Float64 = 1.0e-12,
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    # (du,u,h,p,t)
    Ham_S, Ham_int, H_lambda, H_S, H_Sinv, W0, W0_bath, agg, FCProd, aggIndices, vibindices, groundState, elementtype = p
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    tmp3 = copy(W0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) =
        dmaster_ansatz2(t, rho, drho, history_fun, tmp1, tmp2, tmp3, p, history_fun)
    # tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    tspan_ = deepcopy(tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    OpenQuantumSystems.integrate_delayed(
        tspan_,
        dmaster_,
        history_fun,
        x0,
        state,
        dstate,
        fout;
        p = p,
        reltol = reltol,
        abstol = abstol,
        alg = alg,
        kwargs...,
    )
end

function dmaster_ansatz2(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp1::Array,
    tmp2::Array,
    tmp3::Array,
    p,
    h
) where {B<:Basis,T<:Operator{B,B}}
    println(t)
    Ham_S, Ham_int, H_lambda, H_S, H_Sinv, W0, W0_bath, agg, FCProd, aggIndices, vibindices, groundState, elementtype = p
    # Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    Ham_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
    kernel_integrated, err = QuadGK.quadgk(
        s -> MemoryKernel2(t, s, tmp1, tmp2, tmp3, h, p, Ham_II_t),
        0,
        t,
        rtol = 1e-8,
    )

    LinearAlgebra.mul!(tmp1, Ham_II_t, W0_bath.data, -elementtype(im), zero(elementtype))
    LinearAlgebra.mul!(tmp1, W0_bath.data, Ham_II_t, elementtype(im), one(elementtype))
    K_traced = trace_bath(tmp1, agg, FCProd, aggIndices, vibindices; groundState = groundState)
    rho0 = trace_bath(W0, agg, FCProd, aggIndices, vibindices; groundState = groundState)
    println(t)
    println(K_traced)
    println(K_traced .* rho0.data)

    # println(data)
    LinearAlgebra.mul!(drho.data, one(elementtype), K_traced, -elementtype(im), zero(elementtype))
    data = kernel_integrated .* rho.data
    LinearAlgebra.mul!(drho.data, one(elementtype), data, -elementtype(1), one(elementtype))
    # println(-1im * (K_traced .* rho.data))
    # println(kernel_integrated .* rho.data)

    return drho
end

function MemoryKernel2(t, s, tmp1, tmp2, tmp3, h, p, Ham_II_t)
    Ham_S, Ham_int, H_lambda, H_S, H_Sinv, W0, W0_bath, agg, FCProd, aggIndices, vibindices, groundState, elementtype = p
    Ham_II_s = getInteractionHamIPictureA(Ham_S, Ham_int, s)
    Ham = Ham_S + Ham_int
    LinearAlgebra.mul!(tmp1, Ham.data, W0_bath.data, -elementtype(im)*s, zero(elementtype))
    LinearAlgebra.mul!(tmp1, W0_bath.data, Ham.data, elementtype(im)*s, one(elementtype))
    LinearAlgebra.mul!(tmp1, one(elementtype), W0_bath.data, elementtype(1), one(elementtype))
    # tmp1 += W0_bath.data


    # LinearAlgebra.mul!(tmp1, H_S, exp(-1im * H_lambda * s), one(elementtype), zero(elementtype))
    # LinearAlgebra.mul!(tmp3, tmp1, H_Sinv, one(elementtype), zero(elementtype))

    # LinearAlgebra.mul!(tmp1, adjoint(tmp3), Ham_int.data, one(elementtype), zero(elementtype))
    # LinearAlgebra.mul!(tmp2, tmp1, tmp3, one(elementtype), zero(elementtype))
    # println(tmp2)
    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end
    tmp2 .= ad(rho, tmp1, agg, FCProd, aggIndices, vibindices; groundState = groundState)
    tmp3 .= getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, s)

    QuantumOpticsBase.mul!(tmp1, tmp3, tmp2, elementtype(1), zero(elementtype))
    QuantumOpticsBase.mul!(tmp1, tmp2, tmp3, -elementtype(1), one(elementtype))

    QuantumOpticsBase.mul!(tmp2, Ham_II_t, tmp1, elementtype(1), zero(elementtype))
    QuantumOpticsBase.mul!(tmp2, tmp1, Ham_II_t, -elementtype(1), one(elementtype))
    MK_traced = trace_bath(tmp2, agg, FCProd, aggIndices, vibindices; groundState = groundState)
    println(MK_traced)
    return MK_traced
end

HR = 0.01
shift = (2.0 * HR)
modes = [Mode(180., shift)]
mols = [
    Molecule([Mode(300., shift)], 2, [12500., 12670.]),
    Molecule([Mode(300., shift)], 2, [12510., 12730.])
]

agg = Aggregate(mols)
for mol_i in 2:length(agg.molecules)
    agg.coupling[mol_i, mol_i+1] = 200
    agg.coupling[mol_i+1, mol_i] = 200
end
# agg.coupling[1, 3] = 50
# agg.coupling[3, 1] = 50
# agg.coupling[:, :] = rand(Float64, (length(mols)+1, length(mols)+1)) * 200

println("getIndices")
@time aggIndices = getIndices(agg; groundState = true)
aggInds = getIndices(agg; groundState = true)
vibindices = getVibIndices(agg, aggInds)
aggIndsLen = length(aggInds)
basis = GenericBasis([aggIndsLen])
FCFact = getFranckCondonFactors(agg, aggInds; groundState = true)
FCProd = getFCProd(agg, FCFact, aggIndices, vibindices; groundState = true)
Ham = getAggHamiltonian(agg, aggInds, FCFact; groundState = true, groundEnergy=false)

Ham_bath = getAggHamiltonianBath(agg)
Ham_sys = getAggHamiltonianSystem(agg; groundState = false)
b_sys = GenericBasis([size(Ham_sys, 1)])
b_bath = GenericBasis([size(Ham_bath, 1)])

Ham_int = getAggHamiltonianInteraction(agg, aggIndices, FCFact; groundState=true)
Ham_S = Ham - Ham_int
println(agg.coupling)

println(Ham)

println("###########")

t_max = 0.05
t_count = 200
t0 = 0.
t_step = (t_max - t0) / (t_count)
tspan = [t0:t_step:t_max;]

T = 1e-1
mu_array = [[2, 1]]
W0 = thermal_state(T, mu_array, Ham, vibindices, aggIndices; diagonalize = true)
println(W0)
W0_bath = get_rho_bath(W0, agg, FCProd, aggIndices, vibindices; groundState=true)

H_lambda, H_S = eigen(Ham_S.data)
H_Sinv = inv(H_S)
H_lambda = diagm(H_lambda)
println("")
