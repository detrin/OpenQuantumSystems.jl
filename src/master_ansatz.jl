
using OpenQuantumSystems
import QuantumOpticsBase, LinearAlgebra, OrdinaryDiffEq, QuadGK, DelayDiffEq

function master_ansatz(
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
    Ham_S, Ham_int, H_lambda, H_S, H_Sinv, Ham_B, W0, W0_bath, agg, FCProd, aggIndices, vibindices, elementtype = p
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    tmp3 = copy(W0.data)
    dmaster_(t, rho::T, drho::T, history_fun, p) =
        dmaster_ansatz(t, rho, drho, history_fun, tmp1, tmp2, tmp3, p, history_fun)
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

function dmaster_ansatz(
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
    Ham_S, Ham_int, H_lambda, H_S, H_Sinv, Ham_B, W0, W0_bath, agg, FCProd, aggIndices, vibindices, elementtype = p
    # Ham_II_t = getInteractionHamIPicture(Ham_S, Ham_int, t)
    Ham_II_t = getInteractionHamIPictureA(Ham_int.data, H_lambda, H_S, H_Sinv, t)
    kernel_integrated, err = QuadGK.quadgk(
        s -> MemoryKernel(t, s, tmp1, tmp2, tmp3, h, p, Ham_II_t),
        0,
        t,
        rtol = 1e-3,
        atol = 1e-3
    )
    
    LinearAlgebra.mul!(tmp1, Ham_II_t, W0_bath.data, -elementtype(im), zero(elementtype))
    LinearAlgebra.mul!(tmp1, W0_bath.data, Ham_II_t, elementtype(im), one(elementtype))
    K_traced = trace_bath(tmp1, agg, FCProd, aggIndices, vibindices)
    # rho0 = trace_bath(W0, agg, FCProd, aggIndices, vibindices)
    # println(t)
    # println(K_traced)
    # println(K_traced .* rho0.data)

    # println(data)
    LinearAlgebra.mul!(drho.data, one(elementtype), K_traced, -elementtype(im), zero(elementtype))
    data = kernel_integrated .* rho.data
    LinearAlgebra.mul!(drho.data, one(elementtype), data, -elementtype(1), one(elementtype))
    # println(-1im * (K_traced .* rho.data))
    # println(kernel_integrated .* rho.data)
    
    return drho
end

function MemoryKernel(t, s, tmp1, tmp2, tmp3, h, p, Ham_II_t)
    Ham_0, Ham_I, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, Ham_B, W0, W0_bath, agg, FCProd, aggIndices, vibindices = p
    # sprintln(t)
    # Ham_II_s = getInteractionHamIPictureA(Ham_S, Ham_int, s)
    tmp3 .= getInteractionHamIPictureA(Ham_I.data, Ham_0_lambda, Ham_0_S, Ham_0_Sinv, s)
    # LinearAlgebra.mul!(tmp1, H_S, exp(-1im * H_lambda * s), one(elementtype), zero(elementtype))
    # LinearAlgebra.mul!(tmp3, tmp1, H_Sinv, one(elementtype), zero(elementtype))
    
    # LinearAlgebra.mul!(tmp1, adjoint(tmp3), Ham_int.data, one(elementtype), zero(elementtype))
    # LinearAlgebra.mul!(tmp2, tmp1, tmp3, one(elementtype), zero(elementtype))
    # println(tmp2)
    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end
    # println("rho", rho)
    U_0_op = evolutionOperator(Ham_0, s)
    W0_int_s = U_0_op' * W0_bath * U_0_op
    tmp2 .= ad(rho, W0_int_s.data, agg, FCProd, aggIndices, vibindices)
    QuantumOpticsBase.mul!(tmp1, tmp3, tmp2, elementtype(1), zero(elementtype))
    QuantumOpticsBase.mul!(tmp1, tmp2, tmp3, -elementtype(1), one(elementtype))

    QuantumOpticsBase.mul!(tmp2, Ham_II_t, tmp1, elementtype(1), zero(elementtype))
    QuantumOpticsBase.mul!(tmp2, tmp1, Ham_II_t, -elementtype(1), one(elementtype))
    MK_traced = trace_bath(tmp2, agg, FCProd, aggIndices, vibindices)
    # println(MK_traced)
    return MK_traced
end