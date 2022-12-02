using Test
using OpenQuantumSystems
using Random, SparseArrays, LinearAlgebra, StableRNGs

@testset "rate_constant" begin

    Random.seed!(StableRNG(0), 1)

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))
    D(x1::StateVector, x2::StateVector) = norm(x2 - x1)
    D(op1::AbstractOperator, op2::AbstractOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))
    D(op1::AbstractSuperOperator, op2::AbstractSuperOperator) =
        abs(tracedistance_nh(dense(op1), dense(op2)))


    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 3, [2.0, 200.0])
    mol2 = Molecule([mode2], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    agg = setupAggregate(aggCore)
    aggTools = agg.tools
    aggOperators = agg.operators

    Ham_B = agg.operators.Ham_B
    Ham_I = agg.operators.Ham_I
    Ham_0 = agg.operators.Ham_0
    Ham = agg.operators.Ham

    Ham_0_lambda, Ham_0_S = eigen(Ham_0.data)
    Ham_0_Sinv = inv(Ham_0_S)
    Ham_0_lambda = diagm(Ham_0_lambda)

    basis = agg.tools.basis
    indicesLen = agg.tools.bSize
    indices = agg.tools.indices
    indicesMap = agg.tools.indicesMap
    FCFact = agg.tools.FCfactors
    FCProd = agg.tools.FCproduct

    tspan = get_tspan(0., 0.02, 100)
    W0, rho0, W0_bath = ultrafast_laser_excitation(10., [0.0, 0.3, 0.7], agg)

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))

    t = 0.01
    

end

function W_abcd_1_bath_core(t, t1, t2, p, tmp1, tmp2; bath_evolution=:none)
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, _ = p
    
    # W_1_bath = deepcopy(W0_bath.data)
    # elLen = aggCore.molCount+1
    # indicesMap = aggTools.indicesMap
    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I

    Ham_II_t1 = getInteractionHamIPicture(Ham_0, Ham_I, t1)
    Ham_II_t2 = getInteractionHamIPicture(Ham_0, Ham_I, t2)
    
    W_1_bath = zeros(ComplexF64, aggTools.bSize, aggTools.bSize)
    rho_t = OpenQuantumSystems.interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    rho_t2 = OpenQuantumSystems.interpolate_with_tspan(rho_0_int_t_itp, tspan, t2)
    W_bath_t2 = OpenQuantumSystems.interpolate_with_tspan(W_0_bath_t_itp, tspan, t2)

    tmp1[:, :] = ad(rho_t2, W_bath_t2, aggCore, aggTools)
    # tmp1[:, :] = W_bath_t2[:, :]
    tmp2[:, :] = Ham_II_t2.data * tmp1 - tmp1 * Ham_II_t2.data
    W_1_bath[:, :] = Ham_II_t1.data * tmp2 - tmp2 * Ham_II_t1.data
    
    return W_1_bath
end

function W_1_bath(
        t, p, tmp1, tmp2; 
        bath_evolution=:none, bath_ansatz=:population, normalize=false,
        W_1_rtol=1e-12, W_1_atol=1e-12
    )
    aggCore, aggTools, aggOperators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, _ = p
    elLen = aggCore.molCount+1
    indicesMap = aggTools.indicesMap
    
    W_1_diff, err = QuadGK.quadgk(
        t1 -> begin
            W_1_diff_t1, err = QuadGK.quadgk(
                t2 -> W_abcd_1_bath_core(t, t1, t2, p, tmp1, tmp2; bath_evolution=bath_evolution),
                0,
                t1,
                rtol = W_1_rtol,
                atol = W_1_atol,
            );
            return W_1_diff_t1
        end,
        0,
        t,
        rtol = W_1_rtol,
        atol = W_1_atol,
    )
    rho_t = OpenQuantumSystems.interpolate_with_tspan(rho_0_int_t_itp, tspan, t)
    W_bath_t = OpenQuantumSystems.interpolate_with_tspan(W_0_bath_t_itp, tspan, t)
    for a=2:elLen, b=1:elLen
        a1 = indicesMap[a][1]
        a2 = indicesMap[a][end]
        b1 = indicesMap[b][1]
        b2 = indicesMap[b][end]
        if rho_t[a, b] == 0 
            tmp1[a1:a2, b1:b2] = W_1_diff[a1:a2, b1:b2] 
        else
            tmp1[a1:a2, b1:b2] = W_1_diff[a1:a2, b1:b2]  / rho_t[a, b]
        end
    end
    W_1_bath = deepcopy(W_bath_t) - tmp1
    if normalize
        W_1_bath = normalize_bath(W_1_bath, aggCore, aggTools, aggOperators)
    end
    return W_1_bath
end

function W_1_bath(t, W0, W0_bath, agg::Aggregate; W_1_rtol=1e-12, W_1_atol=1e-12, K_rtol=1e-12, K_atol=1e-12)
    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, eltype(W0))
    return W_1_bath(t, p, tmp1, tmp2; W_1_rtol=W_1_rtol, W_1_atol=W_1_atol, K_rtol=K_rtol, K_atol=K_atol)
end

function QME_sI_iterative(
    W0::T,
    rho_0_int_t,
    W_0_bath_t,
    tspan::Array,
    agg::Aggregate;
    bath_evolution=:none, 
    bath_ansatz=:population, 
    normalize=false,
    reltol::AbstractFloat = 1.0e-12,
    abstol::AbstractFloat = 1.0e-12,
    int_reltol::AbstractFloat = 1.0e-4,
    int_abstol::AbstractFloat = 1.0e-4,
    W_1_rtol::AbstractFloat = 1e-12, 
    W_1_atol::AbstractFloat = 1e-12, 
    alg::Any = DelayDiffEq.MethodOfSteps(DelayDiffEq.Vern6()),
    fout::Union{Function,Nothing} = nothing,
    kwargs...,
) where {B<:Basis,T<:Operator{B,B}}
    history_fun(p, t) = T(rho0.basis_l, rho0.basis_r, zeros(ComplexF64, size(rho0.data)))
    rho0 = trace_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)
    W0_bath = get_rho_bath(W0, agg.core, agg.tools; vib_basis=agg.operators.vib_basis)

    tmp1 = copy(W0.data)
    tmp2 = copy(W0.data)
    
    # Calculate and interpolate rho_0_int_t, W_0_bath_t, W_1_bath_t    
    if ndims(rho_0_int_t) == 1 && typeof(rho_0_int_t[1]) <: Operator
        rho_0_int_t = operator_recast(rho_0_int_t)
    end
    if ndims(rho_0_int_t) == 3
        rho_0_int_t = [rho_0_int_t[i, :, :] for i=1:size(rho_0_int_t, 1)]
    end
    rho_0_int_t_itp = Interpolations.interpolate(
        rho_0_int_t,
        Interpolations.BSpline(Interpolations.Linear())
    )
    
    if ndims(W_0_bath_t) == 1 && typeof(W_0_bath_t[1]) <: Operator
        W_0_bath_t = operator_recast(W_0_bath_t)
    end
    if ndims(W_0_bath_t) == 3
        W_0_bath_t = [W_0_bath_t[i, :, :] for i=1:size(W_0_bath_t, 1)]
    end
    W_0_bath_t_itp = Interpolations.interpolate(
        W_0_bath_t, 
        Interpolations.BSpline(Interpolations.Linear())
    )
    println("running W_1_bath")
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, tspan, rho_0_int_t_itp, W_0_bath_t_itp, eltype(W0))
    elLen = agg.core.molCount+1
    W_1_bath_t = []
    for t_i=1:length(tspan)
        t = tspan[t_i]
        W_1_bath_ = W_1_bath(
            t, p, tmp1, tmp2; 
            bath_evolution=bath_evolution, bath_ansatz=bath_ansatz, normalize=normalize,
            W_1_rtol=W_1_rtol, W_1_atol=W_1_atol
        )
        push!(W_1_bath_t, W_1_bath_)
    end
    W_1_bath_itp = Interpolations.interpolate(W_1_bath_t, Interpolations.BSpline(Interpolations.Linear()))
    
    println("running dQME_sI_iterative")
    p = (agg.core, agg.tools, agg.operators, W0, W0_bath, rho_0_int_t_itp, W_1_bath_itp, tspan, eltype(W0))
    dmaster_(t, rho, drho, history_fun, p) = dQME_sI_iterative(
        t,
        rho,
        drho,
        history_fun,
        tmp1,
        tmp2,
        p,
        int_reltol,
        int_abstol,
    )
    tspan_ = convert(Vector{float(eltype(tspan))}, tspan)
    x0 = rho0.data
    state = T(rho0.basis_l, rho0.basis_r, rho0.data)
    dstate = T(rho0.basis_l, rho0.basis_r, rho0.data)
    tspan, rho_int_1_t = OpenQuantumSystems.integrate_delayed(
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
    return tspan, rho_int_1_t, W_1_bath_t
end

function dQME_sI_iterative(
    t::AbstractFloat,
    rho::T,
    drho::T,
    history_fun,
    tmp1::Array,
    tmp2::Array,
    p,
    int_reltol::AbstractFloat,
    int_abstol::AbstractFloat,
) where {B<:Basis,T<:Operator{B,B}}
    aggCore, aggTools, aggOperators, W0, _, _, _, _, elementtype = p
        
    Ham_II_t = getInteractionHamIPicture(aggOperators.Ham_0, aggOperators.Ham_I, t)
    K = Ham_II_t.data * W0.data - W0.data * Ham_II_t.data
    K_traced = trace_bath(K, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
    # println("dQME_sI_iterative", " ", t)
    kernel_integrated_traced, err = QuadGK.quadgk(
        s -> kernel_sI_iterative(t, s, history_fun, p, tmp1, tmp2, Ham_II_t),
        0,
        t,
        rtol = int_reltol,
        atol = int_abstol,
    )    
    drho.data[:, :] = -elementtype(im) * K_traced - kernel_integrated_traced

    return drho
end

function kernel_sI_iterative(t, s, h, p, tmp1, tmp2, Ham_II_t)
    aggCore, aggTools, aggOperators, W0, W0_bath, rho_0_int_t_itp, W_1_bath_itp, tspan, _ = p

    rho = h(p, s)

    if (typeof(rho) <: Operator)
        rho = rho.data
    end

    Ham_0 = aggOperators.Ham_0
    Ham_I = aggOperators.Ham_I
    Ham_II_s = getInteractionHamIPicture(Ham_0, Ham_I, s)
    W_bath_s = OpenQuantumSystems.interpolate_with_tspan(W_1_bath_itp, tspan, s)

    tmp1[:, :] = ad(rho, W_bath_s, aggCore, aggTools)
    tmp2[:, :] = Ham_II_s.data * tmp1 - tmp1 * Ham_II_s.data
    tmp1[:, :] = Ham_II_t.data * tmp2 - tmp2 * Ham_II_t.data

    return trace_bath(tmp1, aggCore, aggTools; vib_basis=aggOperators.vib_basis)
end