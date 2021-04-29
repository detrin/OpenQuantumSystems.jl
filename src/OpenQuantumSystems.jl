module OpenQuantumSystems

using SparseArrays, LinearAlgebra
import LinearAlgebra: mul!, rmul!

export bases,
    Basis,
    GenericBasis,
    CompositeBasis,
    basis,
    tensor,
    ⊗,
    permutesystems,
    @samebases,

    # states.jl
    states,
    StateVector,
    Bra,
    Ket,
    basisstate,
    norm,
    dagger,
    normalize,
    normalize!,

    # operators.jl
    operators,
    AbstractOperator,
    DataOperator,
    expect,
    variance,
    identityoperator,
    ptrace,
    embed,
    dense,
    tr,
    sparse,

    # operators_dense.jl
    operators_dense,
    Operator,
    DenseOperator,
    DenseOpType,
    projector,
    dm,

    # operators_sparse.jl
    operators_sparse,
    SparseOperator,
    diagonaloperator,
    SparseOpType,

    # superoperators.jl
    superoperators,
    SuperOperator,
    DenseSuperOperator,
    DenseSuperOpType,
    SparseSuperOperator,
    SparseSuperOpType,
    spre,
    spost,
    liouvillian,

    # fock.jl
    fock,
    FockBasis,
    number,
    destroy,
    create,
    fockstate,
    coherentstate,
    coherentstate!,
    displace,

    # spin.jl
    spin,
    SpinBasis,
    sigmax,
    sigmay,
    sigmaz,
    sigmap,
    sigmam,
    spinup,
    spindown,

    # state_definitions.jl
    randstate,
    randoperator,
    thermalstate,
    coherentthermalstate,
    phase_average,
    passive_state,

    # nlevel.jl
    nlevel,
    NLevelBasis,
    transition,
    nlevelstate,

    # metrics.jl
    metrics,
    tracenorm,
    tracenorm_h,
    tracenorm_nh,
    tracedistance,
    tracedistance_h,
    tracedistance_nh,
    entropy_vn,
    entropy_renyi,
    fidelity,
    ptranspose,
    PPT,
    negativity,
    logarithmic_negativity,
    entanglement_entropy,
    PauliBasis,
    PauliTransferMatrix,
    DensePauliTransferMatrix,
    ChiMatrix,
    DenseChiMatrix,
    avg_gate_fidelity,

    # molecules.jl
    Mode,
    franckCondonFactors,
    vibrationalIndices,
    electronicIndices,
    Molecule,
    getMolStateEnergy,
    getMolStateFC,

    # aggregate.jl
    Aggregate,
    getNvib,
    vibrationalIndices,
    electronicIndices,
    getIndices,
    getFranckCondonFactors,
    getAggStateEnergy,
    getAggHamiltonian,
    getFranckCondonFactorsSparse,
    getAggHamiltonianSparse

include("core.jl")
include("sortedindices.jl")
include("bases.jl")
include("states.jl")
include("operators.jl")
include("operators_dense.jl")
include("sparsematrix.jl")
include("operators_sparse.jl")
include("superoperators.jl")
include("spin.jl")
include("fock.jl")
include("state_definitions.jl")
include("pauli.jl")
include("metrics.jl")
include("nlevel.jl")
include("molecules.jl")
include("aggregate.jl")

end
