using Test
using OpenQuantumSystems
using LinearAlgebra
using SparseArrays

# https://github.com/JuliaIO/Suppressor.jl/blob/master/src/Suppressor.jl
macro suppress(block)
    quote
        if ccall(:jl_generating_output, Cint, ()) == 0
            original_stdout = stdout
            out_rd, out_wr = redirect_stdout()
            out_reader = @async read(out_rd, String)

            original_stderr = stderr
            err_rd, err_wr = redirect_stderr()
            err_reader = @async read(err_rd, String)

            # approach adapted from https://github.com/JuliaLang/IJulia.jl/pull/667/files
            logstate = Base.CoreLogging._global_logstate
            logger = logstate.logger
            if logger.stream == original_stderr
                new_logstate =
                    Base.CoreLogging.LogState(typeof(logger)(err_wr, logger.min_level))
                Core.eval(Base.CoreLogging, Expr(:(=), :(_global_logstate), new_logstate))
            end
        end

        try
            $(esc(block))
        finally
            if ccall(:jl_generating_output, Cint, ()) == 0
                redirect_stdout(original_stdout)
                close(out_wr)

                redirect_stderr(original_stderr)
                close(err_wr)

                if logger.stream == stderr
                    Core.eval(Base.CoreLogging, Expr(:(=), :(_global_logstate), logstate))
                end
            end
        end
    end
end

@testset "aggregate" begin

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))


    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.4, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1, mode2], 3, Energy)
    mol2 = Molecule([mode1, mode2], 3, Energy)
    agg = Aggregate([mol1, mol2])

    @test getNvib(agg) == [[3, 3], [3, 3]]
    @test electronicIndices(agg) == [[1, 1], [2, 1], [1, 2]]

    vibInds = [
        [[1, 1], [1, 1]],
        [[1, 1], [1, 2]],
        [[1, 1], [1, 3]],
        [[1, 1], [2, 1]],
        [[1, 1], [2, 2]],
        [[1, 1], [2, 3]],
        [[1, 1], [3, 1]],
        [[1, 1], [3, 2]],
        [[1, 1], [3, 3]],
        [[1, 2], [1, 1]],
        [[1, 2], [1, 2]],
        [[1, 2], [1, 3]],
        [[1, 2], [2, 1]],
        [[1, 2], [2, 2]],
        [[1, 2], [2, 3]],
        [[1, 2], [3, 1]],
        [[1, 2], [3, 2]],
        [[1, 2], [3, 3]],
        [[1, 3], [1, 1]],
        [[1, 3], [1, 2]],
        [[1, 3], [1, 3]],
        [[1, 3], [2, 1]],
        [[1, 3], [2, 2]],
        [[1, 3], [2, 3]],
        [[1, 3], [3, 1]],
        [[1, 3], [3, 2]],
        [[1, 3], [3, 3]],
        [[2, 1], [1, 1]],
        [[2, 1], [1, 2]],
        [[2, 1], [1, 3]],
        [[2, 1], [2, 1]],
        [[2, 1], [2, 2]],
        [[2, 1], [2, 3]],
        [[2, 1], [3, 1]],
        [[2, 1], [3, 2]],
        [[2, 1], [3, 3]],
        [[2, 2], [1, 1]],
        [[2, 2], [1, 2]],
        [[2, 2], [1, 3]],
        [[2, 2], [2, 1]],
        [[2, 2], [2, 2]],
        [[2, 2], [2, 3]],
        [[2, 2], [3, 1]],
        [[2, 2], [3, 2]],
        [[2, 2], [3, 3]],
        [[2, 3], [1, 1]],
        [[2, 3], [1, 2]],
        [[2, 3], [1, 3]],
        [[2, 3], [2, 1]],
        [[2, 3], [2, 2]],
        [[2, 3], [2, 3]],
        [[2, 3], [3, 1]],
        [[2, 3], [3, 2]],
        [[2, 3], [3, 3]],
        [[3, 1], [1, 1]],
        [[3, 1], [1, 2]],
        [[3, 1], [1, 3]],
        [[3, 1], [2, 1]],
        [[3, 1], [2, 2]],
        [[3, 1], [2, 3]],
        [[3, 1], [3, 1]],
        [[3, 1], [3, 2]],
        [[3, 1], [3, 3]],
        [[3, 2], [1, 1]],
        [[3, 2], [1, 2]],
        [[3, 2], [1, 3]],
        [[3, 2], [2, 1]],
        [[3, 2], [2, 2]],
        [[3, 2], [2, 3]],
        [[3, 2], [3, 1]],
        [[3, 2], [3, 2]],
        [[3, 2], [3, 3]],
        [[3, 3], [1, 1]],
        [[3, 3], [1, 2]],
        [[3, 3], [1, 3]],
        [[3, 3], [2, 1]],
        [[3, 3], [2, 2]],
        [[3, 3], [2, 3]],
        [[3, 3], [3, 1]],
        [[3, 3], [3, 2]],
        [[3, 3], [3, 3]],
    ]
    @test vibrationalIndices(agg) == vibInds

    aggInds = [
        [[1, 1], [[1, 1], [1, 1]]],
        [[1, 1], [[1, 1], [1, 2]]],
        [[1, 1], [[1, 1], [1, 3]]],
        [[1, 1], [[1, 1], [2, 1]]],
        [[1, 1], [[1, 1], [2, 2]]],
        [[1, 1], [[1, 1], [2, 3]]],
        [[1, 1], [[1, 1], [3, 1]]],
        [[1, 1], [[1, 1], [3, 2]]],
        [[1, 1], [[1, 1], [3, 3]]],
        [[1, 1], [[1, 2], [1, 1]]],
        [[1, 1], [[1, 2], [1, 2]]],
        [[1, 1], [[1, 2], [1, 3]]],
        [[1, 1], [[1, 2], [2, 1]]],
        [[1, 1], [[1, 2], [2, 2]]],
        [[1, 1], [[1, 2], [2, 3]]],
        [[1, 1], [[1, 2], [3, 1]]],
        [[1, 1], [[1, 2], [3, 2]]],
        [[1, 1], [[1, 2], [3, 3]]],
        [[1, 1], [[1, 3], [1, 1]]],
        [[1, 1], [[1, 3], [1, 2]]],
        [[1, 1], [[1, 3], [1, 3]]],
        [[1, 1], [[1, 3], [2, 1]]],
        [[1, 1], [[1, 3], [2, 2]]],
        [[1, 1], [[1, 3], [2, 3]]],
        [[1, 1], [[1, 3], [3, 1]]],
        [[1, 1], [[1, 3], [3, 2]]],
        [[1, 1], [[1, 3], [3, 3]]],
        [[1, 1], [[2, 1], [1, 1]]],
        [[1, 1], [[2, 1], [1, 2]]],
        [[1, 1], [[2, 1], [1, 3]]],
        [[1, 1], [[2, 1], [2, 1]]],
        [[1, 1], [[2, 1], [2, 2]]],
        [[1, 1], [[2, 1], [2, 3]]],
        [[1, 1], [[2, 1], [3, 1]]],
        [[1, 1], [[2, 1], [3, 2]]],
        [[1, 1], [[2, 1], [3, 3]]],
        [[1, 1], [[2, 2], [1, 1]]],
        [[1, 1], [[2, 2], [1, 2]]],
        [[1, 1], [[2, 2], [1, 3]]],
        [[1, 1], [[2, 2], [2, 1]]],
        [[1, 1], [[2, 2], [2, 2]]],
        [[1, 1], [[2, 2], [2, 3]]],
        [[1, 1], [[2, 2], [3, 1]]],
        [[1, 1], [[2, 2], [3, 2]]],
        [[1, 1], [[2, 2], [3, 3]]],
        [[1, 1], [[2, 3], [1, 1]]],
        [[1, 1], [[2, 3], [1, 2]]],
        [[1, 1], [[2, 3], [1, 3]]],
        [[1, 1], [[2, 3], [2, 1]]],
        [[1, 1], [[2, 3], [2, 2]]],
        [[1, 1], [[2, 3], [2, 3]]],
        [[1, 1], [[2, 3], [3, 1]]],
        [[1, 1], [[2, 3], [3, 2]]],
        [[1, 1], [[2, 3], [3, 3]]],
        [[1, 1], [[3, 1], [1, 1]]],
        [[1, 1], [[3, 1], [1, 2]]],
        [[1, 1], [[3, 1], [1, 3]]],
        [[1, 1], [[3, 1], [2, 1]]],
        [[1, 1], [[3, 1], [2, 2]]],
        [[1, 1], [[3, 1], [2, 3]]],
        [[1, 1], [[3, 1], [3, 1]]],
        [[1, 1], [[3, 1], [3, 2]]],
        [[1, 1], [[3, 1], [3, 3]]],
        [[1, 1], [[3, 2], [1, 1]]],
        [[1, 1], [[3, 2], [1, 2]]],
        [[1, 1], [[3, 2], [1, 3]]],
        [[1, 1], [[3, 2], [2, 1]]],
        [[1, 1], [[3, 2], [2, 2]]],
        [[1, 1], [[3, 2], [2, 3]]],
        [[1, 1], [[3, 2], [3, 1]]],
        [[1, 1], [[3, 2], [3, 2]]],
        [[1, 1], [[3, 2], [3, 3]]],
        [[1, 1], [[3, 3], [1, 1]]],
        [[1, 1], [[3, 3], [1, 2]]],
        [[1, 1], [[3, 3], [1, 3]]],
        [[1, 1], [[3, 3], [2, 1]]],
        [[1, 1], [[3, 3], [2, 2]]],
        [[1, 1], [[3, 3], [2, 3]]],
        [[1, 1], [[3, 3], [3, 1]]],
        [[1, 1], [[3, 3], [3, 2]]],
        [[1, 1], [[3, 3], [3, 3]]],
        [[2, 1], [[1, 1], [1, 1]]],
        [[2, 1], [[1, 1], [1, 2]]],
        [[2, 1], [[1, 1], [1, 3]]],
        [[2, 1], [[1, 1], [2, 1]]],
        [[2, 1], [[1, 1], [2, 2]]],
        [[2, 1], [[1, 1], [2, 3]]],
        [[2, 1], [[1, 1], [3, 1]]],
        [[2, 1], [[1, 1], [3, 2]]],
        [[2, 1], [[1, 1], [3, 3]]],
        [[2, 1], [[1, 2], [1, 1]]],
        [[2, 1], [[1, 2], [1, 2]]],
        [[2, 1], [[1, 2], [1, 3]]],
        [[2, 1], [[1, 2], [2, 1]]],
        [[2, 1], [[1, 2], [2, 2]]],
        [[2, 1], [[1, 2], [2, 3]]],
        [[2, 1], [[1, 2], [3, 1]]],
        [[2, 1], [[1, 2], [3, 2]]],
        [[2, 1], [[1, 2], [3, 3]]],
        [[2, 1], [[1, 3], [1, 1]]],
        [[2, 1], [[1, 3], [1, 2]]],
        [[2, 1], [[1, 3], [1, 3]]],
        [[2, 1], [[1, 3], [2, 1]]],
        [[2, 1], [[1, 3], [2, 2]]],
        [[2, 1], [[1, 3], [2, 3]]],
        [[2, 1], [[1, 3], [3, 1]]],
        [[2, 1], [[1, 3], [3, 2]]],
        [[2, 1], [[1, 3], [3, 3]]],
        [[2, 1], [[2, 1], [1, 1]]],
        [[2, 1], [[2, 1], [1, 2]]],
        [[2, 1], [[2, 1], [1, 3]]],
        [[2, 1], [[2, 1], [2, 1]]],
        [[2, 1], [[2, 1], [2, 2]]],
        [[2, 1], [[2, 1], [2, 3]]],
        [[2, 1], [[2, 1], [3, 1]]],
        [[2, 1], [[2, 1], [3, 2]]],
        [[2, 1], [[2, 1], [3, 3]]],
        [[2, 1], [[2, 2], [1, 1]]],
        [[2, 1], [[2, 2], [1, 2]]],
        [[2, 1], [[2, 2], [1, 3]]],
        [[2, 1], [[2, 2], [2, 1]]],
        [[2, 1], [[2, 2], [2, 2]]],
        [[2, 1], [[2, 2], [2, 3]]],
        [[2, 1], [[2, 2], [3, 1]]],
        [[2, 1], [[2, 2], [3, 2]]],
        [[2, 1], [[2, 2], [3, 3]]],
        [[2, 1], [[2, 3], [1, 1]]],
        [[2, 1], [[2, 3], [1, 2]]],
        [[2, 1], [[2, 3], [1, 3]]],
        [[2, 1], [[2, 3], [2, 1]]],
        [[2, 1], [[2, 3], [2, 2]]],
        [[2, 1], [[2, 3], [2, 3]]],
        [[2, 1], [[2, 3], [3, 1]]],
        [[2, 1], [[2, 3], [3, 2]]],
        [[2, 1], [[2, 3], [3, 3]]],
        [[2, 1], [[3, 1], [1, 1]]],
        [[2, 1], [[3, 1], [1, 2]]],
        [[2, 1], [[3, 1], [1, 3]]],
        [[2, 1], [[3, 1], [2, 1]]],
        [[2, 1], [[3, 1], [2, 2]]],
        [[2, 1], [[3, 1], [2, 3]]],
        [[2, 1], [[3, 1], [3, 1]]],
        [[2, 1], [[3, 1], [3, 2]]],
        [[2, 1], [[3, 1], [3, 3]]],
        [[2, 1], [[3, 2], [1, 1]]],
        [[2, 1], [[3, 2], [1, 2]]],
        [[2, 1], [[3, 2], [1, 3]]],
        [[2, 1], [[3, 2], [2, 1]]],
        [[2, 1], [[3, 2], [2, 2]]],
        [[2, 1], [[3, 2], [2, 3]]],
        [[2, 1], [[3, 2], [3, 1]]],
        [[2, 1], [[3, 2], [3, 2]]],
        [[2, 1], [[3, 2], [3, 3]]],
        [[2, 1], [[3, 3], [1, 1]]],
        [[2, 1], [[3, 3], [1, 2]]],
        [[2, 1], [[3, 3], [1, 3]]],
        [[2, 1], [[3, 3], [2, 1]]],
        [[2, 1], [[3, 3], [2, 2]]],
        [[2, 1], [[3, 3], [2, 3]]],
        [[2, 1], [[3, 3], [3, 1]]],
        [[2, 1], [[3, 3], [3, 2]]],
        [[2, 1], [[3, 3], [3, 3]]],
        [[1, 2], [[1, 1], [1, 1]]],
        [[1, 2], [[1, 1], [1, 2]]],
        [[1, 2], [[1, 1], [1, 3]]],
        [[1, 2], [[1, 1], [2, 1]]],
        [[1, 2], [[1, 1], [2, 2]]],
        [[1, 2], [[1, 1], [2, 3]]],
        [[1, 2], [[1, 1], [3, 1]]],
        [[1, 2], [[1, 1], [3, 2]]],
        [[1, 2], [[1, 1], [3, 3]]],
        [[1, 2], [[1, 2], [1, 1]]],
        [[1, 2], [[1, 2], [1, 2]]],
        [[1, 2], [[1, 2], [1, 3]]],
        [[1, 2], [[1, 2], [2, 1]]],
        [[1, 2], [[1, 2], [2, 2]]],
        [[1, 2], [[1, 2], [2, 3]]],
        [[1, 2], [[1, 2], [3, 1]]],
        [[1, 2], [[1, 2], [3, 2]]],
        [[1, 2], [[1, 2], [3, 3]]],
        [[1, 2], [[1, 3], [1, 1]]],
        [[1, 2], [[1, 3], [1, 2]]],
        [[1, 2], [[1, 3], [1, 3]]],
        [[1, 2], [[1, 3], [2, 1]]],
        [[1, 2], [[1, 3], [2, 2]]],
        [[1, 2], [[1, 3], [2, 3]]],
        [[1, 2], [[1, 3], [3, 1]]],
        [[1, 2], [[1, 3], [3, 2]]],
        [[1, 2], [[1, 3], [3, 3]]],
        [[1, 2], [[2, 1], [1, 1]]],
        [[1, 2], [[2, 1], [1, 2]]],
        [[1, 2], [[2, 1], [1, 3]]],
        [[1, 2], [[2, 1], [2, 1]]],
        [[1, 2], [[2, 1], [2, 2]]],
        [[1, 2], [[2, 1], [2, 3]]],
        [[1, 2], [[2, 1], [3, 1]]],
        [[1, 2], [[2, 1], [3, 2]]],
        [[1, 2], [[2, 1], [3, 3]]],
        [[1, 2], [[2, 2], [1, 1]]],
        [[1, 2], [[2, 2], [1, 2]]],
        [[1, 2], [[2, 2], [1, 3]]],
        [[1, 2], [[2, 2], [2, 1]]],
        [[1, 2], [[2, 2], [2, 2]]],
        [[1, 2], [[2, 2], [2, 3]]],
        [[1, 2], [[2, 2], [3, 1]]],
        [[1, 2], [[2, 2], [3, 2]]],
        [[1, 2], [[2, 2], [3, 3]]],
        [[1, 2], [[2, 3], [1, 1]]],
        [[1, 2], [[2, 3], [1, 2]]],
        [[1, 2], [[2, 3], [1, 3]]],
        [[1, 2], [[2, 3], [2, 1]]],
        [[1, 2], [[2, 3], [2, 2]]],
        [[1, 2], [[2, 3], [2, 3]]],
        [[1, 2], [[2, 3], [3, 1]]],
        [[1, 2], [[2, 3], [3, 2]]],
        [[1, 2], [[2, 3], [3, 3]]],
        [[1, 2], [[3, 1], [1, 1]]],
        [[1, 2], [[3, 1], [1, 2]]],
        [[1, 2], [[3, 1], [1, 3]]],
        [[1, 2], [[3, 1], [2, 1]]],
        [[1, 2], [[3, 1], [2, 2]]],
        [[1, 2], [[3, 1], [2, 3]]],
        [[1, 2], [[3, 1], [3, 1]]],
        [[1, 2], [[3, 1], [3, 2]]],
        [[1, 2], [[3, 1], [3, 3]]],
        [[1, 2], [[3, 2], [1, 1]]],
        [[1, 2], [[3, 2], [1, 2]]],
        [[1, 2], [[3, 2], [1, 3]]],
        [[1, 2], [[3, 2], [2, 1]]],
        [[1, 2], [[3, 2], [2, 2]]],
        [[1, 2], [[3, 2], [2, 3]]],
        [[1, 2], [[3, 2], [3, 1]]],
        [[1, 2], [[3, 2], [3, 2]]],
        [[1, 2], [[3, 2], [3, 3]]],
        [[1, 2], [[3, 3], [1, 1]]],
        [[1, 2], [[3, 3], [1, 2]]],
        [[1, 2], [[3, 3], [1, 3]]],
        [[1, 2], [[3, 3], [2, 1]]],
        [[1, 2], [[3, 3], [2, 2]]],
        [[1, 2], [[3, 3], [2, 3]]],
        [[1, 2], [[3, 3], [3, 1]]],
        [[1, 2], [[3, 3], [3, 2]]],
        [[1, 2], [[3, 3], [3, 3]]],
    ]
    @test getIndices(agg) == aggInds

    mode1 = Mode(0.2, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 2, Energy)
    mol2 = Molecule([mode1], 2, Energy)
    agg = Aggregate([mol1, mol2])
    aggInds = getIndices(agg)
    FC_part = getFranckCondonFactors(agg; groundState = false)
    FC1 = getFranckCondonFactors(agg)
    FC2 = getFranckCondonFactors(agg, aggInds)

    @test size(FC1) == (12, 12)
    @test size(FC_part) == (8, 8)
    @test FC1 == FC2
    @test FC_part == FC1[5:12, 5:12]
    FC = [
        1.0 0.0 0.0 0.0 0.6065306597126334 0.42888194248035344 -0.42888194248035333 -0.3032653298563167
        0.0 1.0 0.0 0.0 -0.42888194248035333 0.30326532985631666 0.30326532985631666 -0.21444097124017664
        0.0 0.0 1.0 0.0 0.42888194248035344 0.30326532985631677 0.30326532985631666 0.2144409712401767
        0.0 0.0 0.0 1.0 -0.3032653298563167 0.2144409712401767 -0.21444097124017664 0.15163266492815833
        0.6065306597126334 -0.42888194248035333 0.42888194248035344 -0.3032653298563167 1.0 0.0 0.0 0.0
        0.42888194248035344 0.30326532985631666 0.30326532985631677 0.2144409712401767 0.0 1.0 0.0 0.0
        -0.42888194248035333 0.30326532985631666 0.30326532985631666 -0.21444097124017664 0.0 0.0 1.0 0.0
        -0.3032653298563167 -0.21444097124017664 0.2144409712401767 0.15163266492815833 0.0 0.0 0.0 1.0
    ]
    @test 1e-12 > D(FC_part, FC)

    @test getAggStateEnergy(agg, [1, 1], [[1], [2]]) == 0.4

    @test getAggStateEnergy(agg, [1, 1], [[1], [2]]) == 0.4

    @test OpenQuantumSystems.elIndOrder([1, 1, 1]) == 1
    @test OpenQuantumSystems.elIndOrder([1, 1, 2]) == 4
    @test OpenQuantumSystems.elIndOrder([2, 1, 1]) == 2

    aggInds = getIndices(agg; groundState = false)
    agg.coupling[2, 3] = 200
    agg.coupling[3, 2] = 200
    Ham_ref = [
        0.0 0.0 0.0 0.0 121.30613194252669 85.7763884960707 -85.77638849607067 -60.653065971263345
        0.0 0.20000000000001705 0.0 0.0 -85.77638849607067 60.65306597126333 60.65306597126333 -42.888194248035326
        0.0 0.0 0.20000000000001705 0.0 85.7763884960707 60.65306597126335 60.65306597126333 42.88819424803534
        0.0 0.0 0.0 0.4000000000000341 -60.653065971263345 42.88819424803534 -42.888194248035326 30.326532985631665
        121.30613194252669 -85.77638849607067 85.7763884960707 -60.653065971263345 0.0 0.0 0.0 0.0
        85.7763884960707 60.65306597126333 60.65306597126335 42.88819424803534 0.0 0.20000000000001705 0.0 0.0
        -85.77638849607067 60.65306597126333 60.65306597126333 -42.888194248035326 0.0 0.0 0.20000000000001705 0.0
        -60.653065971263345 -42.888194248035326 42.88819424803534 30.326532985631665 0.0 0.0 0.0 0.4000000000000341
    ]
    Ham1 = getAggHamiltonian(agg, aggInds, FC_part; groundState = false)
    Ham2 = getAggHamiltonian(agg, aggInds; groundState = false)
    Ham3 = getAggHamiltonian(agg; groundState = false)

    @test 1e-12 > D(Ham_ref, Ham1.data)
    @test 1e-12 > D(Ham_ref, Ham2.data)
    @test 1e-12 > D(Ham_ref, Ham3.data)

    Ham_ref2 = [
        200.2 0.0 0.0 0.0 121.30613194252669 85.7763884960707 -85.77638849607067 -60.653065971263345
        0.0 200.4 0.0 0.0 -85.77638849607067 60.65306597126333 60.65306597126333 -42.888194248035326
        0.0 0.0 200.4 0.0 85.7763884960707 60.65306597126335 60.65306597126333 42.88819424803534
        0.0 0.0 0.0 200.60000000000002 -60.653065971263345 42.88819424803534 -42.888194248035326 30.326532985631665
        121.30613194252669 -85.77638849607067 85.7763884960707 -60.653065971263345 200.2 0.0 0.0 0.0
        85.7763884960707 60.65306597126333 60.65306597126335 42.88819424803534 0.0 200.4 0.0 0.0
        -85.77638849607067 60.65306597126333 60.65306597126333 -42.888194248035326 0.0 0.0 200.4 0.0
        -60.653065971263345 -42.888194248035326 42.88819424803534 30.326532985631665 0.0 0.0 0.0 200.60000000000002
    ]
    Ham4 = getAggHamiltonian(agg, aggInds, FC_part; groundEnergy = true)
    @test 1e-12 > D(Ham_ref2, Ham4.data)

    FCSparse = getFranckCondonFactorsSparse(agg, aggInds; groundState = false)
    @test 1e-12 > D(FC, Matrix(FCSparse))

    FCSparse = getFranckCondonFactorsSparse(agg; groundState = false)
    @test 1e-12 > D(FC, Matrix(FCSparse))

    HamSparse1 = getAggHamiltonianSparse(agg, aggInds, FCSparse; groundState = false)
    HamSparse2 = getAggHamiltonianSparse(agg, aggInds; groundState = false)
    HamSparse3 = getAggHamiltonianSparse(agg; groundState = false)
    @test 1e-12 > D(Ham_ref, Matrix(HamSparse1.data))
    @test 1e-12 > D(Ham_ref, Matrix(HamSparse2.data))
    @test 1e-12 > D(Ham_ref, Matrix(HamSparse3.data))


    agg = Aggregate([mol1, mol2])
    agg.coupling[2, 3] = 200
    agg.coupling[3, 2] = 200
    agg2 = Aggregate([mol1, mol2], [0.0 0.0 0.0; 0.0 0.0 200.0; 0.0 200.0 0.0])
    @test agg2.coupling == agg.coupling

    modes = [Mode(2.0, 2.0), Mode(2.0, 2.0)]
    mols = [
        Molecule(modes, 2, [0.0, 200.0]),
        Molecule(modes, 2, [0.0, 300.0]),
        Molecule(modes, 2, [0.0, 400.0]),
    ]

    agg = Aggregate(mols)
    agg.coupling[2, 3] = 100
    agg.coupling[3, 2] = 100
    agg.coupling[3, 4] = 100
    agg.coupling[4, 3] = 100

    Ham_sys = getAggHamiltonianSystem(agg)
    Ham_sys_ref = [200.0 100.0 0.0; 100.0 300.0 100.0; 0.0 100.0 400.0]
    @test 1e12 > D(Ham_sys.data, Ham_sys_ref)

    Ham_sys = getAggHamiltonianSystem(agg; groundEnergy = false)
    Ham_sys_ref = [0.0 100.0 0.0; 100.0 100.0 100.0; 0.0 100.0 200.0]
    @test 1e12 > D(Ham_sys.data, Ham_sys_ref)

    Ham_sys_ref =
        [0.0 0.0 0.0 0.0; 0.0 200.0 100.0 0.0; 0.0 100.0 300.0 100.0; 0.0 0.0 100.0 400.0]
    Ham_sys = getAggHamiltonianSystem(agg; groundState = true)
    @test 1e12 > D(Ham_sys.data, Ham_sys_ref)

    Ham_sys = getAggHamiltonianSystem(agg; groundState = true, groundEnergy = false)
    Ham_sys_ref = [0.0 0.0 0.0 0.0; 0.0 200.0 100.0 0.0; 0.0 100.0 300.0 100.0; 0.0 0.0 100.0 400.0]
    @test 1e12 > D(Ham_sys.data, Ham_sys_ref)

    modes = [Mode(2.0, 2.0)]
    mols = [Molecule(modes, 2, [0.0, 200.0]), Molecule(modes, 2, [0.0, 300.0])]

    agg = Aggregate(mols)
    agg.coupling[2, 3] = 100
    agg.coupling[3, 2] = 100

    Ham_bath_ref = [2.0 0.0 0.0 0.0; 0.0 4.0 0.0 0.0; 0.0 0.0 4.0 0.0; 0.0 0.0 0.0 6.0]
    Ham_bath = getAggHamiltonianBath(agg)
    @test 1e12 > D(Ham_bath.data, Ham_bath_ref)

    Ham_bath = getAggHamiltonianBath(agg)
    Ham_sys = getAggHamiltonianSystem(agg)
    b_sys = GenericBasis([size(Ham_sys, 1)])
    b_bath = GenericBasis([size(Ham_bath, 1)])

    Ham_S =
        tensor(OneDenseOperator(b_bath), Ham_sys) +
        tensor(Ham_bath, OneDenseOperator(b_sys))
    Ham_int = getAggHamiltonianInteraction(agg)
    @test 1e12 > D(Ham_ref, Ham_int.data + Ham_S.data)

    aggInds = getIndices(agg)
    vibIndices = getVibIndices(agg, aggInds)
    @test vibIndices == [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]

    aggInds = getIndices(agg; groundState = false)
    FCfact = getFranckCondonFactors(agg; groundState = false)
    Ham_S2 = getAggHamSysBath(agg, aggInds; groundState = false, groundEnergy = true)
    @test 1e12 > D(Ham_S.data, Ham_S2.data)

    Ham_S3 = getAggHamSysBath2(agg, aggInds; groundState = false, groundEnergy = true)
    @test 1e12 > D(Ham_S.data, Ham_S3.data)


    aggInds_ref = getIndices(agg; groundState = false)
    vibindices_ref = getVibIndices(agg, aggInds_ref)
    aggIndLen_ref = length(aggInds_ref)
    basis_ref = GenericBasis([aggIndLen_ref])
    FCFact_ref = getFranckCondonFactors(agg, aggInds_ref; groundState = false)
    FCProd_ref =
        getFCProd(agg, FCFact_ref, aggInds_ref, vibindices_ref; groundState = false)
    Ham_ref = getAggHamiltonian(
        agg,
        aggInds_ref,
        FCFact_ref;
        groundState = false,
        groundEnergy = true,
    )
    Ham_0_ref = getAggHamSysBath(agg, aggInds_ref; groundState = false, groundEnergy = true)
    Ham_I_ref = Ham_ref - Ham_0_ref

    aggInds, vibindices, aggIndLen, basis, FCFact, FCProd, Ham, Ham_0, Ham_I =
        setupAggregate(agg; groundState = false, groundEnergy = true, verbose = false)
    @test aggInds_ref == aggInds
    @test vibindices_ref == vibindices
    @test aggIndLen_ref == aggIndLen
    @test basis_ref == basis
    @test FCFact_ref == FCFact
    @test FCProd_ref == FCProd
    @test Ham_ref == Ham
    @test Ham_0_ref == Ham_0
    @test Ham_I_ref == Ham_I

    @suppress aggInds, vibindices, aggIndLen, basis, FCFact, FCProd, Ham, Ham_0, Ham_I =
        setupAggregate(agg; groundState = false, groundEnergy = true, verbose = true)
    @test aggInds_ref == aggInds
    @test vibindices_ref == vibindices
    @test aggIndLen_ref == aggIndLen
    @test basis_ref == basis
    @test FCFact_ref == FCFact
    @test FCProd_ref == FCProd
    @test Ham_ref == Ham
    @test Ham_0_ref == Ham_0
    @test Ham_I_ref == Ham_I
end
