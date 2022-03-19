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

@testset "aggregateCore" begin

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1, mode2], 3, [2.0, 200.0])
    mol2 = Molecule([mode1], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50

    @test aggCore.molecules == [mol1, mol2]
    @test aggCore.coupling == [0.0 0.0 0.0; 0.0 0.0 50.0; 0.0 50.0 0.0]
    @test aggCore.molCount == 2

    @test getNvib(aggCore) == [[3, 3], [3]]
    @test getShifts(aggCore) == [[1.0, 2.0], [1.0]]
    @test getFrequencies(aggCore) == [[0.2, 0.3], [0.2]]

    @test getAggStateEnergy(aggCore, [1, 1], [[1, 1], [1]]) == 5.35
    @test getAggStateEnergy(aggCore, [1, 1], [[1, 2], [1]]) == 5.65
    @test getAggStateEnergy(aggCore, [1, 1], [[2, 1], [1]]) == 5.55
    @test getAggStateEnergy(aggCore, [1, 1], [[1, 1], [2]]) == 5.55
    @test getAggStateEnergy(aggCore, [2, 1], [[1, 1], [1]]) == 203.35

    @test OpenQuantumSystems.elIndOrder([1, 1, 1]) == 1
    @test OpenQuantumSystems.elIndOrder([2, 1, 1]) == 2
    @test OpenQuantumSystems.elIndOrder([1, 2, 1]) == 3

    aggCore_ = AggregateCore([mol1, mol2], [0.0 0.0 0.0; 0.0 0.0 50.0; 0.0 50.0 0.0])
    @test aggCore_ == aggCore
end
