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

@testset "aggregateTools" begin

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1, mode2], 3, [2.0, 200.0])
    mol2 = Molecule([mode1], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    aggTools = AggregateTools(aggCore)
    
    @test aggTools.elIndices == [[1, 1], [2, 1], [1, 2]]
    @test aggTools.vibIndices == [[[1, 1], [1]], [[1, 1], [2]], [[1, 1], [3]], [[1, 2], [1]], [[1, 2], [2]], [[1, 2], [3]], [[1, 3], [1]], [[1, 3], [2]], [[1, 3], [3]], [[2, 1], [1]], [[2, 1], [2]], [[2, 1], [3]], [[2, 2], [1]], [[2, 2], [2]], [[2, 2], [3]], [[2, 3], [1]], [[2, 3], [2]], [[2, 3], [3]], [[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]], [[3, 2], [1]], [[3, 2], [2]], [[3, 2], [3]], [[3, 3], [1]], [[3, 3], [2]], [[3, 3], [3]]]
    @test aggTools.indices == [[[1, 1], [[1, 1], [1]]], [[1, 1], [[1, 1], [2]]], [[1, 1], [[1, 1], [3]]], [[1, 1], [[1, 2], [1]]], [[1, 1], [[1, 2], [2]]], [[1, 1], [[1, 2], [3]]], [[1, 1], [[1, 3], [1]]], [[1, 1], [[1, 3], [2]]], [[1, 1], [[1, 3], [3]]], [[1, 1], [[2, 1], [1]]], [[1, 1], [[2, 1], [2]]], [[1, 1], [[2, 1], [3]]], [[1, 1], [[2, 2], [1]]], [[1, 1], [[2, 2], [2]]], [[1, 1], [[2, 2], [3]]], [[1, 1], [[2, 3], [1]]], [[1, 1], [[2, 3], [2]]], [[1, 1], [[2, 3], [3]]], [[1, 1], [[3, 1], [1]]], [[1, 1], [[3, 1], [2]]], [[1, 1], [[3, 1], [3]]], [[1, 1], [[3, 2], [1]]], [[1, 1], [[3, 2], [2]]], [[1, 1], [[3, 2], [3]]], [[1, 1], [[3, 3], [1]]], [[1, 1], [[3, 3], [2]]], [[1, 1], [[3, 3], [3]]], [[2, 1], [[1, 1], [1]]], [[2, 1], [[1, 1], [2]]], [[2, 1], [[1, 1], [3]]], [[2, 1], [[1, 2], [1]]], [[2, 1], [[1, 2], [2]]], [[2, 1], [[1, 2], [3]]], [[2, 1], [[1, 3], [1]]], [[2, 1], [[1, 3], [2]]], [[2, 1], [[1, 3], [3]]], [[2, 1], [[2, 1], [1]]], [[2, 1], [[2, 1], [2]]], [[2, 1], [[2, 1], [3]]], [[2, 1], [[2, 2], [1]]], [[2, 1], [[2, 2], [2]]], [[2, 1], [[2, 2], [3]]], [[2, 1], [[2, 3], [1]]], [[2, 1], [[2, 3], [2]]], [[2, 1], [[2, 3], [3]]], [[2, 1], [[3, 1], [1]]], [[2, 1], [[3, 1], [2]]], [[2, 1], [[3, 1], [3]]], [[2, 1], [[3, 2], [1]]], [[2, 1], [[3, 2], [2]]], [[2, 1], [[3, 2], [3]]], [[2, 1], [[3, 3], [1]]], [[2, 1], [[3, 3], [2]]], [[2, 1], [[3, 3], [3]]], [[1, 2], [[1, 1], [1]]], [[1, 2], [[1, 1], [2]]], [[1, 2], [[1, 1], [3]]], [[1, 2], [[1, 2], [1]]], [[1, 2], [[1, 2], [2]]], [[1, 2], [[1, 2], [3]]], [[1, 2], [[1, 3], [1]]], [[1, 2], [[1, 3], [2]]], [[1, 2], [[1, 3], [3]]], [[1, 2], [[2, 1], [1]]], [[1, 2], [[2, 1], [2]]], [[1, 2], [[2, 1], [3]]], [[1, 2], [[2, 2], [1]]], [[1, 2], [[2, 2], [2]]], [[1, 2], [[2, 2], [3]]], [[1, 2], [[2, 3], [1]]], [[1, 2], [[2, 3], [2]]], [[1, 2], [[2, 3], [3]]], [[1, 2], [[3, 1], [1]]], [[1, 2], [[3, 1], [2]]], [[1, 2], [[3, 1], [3]]], [[1, 2], [[3, 2], [1]]], [[1, 2], [[3, 2], [2]]], [[1, 2], [[3, 2], [3]]], [[1, 2], [[3, 3], [1]]], [[1, 2], [[3, 3], [2]]], [[1, 2], [[3, 3], [3]]]]
    @test aggTools.indicesMap == [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27], [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54], [55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81]]
    
    mode1 = Mode(0.2, 1.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 2, [2.0, 200.0])
    mol2 = Molecule([mode1], 2, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    aggTools = AggregateTools(aggCore)

    ref = [1.0 0.0 0.0 0.0 0.7788007830714049 0.0 0.5506953149031838 0.0 0.7788007830714049 0.5506953149031838 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.7788007830714049 0.0 0.5506953149031838 -0.5506953149031837 0.3894003915357024 -0.0 0.0; 0.0 0.0 1.0 0.0 -0.5506953149031837 -0.0 0.3894003915357024 0.0 0.0 0.0 0.7788007830714049 0.5506953149031838; 0.0 0.0 0.0 1.0 -0.0 -0.5506953149031837 0.0 0.3894003915357024 -0.0 0.0 -0.5506953149031837 0.3894003915357024; 0.7788007830714049 0.0 -0.5506953149031837 -0.0 1.0 0.0 0.0 0.0 0.6065306597126334 0.42888194248035344 -0.42888194248035333 -0.3032653298563167; 0.0 0.7788007830714049 -0.0 -0.5506953149031837 0.0 1.0 0.0 0.0 -0.42888194248035333 0.30326532985631666 0.30326532985631666 -0.21444097124017664; 0.5506953149031838 0.0 0.3894003915357024 0.0 0.0 0.0 1.0 0.0 0.42888194248035344 0.30326532985631677 0.30326532985631666 0.2144409712401767; 0.0 0.5506953149031838 0.0 0.3894003915357024 0.0 0.0 0.0 1.0 -0.3032653298563167 0.2144409712401767 -0.21444097124017664 0.15163266492815833; 0.7788007830714049 -0.5506953149031837 0.0 -0.0 0.6065306597126334 -0.42888194248035333 0.42888194248035344 -0.3032653298563167 1.0 0.0 0.0 0.0; 0.5506953149031838 0.3894003915357024 0.0 0.0 0.42888194248035344 0.30326532985631666 0.30326532985631677 0.2144409712401767 0.0 1.0 0.0 0.0; 0.0 -0.0 0.7788007830714049 -0.5506953149031837 -0.42888194248035333 0.30326532985631666 0.30326532985631666 -0.21444097124017664 0.0 0.0 1.0 0.0; 0.0 0.0 0.5506953149031838 0.3894003915357024 -0.3032653298563167 -0.21444097124017664 0.2144409712401767 0.15163266492815833 0.0 0.0 0.0 1.0]
    @test 1e-14 > D(aggTools.FCfactors, ref)
    
    @test aggTools.bSystemSize == 3
    @test aggTools.bBathSize == 4
    @test aggTools.bSize == 12

    @test aggTools.basisSystem == GenericBasis([3])
    @test aggTools.basisBath == GenericBasis([4])
    @test aggTools.basis == GenericBasis([12])

end