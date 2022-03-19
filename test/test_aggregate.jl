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

    # D(op1::Array, op2::Array) = abs(norm(op1 - op2))

    mode1 = Mode(0.2, 1.0)
    mode2 = Mode(0.3, 2.0)
    Energy = [0.0, 200.0]
    mol1 = Molecule([mode1], 3, [2.0, 200.0])
    mol2 = Molecule([mode2], 3, [3.0, 300.0])
    aggCore = AggregateCore([mol1, mol2])
    aggCore.coupling[2, 3] = 50
    aggCore.coupling[3, 2] = 50
    aggTools = AggregateTools(aggCore)
    aggOperators = AggregateOperators(aggCore, aggTools; groundEnergy=true)

    agg = setupAggregate(aggCore; groundEnergy=true)

    @test agg.core == aggCore
    @test agg.tools == aggTools
    @test agg.operators == aggOperators

    agg = Aggregate(aggCore, nothing, nothing)
    agg_ = Aggregate(aggCore)
    @test agg_ == agg

    agg = setupAggregate(aggCore)
    aggTools_ = AggregateTools(agg)
    @test aggTools_ == aggTools

    aggOperators_ = AggregateOperators(agg)
    @test aggOperators_ == aggOperators

    agg = setupAggregate(aggCore)
    agg_ = Aggregate(aggCore, nothing, nothing)
    setupAggregate!(agg_)
    @test agg_ == agg


end
