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

@testset "aggregateOperators" begin

    D(op1::Array, op2::Array) = abs(norm(op1 - op2))

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
    
    @test 1e-14 > D(aggOperators.Ham_sys.data, [5.0 0.0 0.0; 0.0 203.5 50.0; 0.0 50.0 306.0])

    ref = [0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.55 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.85 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.45 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.75 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.05 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.65 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.95 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.25]
    @test 1e-14 > D(aggOperators.Ham_bath.data, ref)

    ref = [5.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 5.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0; 0.0 0.0 5.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0; 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0; 0.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0; 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0; 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 203.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 14.325239843009498 20.25894847023148 20.25894847023148 -10.12947423511573 -14.325239843009511 -14.325239843009511 5.064737117557862 7.162619921504751 7.162619921504751; 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 203.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -20.258948470231473 -14.325239843009513 1.6845980341844322e-14 14.325239843009507 10.129474235115739 -1.1911906935453394e-14 -7.162619921504749 -5.064737117557866 5.955953467726693e-15; 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 0.0 203.5 0.0 0.0 0.0 0.0 0.0 0.0 20.258948470231466 -8.890346519833037e-15 -14.325239843009516 -14.325239843009502 6.286424311272163e-15 10.12947423511574 7.162619921504747 -3.1432121556360795e-15 -5.064737117557867; 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 203.5 0.0 0.0 0.0 0.0 0.0 10.12947423511573 14.325239843009513 14.325239843009513 7.162619921504748 10.129474235115739 10.129474235115739 -10.743929882257122 -15.194211352673607 -15.194211352673607; 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 203.5 0.0 0.0 0.0 0.0 -14.325239843009511 -10.12947423511574 1.1911906935453396e-14 -10.129474235115737 -7.1626199215047555 8.422990170922161e-15 15.194211352673603 10.743929882257133 -1.2634485256383239e-14; 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 -0.0 -0.0 0.0 0.0 0.0 0.0 0.0 203.5 0.0 0.0 0.0 14.325239843009504 -6.286424311272164e-15 -10.129474235115742 10.129474235115731 -4.445173259916518e-15 -7.162619921504756 -15.194211352673598 6.667759889874776e-15 10.743929882257133; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 203.5 0.0 0.0 5.064737117557865 7.1626199215047555 7.1626199215047555 10.743929882257127 15.194211352673614 15.194211352673614 1.7906549803761975 2.5323685587789493 2.5323685587789493; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 203.5 0.0 -7.162619921504754 -5.064737117557869 5.955953467726697e-15 -15.194211352673612 -10.743929882257138 1.2634485256383247e-14 -2.532368558778949 -1.790654980376199 2.1057475427305525e-15; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 203.5 7.162619921504751 -3.1432121556360815e-15 -5.06473711755787 15.194211352673607 -6.66775988987478e-15 -10.743929882257138 2.532368558778948 -1.111293314979136e-15 -1.7906549803761995; 0.0 -0.0 0.0 0.0 -0.0 0.0 0.0 -0.0 0.0 14.325239843009498 -20.258948470231473 20.258948470231466 10.12947423511573 -14.325239843009511 14.325239843009504 5.064737117557865 -7.162619921504754 7.162619921504751 306.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 20.25894847023148 -14.325239843009513 -8.890346519833037e-15 14.325239843009513 -10.12947423511574 -6.286424311272164e-15 7.1626199215047555 -5.064737117557869 -3.1432121556360815e-15 0.0 306.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -0.0 0.0 0.0 -0.0 0.0 0.0 -0.0 20.25894847023148 1.6845980341844322e-14 -14.325239843009516 14.325239843009513 1.1911906935453396e-14 -10.129474235115742 7.1626199215047555 5.955953467726697e-15 -5.06473711755787 0.0 0.0 306.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -0.0 0.0 0.0 -0.0 0.0 0.0 -0.0 0.0 -10.12947423511573 14.325239843009507 -14.325239843009502 7.162619921504748 -10.129474235115737 10.129474235115731 10.743929882257127 -15.194211352673612 15.194211352673607 0.0 0.0 0.0 306.0 0.0 0.0 0.0 0.0 0.0; 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 -14.325239843009511 10.129474235115739 6.286424311272163e-15 10.129474235115739 -7.1626199215047555 -4.445173259916518e-15 15.194211352673614 -10.743929882257138 -6.66775988987478e-15 0.0 0.0 0.0 0.0 306.0 0.0 0.0 0.0 0.0; 0.0 0.0 -0.0 0.0 0.0 -0.0 0.0 0.0 -0.0 -14.325239843009511 -1.1911906935453394e-14 10.12947423511574 10.129474235115739 8.422990170922161e-15 -7.162619921504756 15.194211352673614 1.2634485256383247e-14 -10.743929882257138 0.0 0.0 0.0 0.0 0.0 306.0 0.0 0.0 0.0; 0.0 -0.0 0.0 0.0 -0.0 0.0 0.0 -0.0 0.0 5.064737117557862 -7.162619921504749 7.162619921504747 -10.743929882257122 15.194211352673603 -15.194211352673598 1.7906549803761975 -2.532368558778949 2.532368558778948 0.0 0.0 0.0 0.0 0.0 0.0 306.0 0.0 0.0; 0.0 -0.0 -0.0 0.0 -0.0 -0.0 0.0 -0.0 -0.0 7.162619921504751 -5.064737117557866 -3.1432121556360795e-15 -15.194211352673607 10.743929882257133 6.667759889874776e-15 2.5323685587789493 -1.790654980376199 -1.111293314979136e-15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 306.0 0.0; 0.0 0.0 -0.0 0.0 0.0 -0.0 0.0 0.0 -0.0 7.162619921504751 5.955953467726693e-15 -5.064737117557867 -15.194211352673607 -1.2634485256383239e-14 10.743929882257133 2.5323685587789493 2.1057475427305525e-15 -1.7906549803761995 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 306.0]
    @test 1e-14 > D(aggOperators.Ham_S.data, ref)

    ref = [0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.55 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.85 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.45 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.75 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.05 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.65 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.95 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.55 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.85 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.45 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.75 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.05 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.65 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.95 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.55 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.85 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.45 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.75 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.05 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.65 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.95 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.25]
    @test 1e-14 > D(aggOperators.Ham_B.data, ref)

    ref = aggOperators.Ham_S.data + aggOperators.Ham_B.data
    @test 1e-14 > D(aggOperators.Ham_0.data, ref)

    ref = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7071067811865475 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7071067811865475 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7071067811865475 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7071067811865475 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7071067811865475 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7071067811865475 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.414213562373095 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.414213562373095 0.0 -2.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.414213562373095 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.414213562373095 0.0 -2.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.414213562373095 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.414213562373095 0.0 -2.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0 0.0]
    @test 1e-9 > D(aggOperators.Ham_I.data, ref)

    ref = aggOperators.Ham_0.data + aggOperators.Ham_I.data
    @test 1e-14 > D(aggOperators.Ham.data, ref)

end
