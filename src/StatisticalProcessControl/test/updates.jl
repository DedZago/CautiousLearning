using StatisticalProcessControl

@testset "Parameter updates" begin
    L = 3.0
    C = 1.5
    @testset "Double sided" begin
        @testset "AdaptiveEstimator" begin
            UM = AdaptiveEstimator()
            CH = EWMA(L=L, C=C)
            @test check_update(CH, UM) == true

            CH = AEWMA(L=L, C=C)
            @test check_update(CH, UM) == true
        end

        @testset "FixedParameter" begin
            UM = FixedParameter()
            CH = EWMA(L=L, C=C)
            @test check_update(CH, UM) == false

            CH = AEWMA(L=L, C=C)
            @test check_update(CH, UM) == false
        end

        @testset "Cautious Learning" begin
            UM = CautiousLearning(L=1.6)
            CH = EWMA(L=L, C=C)
            @test check_update(CH, UM) == true
            CH = AEWMA(L=L, C=C)
            @test check_update(CH, UM) == true

            UM = CautiousLearning(L=1.1)
            CH = EWMA(L=L, C=C)
            @test check_update(CH, UM) == false
            CH = AEWMA(L=L, C=C)
            @test check_update(CH, UM) == false


            Random.seed!(1)
            CH = AEWMA(L=L)
            x = randn(1)
            CH = update_series(CH, first(x))
            UM = CautiousLearning(L=1.1)
            @test check_update(CH, UM) == true
            y = 10.0
            CH = update_series(CH, first(y))
            @test check_update(CH, UM) == false
            @test check_IC(CH) == false
            @test check_OC(CH) == true
        end
    end

    @testset "One sided" begin
        @testset "AdaptiveEstimator" begin
            UM = AdaptiveEstimator()
            CH = signedEWMA(L=L, C=C)
            @test check_update(CH, UM) == true

            CH = signedAEWMA(L=L, C=C)
            @test check_update(CH, UM) == true
        end

        @testset "FixedParameter" begin
            UM = FixedParameter()
            CH = signedEWMA(L=L, C=C)
            @test check_update(CH, UM) == false

            CH = signedAEWMA(L=L, C=C)
            @test check_update(CH, UM) == false
        end

        @testset "Cautious Learning" begin
            UM = CautiousLearning(L=1.6)
            CH = signedEWMA(L=L, C=C)
            @test check_update(CH, UM) == true
            CH = signedAEWMA(L=L, C=C)
            @test check_update(CH, UM) == true

            UM = CautiousLearning(L=1.1)
            CH = signedEWMA(L=L, C=C)
            @test check_update(CH, UM) == false
            CH = signedAEWMA(L=L, C=C)
            @test check_update(CH, UM) == false


            Random.seed!(1)
            CH = signedAEWMA(L=L)
            x = randn(1)
            CH = update_series(CH, first(x))
            UM = CautiousLearning(L=1.1)
            @test check_update(CH, UM) == true
            y = 10.0
            CH = update_series(CH, first(y))
            @test check_update(CH, UM) == false
            @test check_IC(CH) == false
            @test check_OC(CH) == true
        end

        @testset "Restarting" begin
            UM = CautiousLearning(L=1.6, ATS = 0)
            CH = signedEWMA(L=L, C=C)
            @test check_update(CH, UM) == false
            CH = signedAEWMA(L=L, C=C)
            @test check_update(CH, UM) == false

            UM = CautiousLearning(ATS = 0)
            CH = signedEWMA(L=L, C=0.0)
            @test check_update(CH, UM) == true
            CH = signedAEWMA(L=L, C=0.0)
            @test check_update(CH, UM) == true


            Random.seed!(1)
            CH = signedAEWMA(L=L)
            x = randn(1000)
            zeros_chart = []
            updates_UM = []
            UM = CautiousLearning(ATS = 0)
            for i in eachindex(x)
                CH = update_series(CH, x[i])
                if get_value(CH) == 0.0
                    push!(zeros_chart, i)
                end
                if check_update(CH, UM)
                    push!(updates_UM, i)
                end
            end
            @test zeros_chart == updates_UM
        end
    end
end