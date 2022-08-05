using StatisticalProcessControl

@testset "Parameter updates" begin
    L = 3.0
    C = 1.5
    @testset "SelfStarting" begin
        UM = SelfStarting()
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