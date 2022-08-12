using StatisticalProcessControl
using Random

@testset "Control series" begin
    Random.seed!(1)
    x = randn(100)
    C = 0.0
    l = 0.1
    k = 10.0
    L = 3.0
    @testset "Double sided" begin
        ch = EWMA(l=l, L=L, C=C)
        ch2 = AEWMA(l=l, k=k, C=C, L=L)

        @test isa(ch, EWMA)
        @test isa(ch2, AEWMA)

        pm = get_params(ch)
        pm2 = get_params(ch2)

        @test length(pm) == 2
        @test length(pm2) == 3
        @test pm.l == l
        @test pm.L == L
        @test pm.l == pm2.l
        @test pm.L == pm2.L

        x_sm = apply_series(ch, x)
        x_sm2 = apply_series(ch2, x)
        @test length(x) == length(x_sm)
        @test length(x) == length(x_sm2)
        @test isapprox(x_sm, x_sm2)

        ch = update_series(ch, x[1])
        ch2 = update_series(ch2, x[1])

        @test isapprox(get_value(ch), get_value(ch2))

        # Instantiate OC chart
        chOC = EWMA(l = l, L = L, C = L + 1.0)
        @test check_IC(chOC) == false
        @test check_OC(chOC) == true
    end

    @testset "One sided" begin
        ch = signedEWMA(l=l, L=L, C=C)
        ch2 = signedAEWMA(l=l, k=k, C=C, L=L)

        @test isa(ch, signedEWMA)
        @test isa(ch2, signedAEWMA)

        pm = get_params(ch)
        pm2 = get_params(ch2)

        @test length(pm) == 3
        @test length(pm2) == 4
        @test pm.l == l
        @test pm.L == L
        @test pm.l == pm2.l
        @test pm.L == pm2.L
        @test pm.upw == true
        @test pm.upw == pm2.upw

        x_sm = apply_series(ch, x)
        x_sm2 = apply_series(ch2, x)
        @test length(x) == length(x_sm)
        @test length(x) == length(x_sm2)
        @test isapprox(x_sm, x_sm2)

        ch = update_series(ch, x[1])
        ch2 = update_series(ch2, x[1])

        @test isapprox(get_value(ch), get_value(ch2))

        # Instantiate OC chart
        chOC = signedEWMA(l = l, L = L, C = L + 1.0)
        @test check_IC(chOC) == false
        @test check_OC(chOC) == true
    end
end