
include(srcdir("update_parameter.jl"))


@testset "Poisson parameter" begin
    Random.seed!(1)
    nsim = 100
    x = randn(nsim)
    xbar = mean(x)
    xrec = 0.0
    D = Poisson(0.5)
    for i in 1:nsim
        xrec = update_parameter(xrec, x[i], i, D)
    end
    @test isapprox(xrec, xbar)
end