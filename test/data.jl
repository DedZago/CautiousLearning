include(srcdir("generate_data.jl"))


@testset "Poisson generation" begin
    theta = 2.0
    delta = 2.0
    D = Poisson(theta)

    delta_st = standardize_shift(delta, D)
    @test delta_st == delta*sqrt(theta)
    @test get_shifted_distribution(D, delta) == Poisson(theta + delta_st)

    Random.seed!(1)
    n = 10
    nic = 2
    noc = 8
    x = Float64.([rand(D, nic); rand(Poisson(theta + delta_st), noc)])

    Random.seed!(1)
    y = gen_ic_oc_data(n, D, IC=false, tau=nic+1, delta=delta)

    Random.seed!(1)
    z = [gen_data_seq(D, i, IC=false, tau=nic+1, delta=delta) for i in 1:n]

    @test x == y
    @test y == z


    Random.seed!(2)
    theta = 30
    D = Poisson(theta)
    confintVec = [get_confint(rand(D, 100)) for _ in 1:1e05]
    @test isapprox(mean([confintVec[i][1] <= theta <= confintVec[i][2] for i in eachindex(confintVec)]), 0.95, atol=0.02)
end