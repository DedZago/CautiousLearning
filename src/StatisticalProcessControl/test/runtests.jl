using StatisticalProcessControl
using Test


@testset "StatisticalProcessControl.jl" begin
    tests = ["types",
             ]

    for t in tests
        fpath = "$t.jl"
        println("running $fpath ...")
        include(fpath)
    end
end
