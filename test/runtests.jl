using DrWatson
@quickactivate "CautiousLearning"

using Test


@testset "CautiousLearning" begin
    tests = ["data",
                "parameter-updates"
             ]

    for t in tests
        fpath = "$t.jl"
        println("running $fpath ...")
        include(fpath)
    end
end
