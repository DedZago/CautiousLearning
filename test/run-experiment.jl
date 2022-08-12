using DrWatson
@quickactivate "CautiousLearning"
include(srcdir("cfg.jl"))

@testset "Run experiment" begin
    cfg = SimulationSettings(Arl0 = 5, um = CautiousLearning(ATS=0), seed = 2022-08-12)
    runExperiment(cfg, verbose=true)
end