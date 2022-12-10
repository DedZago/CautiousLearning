using DrWatson
@quickactivate "CautiousLearning"
include(srcdir("cfg.jl"))

@testset "Run experiment" begin
    ch = signedAEWMA(L=0.5)
    cfg = SimulationSettings(ch=ch, ncond=10, Arl0=50, um=CautiousLearning(ATS=0), seed=2022-08-23)
    ncond   = cfg.ncond
    ch      = cfg.ch
    um      = cfg.um
    D       = cfg.D
    m       = cfg.m
    Arl0    = cfg.Arl0
    beta    = cfg.beta
    IC      = cfg.IC
    tau     = cfg.tau
    delta   = cfg.delta
    maxrl   = cfg.maxrl
    simulation = cfg.simulation
        
    # Set simulation seed
    Random.seed!(cfg.seed + simulation)

    yinit = rand(D, m)
    thetaHat = mean(yinit)
    @testset "RunSimulation" begin
        icRun = runSimulation(ch, um, thetaHat, D, m, IC=true, tau=1, delta=0.0, maxrl=maxrl)
        @test isa(icRun, NamedTuple{(:t_alarm,), Tuple{Int64}})
        @test icRun[:t_alarm] == 105
    end
    
    @testset "RunGICP" begin
        out = runGICP(cfg, verbose=true)
    end

    @testset "RunGICP" begin
        res = runConditionalSimulations(out, verbose=true)
    end

    @testset "CautiousLearning C&M" begin
        um = CautiousLearningCM()
        icRun = runSimulation(ch, um, thetaHat, D, m, IC=true, tau=1, delta=0.0, maxrl=maxrl)
        @test isa(icRun, NamedTuple{(:t_alarm,), Tuple{Int64}})
        @test icRun[:t_alarm] == 105
    end
end


um = CautiousLearning(ATS=0)
um2 = AdaptiveEstimator()
cfg = typeof(cfg)(cfg, um=um)
cfg2 = typeof(cfg)(cfg, um=um2)

using BenchmarkTools
# @btime runExperiment(cfg, verbose=false)
# @btime runExperiment(cfg2, verbose=false)
# @btime runSimulation(ch, um2, thetaHat, D, m, IC=true, tau=1, delta=0.0, maxrl=maxrl, seed=seed)
