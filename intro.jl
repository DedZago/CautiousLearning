using Distributed
addprocs(Sys.CPU_THREADS - 1)
@everywhere using DrWatson
@everywhere @quickactivate "CautiousLearning"

@everywhere using Revise
@everywhere using Distributions
@everywhere using Parameters
@everywhere using SharedArrays
@everywhere using Random
@everywhere using DataFrames
@everywhere using StatisticalProcessControl

@everywhere include(srcdir("generate_data.jl"))
@everywhere include(srcdir("update_parameter.jl"))

@everywhere @with_kw struct SimulationSettings
    theta = 5.0
    D = Poisson(theta)
    m = 50
    ch = AEWMA(l=0.2, L=1.2)
    um = CautiousLearning(ATS = 3)
    IC = [true, false]
    delta = [-1.5, -1.25, -1.0, -0.75, -0.5, 0.5, 0.75, 1.0, 1.25, 1.5]
    tau = [1, 50]
    beta = 0.2
    Arl0 = 200
    nsim = 200
    ncond = 10000
    seed = 11223344556677
    maxrl = 10000
    verbose = true
end

@everywhere include(srcdir("simulate_runs.jl"))


DrWatson.default_prefix(e::SimulationSettings) = "SimulationSettings"
DrWatson.default_allowed(::SimulationSettings) = (Real, String, Distribution, UnivariateSeries, SelfStarting, FixedParameter, CautiousLearning)
Base.string(::SelfStarting) = "SelfStarting"
Base.string(::FixedParameter) = "FixedParameter"
Base.string(um::CautiousLearning) = "CautiousLearning(ATS="* string(get_ATS(um))*")"
Base.string(ch::UnivariateSeries) = string(typeof(ch)) * string(get_params(ch))


umVec = [FixedParameter(seed = 2022-08-09)]
# umVec = [FixedParameter(), SelfStarting(), CautiousLearning(ATS=3)]
config = [SimulationSettings(um = um) for um in umVec]
for cfg in config
    svname = datadir("sims", "test", savename(cfg, "jld2"))
    if isfile(svname)
        continue
    else
        out = runExperiment(cfg)
        sv = @strdict config=cfg out
        safesave(datadir("sims", "test", savename(cfg, "jld2")), sv)
    end
end