using DrWatson
@quickactivate "CautiousLearning"

using Revise
using Distributions
using Parameters
using SharedArrays
using Random
using DataFrames
using StatisticalProcessControl

include(srcdir("generate_data.jl"))
include(srcdir("update_parameter.jl"))

@with_kw struct SimulationSettings
    simulation = 1
    D = Poisson(1.0)
    m = 50
    ch = signedEWMA(l=0.2, L=0.8)
    um = CautiousLearning(ATS = 3)
    IC = [true, false]
    delta = [0.25, 0.35, 0.5, 0.75, 1.0, 1.25, 1.5]
    tau = [1, 50, 100]
    beta = 0.1
    Arl0 = 500
    ncond = 10000
    seed = 1234567
    maxrl = 10000
end

include(srcdir("simulate_runs.jl"))


# DrWatson._wsave(s, fig::T) where T <: Plots.Plot{Plots.GRBackend} = savefig(fig, s)
DrWatson.default_prefix(e::SimulationSettings) = "SimulationSettings"
DrWatson.default_allowed(::SimulationSettings) = (Real, String, Distribution, UnivariateSeries, AdaptiveEstimator, FixedParameter, CautiousLearning)
Base.string(::AdaptiveEstimator) = "AdaptiveEstimator"
Base.string(::FixedParameter) = "FixedParameter"
Base.string(um::CautiousLearning) = "CautiousLearning(ATS="* string(get_ATS(um))*")"
Base.string(um::CautiousLearningCM) = "CautiousLearningCM"
Base.string(ch::UnivariateSeries) = string(typeof(ch)) * string(get_params(ch))


# umVec = [FixedParameter(), AdaptiveEstimator(), CautiousLearning(ATS=0)]
umVec = [CautiousLearningCM()]

# seed = 2022-08-23
# Random.seed!(seed)
# maxseed = Int(1e09)
# config = [
#           [SimulationSettings(ch = signedEWMA(l=0.033, L=1.0), D=Poisson(1.0), um = um, seed = rand(1:maxseed)) for um in umVec];
#             [SimulationSettings(ch = signedEWMA(l=0.0230, L=1.0), D=Poisson(4.0), um = um, seed = rand(1:maxseed)) for um in umVec];
#             [SimulationSettings(ch = signedEWMA(l=0.0190, L=1.0), D=Poisson(7.0), um = um, seed = rand(1:maxseed)) for um in umVec];
#            ]

seed = 2022-09-22
Random.seed!(seed)
maxseed = Int(1e09)
lambda = [0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200]
seed_values = [rand(1:maxseed) for _ in lambda]
config = [SimulationSettings(ch = signedEWMA(l=lambda[i], L=1.0), D=Poisson(4.0), um = um, seed = seed_values[i])  for i in 1:length(lambda) for um in umVec]
