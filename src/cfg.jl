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
    theta = 5.0
    D = Poisson(theta)
    m = 50
    ch = signedEWMA(l=0.2, L=0.8)
    um = CautiousLearning(ATS = 3)
    IC = [true, false]
    delta = [0.35, 0.5, 0.75, 1.0, 1.25, 1.5]
    tau = [1, 50]
    beta = 0.1
    Arl0 = 500
    ncond = 10000
    seed = 1234567
    maxrl = 10000
end

include(srcdir("simulate_runs.jl"))


DrWatson.default_prefix(e::SimulationSettings) = "SimulationSettings"
DrWatson.default_allowed(::SimulationSettings) = (Real, String, Distribution, UnivariateSeries, AdaptiveEstimator, FixedParameter, CautiousLearning)
Base.string(::AdaptiveEstimator) = "AdaptiveEstimator"
Base.string(::FixedParameter) = "FixedParameter"
Base.string(um::CautiousLearning) = "CautiousLearning(ATS="* string(get_ATS(um))*")"
Base.string(ch::UnivariateSeries) = string(typeof(ch)) * string(get_params(ch))


# umVec = [CautiousLearning(ATS=0)]
umVec = [FixedParameter(), AdaptiveEstimator(), CautiousLearning(ATS=0), CautiousLearning(ATS=5)]
folder = "test-single-500"
config = [SimulationSettings(um = um, seed = 2022-08-12) for um in umVec]
