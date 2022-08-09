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
    ch = AEWMA(l=0.2, L=1.2)
    um = CautiousLearning(ATS = 3)
    IC = [true, false]
    delta = [-1.5, -1.25, -1.0, -0.75, -0.5, 0.5, 0.75, 1.0, 1.25, 1.5]
    tau = [1, 50]
    beta = 0.2
    Arl0 = 200
    ncond = 10000
    seed = 1234567
    maxrl = 10000
    verbose = true
end

include(srcdir("simulate_runs.jl"))


DrWatson.default_prefix(e::SimulationSettings) = "SimulationSettings"
DrWatson.default_allowed(::SimulationSettings) = (Real, String, Distribution, UnivariateSeries, SelfStarting, FixedParameter, CautiousLearning)
Base.string(::SelfStarting) = "SelfStarting"
Base.string(::FixedParameter) = "FixedParameter"
Base.string(um::CautiousLearning) = "CautiousLearning(ATS="* string(get_ATS(um))*")"
Base.string(ch::UnivariateSeries) = string(typeof(ch)) * string(get_params(ch))
