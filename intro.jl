using Distributed
addprocs(2)
@everywhere using DrWatson
@everywhere @quickactivate "CautiousLearning"

@everywhere using Revise
@everywhere using StatisticalProcessControl
@everywhere using Distributions
@everywhere using SharedArrays

@everywhere include(srcdir("generate_data.jl"))
@everywhere include(srcdir("update_parameter.jl"))
@everywhere include(srcdir("simulate_runs.jl"))

ch = AEWMA(l=0.15, L=0.6)
# um = CautiousLearning(L=0.2)
um = SelfStarting()
# um = FixedParameter()
theta = 5.0
D = Poisson(theta)
m = 50
maxrl = 3000

Ats0 = 5

sa_ats = saControlLimits(ch, um, runSimulation, Ats0, theta, D, m, verbose=false, Amin=0.1, maxiter=1e05, gamma=0.03, adjusted=true)

ch = typeof(ch)(ch, L = sa_ats[:h])
mean([runSimulation(ch, um, theta, D, m, maxrl=maxrl)[:t_alarm] for _ in 1:10000])

ch = EWMA(l=0.2, L=1.2)
um = CautiousLearning(L = 0.12)
um = SelfStarting()
# um = FixedParameter()

IC = [true, false]
delta = [0.5, 0.25]
tau = [1, 50]
beta = 0.2

Arl0 = 100
nsim = 2
ncond = 10000
out = runNestedSimulations(nsim, ncond, ch, um, D, m, Arl0, beta=beta, IC=IC, delta=delta, tau=tau);
#! Take out data generation from runNestedSimulations in order to compute guaranteed cautious learning limits if needed

mean(out["IC"])
mean(out["OC"][1][1])
mean(out["OC"][2][1])
mean(out["OC"][1][2])
mean(out["OC"][2][2])

using Plots
plot(tmp[:chart_values], legend=:outerright, label="C")
hline!([-tmp[:limit_alarm],tmp[:limit_alarm]], fill=true, fillcolour=:orange, fillalpha=0.15, color=:darkred, style=:dash, label=false)
hline!([-tmp[:limit_cautious],tmp[:limit_cautious]], fill=true, fillcolour=:lightgreen, fillalpha=0.35, color=:darkgreen, style=:dash, label=false)
if !IC
    scatter!([tmp[:t_alarm]], [tmp[:chart_values][end]], label=false, color=:red)    
    vline!([tau], style=:dash, colour=:darkblue, label="tau")
end

plot(tmp[:parameter_updates])
if !IC
    vline!([tau], style=:dash, colour=:darkblue, label="tau")
end
plot(tmp[:di])
if !IC
    vline!([tau], style=:dash, colour=:darkblue, label="tau")
end