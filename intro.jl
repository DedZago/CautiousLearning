using DrWatson
@quickactivate "CautiousLearning"

using Revise
using StatisticalProcessControl
using Distributions

include(srcdir("generate_data.jl"))
include(srcdir("update_parameter.jl"))
include(srcdir("simulate_runs.jl"))


ch = AEWMA(l=0.15, L=1.5)
# um = CautiousLearning(L=0.2)
um = SelfStarting()
# um = FixedParameter()
theta = 5.0
D = Poisson(theta)
m = 50
maxrl = 3000

Arl0 = 8

sa_cl = saControlLimits(ch, um, runSimulation, Arl0, theta, D, m, verbose=false, Amin=0.1, maxiter=1e05, z=1.7)

ch = typeof(ch)(ch, L = sa_cl[:h])
mean([runSimulation(ch, um, theta, D, m, maxrl=maxrl)[:t_alarm] for _ in 1:10000])

ch = AEWMA(l=0.15, L=1.2)
um = CautiousLearning(L = sa_cl[:h])
# um = SelfStarting()
# um = FixedParameter()

IC = false
delta = 0.75
tau = 50
tmpList = [runExperiment(ch, um, D, m, IC=IC, delta=delta, tau=tau) for _ in 1:10000]

mean([sim[:t_alarm] for sim in tmpList])


tmp = runExperiment(ch, um, D, m, IC=IC, delta=delta, tau=tau)

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


#! Implement GICP using confidence intervals
# outLow = stochasticApproximationPCC(alpha, runSim, ARL0, model, θlow, direction, Amin=0.001, maxiter=4e04, z = 1.5)
# outUp = stochasticApproximationPCC(alpha, runSim, ARL0, model, θup, direction, Amin=0.001, maxiter=4e04, z = 1.5)

# limits = (outLow[1], outUp[1])
# L = minimum(limits)

# tmp = [runSim(L, model, θ, maxrl = 1e04, direction=direction) for _ in 1:1e04]
# mean(tmp)

# B = 100
# Lvec = zeros(B)
# for b in 1:B
#     println("b: ", b)
#     θsim = rand(get_posterior(model))
#     out = stochasticApproximationPCC(alpha, runSim, ARL0, model, θsim, direction, Amin=0.001, maxiter=4e04, z = 1.5, verbose=false)
#     Lvec[b] = out[1]
# end
# L = quantile(Lvec, beta)
# density(Lvec)
