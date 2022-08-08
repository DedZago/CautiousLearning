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

ch = EWMA(l=0.15, L=1.2)
um = CautiousLearning(L = 0.12)
um = SelfStarting()
# um = FixedParameter()

IC = [true, false]
delta = [0.5]
tau = 1

Arl0 = 30
ncond = 10000
tmp = runExperiment(1, ncond, ch, um, D, m, Arl0, IC=IC, delta=delta, tau=tau);
mean([out[:t_alarm]] for out in tmp)

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