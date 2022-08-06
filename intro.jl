using DrWatson
@quickactivate "CautiousLearning"

using Revise
using StatisticalProcessControl
using Distributions

include(srcdir("generate_data.jl"))
include(srcdir("update_parameter.jl"))

function chart_statistic(x, thetaHat)
    return (x - thetaHat)/sqrt(thetaHat)
end

function runSimulation(ch, um, thetaHat, D, m; IC=true, tau=1, delta=0.0, maxrl=1e04)
    thetaHatVec = zeros(maxrl)
    thetaHatVec[1] = thetaHat
    di = 1
    diVec = Array{Int}(undef, maxrl)
    diVec[1] = di
    thetaHatCaut = thetaHat
    thetaHatCautVec = zeros(maxrl)
    thetaHatCautVec[1] = thetaHatCaut

    valueVec = zeros(maxrl)
    valueVec[1] = get_value(ch)
    i = 1
    while i < maxrl
        thetaHatCaut = thetaHatVec[i - di + 1]
        y = gen_data_seq(D, i, IC=IC, tau=tau, delta=delta)
        ch = update_series(ch, chart_statistic(y, thetaHatCaut))
        thetaHat = update_parameter(thetaHat, y, i + m)
        i += 1
        valueVec[i] = get_value(ch)
        thetaHatVec[i] = thetaHat
        if check_update(ch, um)
            di = 1
        else
            di += 1
        end
        diVec[i] = di
        thetaHatCautVec[i] = thetaHatCaut
        if check_OC(ch)
            return (t_alarm = i, chart_values = valueVec[1:i], limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um), parameter_updates = thetaHatCautVec[1:i], di = diVec[1:i])
        end

    end

    return (t_alarm = maxrl, chart_values = valueVec, limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um), parameter_updates = thetaHatCautVec, di = diVec)
end

function runExperiment(ch, um, D, m; IC=true, tau=1, delta=0.0, maxrl=1e04)
    yinit = rand(D, m)
    thetaHat = mean(yinit)

    return runSimulation(ch, um, thetaHat, D, m, IC=IC, tau=tau, delta=delta, maxrl=maxrl)
end

ch = AEWMA(l=0.15, L=1.2)
um = CautiousLearning(L=0.2)
# um = SelfStarting()
# um = FixedParameter()
D = Poisson(5.0)
m = 50
maxrl = 3000
IC = false
delta = 0.5
tau = 100

tmp = runExperiment(ch, um, D, m, maxrl=maxrl, IC=IC, delta=delta, tau=tau)

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