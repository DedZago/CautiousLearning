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

function runSim(ch, um, D, m; IC=true, tau=1, delta=0.0, maxrl=1e04)
    yinit = rand(D, m)
    thetaHat = mean(yinit)
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
        valueVec[i+1] = get_value(ch)
        if check_OC(ch)
            return (t_alarm = i, chart_values = valueVec[1:i], limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um),parameter_updates = thetaHatCautVec[1:i], di = diVec[1:i])
        end

        thetaHat = update_parameter(thetaHat, y, i + m)
        thetaHatVec[i+1] = thetaHat

        if check_update(ch, um)
            di = 1
        else
            di += 1
        end
        diVec[i+1] = di
        thetaHatCautVec[i+1] = thetaHatCaut
        i += 1
    end

    return (t_alarm = maxrl, chart_values = valueVec, limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um), parameter_updates = thetaHatCautVec, di = diVec)
end

tmp = runSim(EWMA(l = 0.2, L = 1.5), CautiousLearning(0.3),Poisson(5.0), 50, maxrl=500, IC=false, delta=0.5, tau=100)

using Plots

plot(tmp[:chart_values], legend=:outerright)
hline!([-tmp[:limit_alarm],tmp[:limit_alarm]], fill=true, fillcolour=:orange, fillalpha=0.15, color=:darkred, style=:dash)
hline!([-tmp[:limit_cautious],tmp[:limit_cautious]], fill=true, fillcolour=:lightgreen, fillalpha=0.35, color=:darkgreen, style=:dash)
plot(tmp[:parameter_updates])
plot(tmp[:di])