using DrWatson
@quickactivate "CautiousLearning"

using Revise
using StatisticalProcessControl

function chart_statistic(x, thetaHat)
    return (x - thetaHat)/thetaHat
end

function runSim(ch, um, D, m; IC=true, tau=1, delta=0.0, maxrl=1e04)
    yinit = rand(D, m)
    thetaHat = mean(yinit)
    thetaHatVec = [thetaHat]
    di = 1
    thetaHatCaut = thetaHat
    diVec = [1]
    thetaHatCautVec = [thetaHatCaut]

    valueVec = Vector{Float64}()
    i = 1

    while i < maxrl
        thetaHatCaut = thetaHatVec[i - di + 1]
        y = gen_data_seq(D, i, IC=IC, tau=tau, delta=delta)
        ch = update_series(ch, chart_statistic(y, thetaHatCaut))
        if check_OC(ch)
            return (t_alarm = i, chart_values = valueVec, limit = get_limits(ch), parameter_updates = thetaHatCautVec, di = diVec)
        end

        thetaHat = update_parameter(thetaHat, y, i + m)
        push!(thetaHatVec, thetaHat)

        if check_update(ch, um)
            di = 1
        else
            di += 1
        end
        push!(diVec, di)
        push!(thetaHatCautVec, thetaHatCaut)
        push!(valueVec, get_value(ch))
        i += 1
    end

    return (t_alarm = maxrl, chart_values = valueVec, limit = get_limits(ch), parameter_updates = thetaHatCautVec, di = diVec)
end

tmp = runSim(AEWMA(l = 0.1, L = 0.5), CautiousLearning(0.05),Poisson(20.0), 50, maxrl=5000, IC=false, delta=1.0, tau=100)
plot(tmp[:chart_values])
plot(tmp[:parameter_updates])
plot(tmp[:di])