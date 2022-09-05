using DrWatson
@quickactivate "cautiousBootstrap"

include(srcdir("cfg.jl"))

using Random, Distributions, DataFrames, CSV, Dates, StatsBase
using Plots, StatsPlots, LaTeXStrings

function applyChart(ch, um, thetaHat, m, yprosp; seed = 123)
    # Random.seed!(seed)
    maxrl_i = length(yprosp)
    thetaHatVec = zeros(maxrl_i)
    thetaHatVec[1] = thetaHat
    di = 1
    diVec = Array{Int}(undef, maxrl_i)
    diVec[1] = di
    thetaHatCaut = thetaHat
    thetaHatCautVec = zeros(maxrl_i)
    thetaHatCautVec[1] = thetaHatCaut

    t_alarm = zeros(0)

    valueVec = zeros(maxrl_i)
    valueVec[1] = get_value(ch)
    i = 1
    while i < maxrl_i
        thetaHatCaut = thetaHatVec[i - di + 1]
        y = yprosp[i]
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
            push!(t_alarm, i)
        end
    end

    return (t_alarm = t_alarm, dat = yprosp, chart_values = valueVec, limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um), parameter_updates = thetaHatCautVec, di = diVec)
    # return @NamedTuple{t_alarm::Int64}((t_alarm = i, ))
end

function applyChartGICP(ch, um, yinit, yprosp, thetaHat, Arl0; beta::Union{Bool, Float64} = 0.2, maxrl=1e04, verbose=true, seed=Int(rand(1:1e06)))
    m = length(yinit)
    if isa(um, CautiousLearning)
        Ats0 = get_ATS(um)
        if Ats0 != 0
            # Calculate limit if Ats0 != 0, otherwise use zero-restarting chart
            if verbose println("Calculating limits for target ATS...") end
            sa_ats = saControlLimits(ch, AdaptiveEstimator(), runSimulation, Ats0, thetaHat, Poisson(thetaHat),
                                    m, verbose=false, Amin=0.1, maxiter=1e05,
                                    gamma=0.015, adjusted=true, seed=seed)
            um = CautiousLearning(L = sa_ats[:h], ATS = Ats0)
        else
            if verbose println("ATS = 0, skipping limit calculation.") end
        end
        # Estimate cautious learning limit
    end

    if beta == false
        chart = deepcopy(ch)
    else
        chart = adjust_chart_gicp(ch, um, yinit, thetaHat, runSimulation, m, Arl0, beta=beta, verbose=verbose, seed=seed + 1)
    end

    return applyChart(chart, um, thetaHat, m, yprosp)
end


fold = "ICUadmissions"
dat = DataFrame(CSV.File(datadir(fold, "New_York_Forward_COVID-19_Daily_Hospitalization_Summary_by_Region.csv")))
nydat = filter(row -> row.Region == "NEW YORK CITY", dat)
nydat.date .= Date.(nydat[:, 1], dateformat"mm/dd/yyyy")
nydat = filter(row -> year(row.date) == 2020, nydat)

using Plots.PlotMeasures
y2020 = nydat[:, 4]
days2020 = nydat[:, 5]
ic_2020 = 125:175
oc_2020 = (ic_2020[end]+1):length(y2020)
n_ic = length(ic_2020)
n_oc = length(oc_2020)
pl = plot(days2020, y2020, label="", dpi=400, xrotation=45, bottom_margin=3mm)
plot!(pl, days2020[ic_2020], fill(0, n_ic), fillrange=[fill(maximum(y2020), n_ic)], color=:gray, fillcolor = "gray", fillalpha=0.25, alpha=0.0, label="IC")
safesave(plotsdir(fold, "ICU-cases-2020"),pl)

y = y2020[[ic_2020; oc_2020]]
days = days2020[[ic_2020; oc_2020]]

ic_idx = 1:n_ic
oc_idx = (n_ic+1):(n_ic+n_oc)
yIC = y[ic_idx]
daysIC = days[ic_idx]
yOC = y[oc_idx]
daysOC = days[oc_idx]

println(daysIC[[1, end]])

println(days[111])

pl = plot(days, y, label="", dpi=400, xrotation=45, bottom_margin=3mm)
plot!(pl, days[1:n_ic], fill(0, n_ic), fillrange=[fill(maximum(y), n_ic)], color=:gray, fillcolor = "gray", fillalpha=0.25, alpha=0.0, label="IC", legend=:bottomright)
τ = n_ic + 9
vline!([days[τ]], color=:gray, linestyle=:dot, linewidth=2,  label="", markersize=2.5)
safesave(plotsdir(fold, "ICU-IC-OC.png"), pl)

println(days[τ])

thetaHat = mean(yIC)
Arl0 = 500

Random.seed!(2022-08-18)
D = Poisson
fname = plotsdir(fold, "alarms.jld2")

umVec = [CautiousLearning(ATS=0), FixedParameter(), AdaptiveEstimator()]
nms = ["CL-OC", "FP-OC", "AE-OC"]
for i in eachindex(umVec)
    filesave = datadir(fold, nms[i]*".jld2")
    if isfile(filesave)
        res = load(filesave)["res"]
    else
        ch = signedAEWMA(l=0.019, k=9.7699, L = 1.0)
        beta = 0.1
        um = CautiousLearning(ATS=0)
        res = applyChartGICP(ch, um, yIC, yOC, thetaHat, Arl0, beta=beta)
    end
        plotsave = plotsdir(fold, nms[i]*".png")
        pl = plot(res.chart_values[1:25], label=L"C_t", legend=:outerright, dpi=400)
        hline!([res.limit_alarm], style=:dash, colour="red", xlab=L"t", ylab=L"C_t", label="")
        tau = Int(first(res.t_alarm))
        scatter!([tau], [res.chart_values[tau]], colour="red", label="")
        safesave(filesave, @strdict res)
        safesave(plotsave, pl)
        println(tau,"\t", daysOC[tau])
end
