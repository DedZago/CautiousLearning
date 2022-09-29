using DrWatson
@quickactivate "CautiousLearning"

include(srcdir("cfg.jl"))
using Plots, StatsPlots, LaTeXStrings

function runSimulationInfo(ch, um, thetaHat, D, m; IC=true, tau=1, delta=0.0, maxrl=Int(1e04))
    maxrl_i = Int(round(maxrl))
    thetaHatVec = zeros(maxrl_i)
    thetaHatVec[1] = thetaHat
    di = 1
    diVec = Array{Int}(undef, maxrl_i)
    diVec[1] = di
    thetaHatCaut = thetaHat
    thetaHatCautVec = zeros(maxrl_i)
    thetaHatCautVec[1] = thetaHatCaut

    valueVec = zeros(maxrl_i)
    valueVec[1] = get_value(ch)
    i = 1
    while i < maxrl_i
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
            break
        end

    end

    if isa(um, CautiousLearning) && um.ATS == 0.0
        limit_cautious = 0.0
    else
        limit_cautious = get_warning_limit(um)
    end
    return @NamedTuple{t_alarm::Int64, chart_values::Vector{Float64}, limit_alarm::Float64, limit_cautious::Float64, parameter_updates::Vector{Float64}, di::Vector{Int64}}((t_alarm = i, chart_values = valueVec[1:i], limit_alarm = get_limits(ch), limit_cautious = limit_cautious, parameter_updates = thetaHatCautVec[1:i], di = diVec[1:i]))
    # return @NamedTuple{t_alarm::Int64}((t_alarm = i, ))

end

function plot_comparison(res; tau = 0, fill=true)
    tauHat = res[:t_alarm]
    t = 1:tauHat
    R = res[:chart_values]
    S = res[:limit_cautious]
    L = res[:limit_alarm]
    p = plot(t, R,  dpi=400, color="black", xlab=L"t", ylab=L"R_t", label="",legend=:outerright)
    if fill
        # Colour update regions
        if !isinf(S)
            #! For increases only
            hline!(p, [S], label=L"\mathcal{I}_t", fillrange=0.0, fillalpha=0.25, fillcolor="lightgreen", alpha=0.0)
            if S == 0.0
                hline!(p, [S], label="", colour="lightgreen", linewidth=2)
            end
        end
        if S != L
            hline!(p, [L], label=L"\mathcal{S}_t", fillrange=[S], fillalpha=0.15, alpha=0.0, color="orange")
        end
    end
    hline!(p, [L], label="", color = "darkred", linestyle=:dot, linewidth=2)
    scatter!(p, [tauHat], [R[tauHat]], color = :red, label="")
    vline!(p, [tau], linestyle=:dash, color=:blue, label=L"\tau")
    return p
end

function plot_window_of_opportunity(res; jmpsize = 0.1, tau = 0, xlab=L"t", ylab=L"\hat{\theta}_t")
    tauHat = res[:t_alarm]
    t = collect(1:tauHat-1)
    idx = diff(res[:parameter_updates]) .>= jmpsize
    idx_jump = first(t[idx])
    p = plot(res[:parameter_updates], dpi=400, xlab=xlab, ylab=ylab, label="",legend=:bottomright)
    if tau != 0
        vline!(p, [tau], style=:dash, label=L"\tau")
        window_of_opp = collect(tau:idx_jump)
        miny, maxy = ylims(p)
        plot!(p, window_of_opp, fill(miny, length(window_of_opp)), fillrange=fill(maxy, length(window_of_opp)), fillcolor = "gray", fillalpha=0.25, alpha=0.0, label="")
    end
    return p
end


ch = signedEWMA(L=0.5)
cfg = SimulationSettings(ch=ch, D=Poisson(4.0), Arl0=500, um=CautiousLearning(ATS=0), seed=2022-08-23)
ncond   = cfg.ncond
ch      = cfg.ch
um      = cfg.um
D       = cfg.D
m       = cfg.m
Arl0    = cfg.Arl0
beta    = cfg.beta
IC      = cfg.IC
tau     = cfg.tau
delta   = cfg.delta
maxrl   = cfg.maxrl
simulation = cfg.simulation

DrWatson._wsave(s, fig::T) where T <: Plots.Plot{Plots.GRBackend} = savefig(fig, s)
    
# Set simulation seed
Random.seed!(123)
yinit = rand(D, m)
thetaHat = mean(yinit)
L_gicp = runGICP(cfg, verbose=true)

chart = typeof(ch)(ch, L=L_gicp["res"][:L])                        


folder = plotsdir("sims", "window-of-opportunity")
Random.seed!(13)
tau = 150
res = runSimulationInfo(chart, um, thetaHat, D, m, IC=false, tau=tau, delta=0.65, maxrl=maxrl)
plt = plot_comparison(res, tau=tau)
safesave(folder * "/shaded-regions-cl", plt)

plt = plot_window_of_opportunity(res, tau=tau, jmpsize=0.1)
safesave(folder * "/thetahat", plt)

Random.seed!(123)
yinit = rand(D, m)
thetaHat = mean(yinit)
um = AdaptiveEstimator()
cfg = SimulationSettings(cfg, um = um)
L_gicp = runGICP(cfg, verbose=true)

chart = typeof(ch)(ch, L=L_gicp["res"][:L])                        

Random.seed!(13)
tau = 150
res = runSimulationInfo(chart, um, thetaHat, D, m, IC=false, tau=tau, delta=0.5, maxrl=maxrl)
plt = plot_comparison(res, tau=tau)
safesave(folder * "/shaded-regions-cl", plt)

plt = plot_window_of_opportunity(res, tau=tau, jmpsize=0.1)
safesave(folder * "/thetahat", plt)
