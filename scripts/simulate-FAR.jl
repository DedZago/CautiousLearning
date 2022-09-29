using Distributed
addprocs(10)
@everywhere using DrWatson
@everywhere @quickactivate "CautiousLearning"

@everywhere using DataFrames, SharedArrays
@everywhere using Plots, StatsPlots, GLM, LaTeXStrings, RCall, Random
@everywhere include(srcdir("cfg.jl"))

@everywhere function simulate_FAR(ch, θ, m; Nrep = 10000, delta=0.0, tau=0, exact = false)
    # Average update delay when process is IC
    IC = delta == 0.0 ? true : false

    delta_sh = standardize_shift(delta, Poisson(θ))
    y0 = gen_ic_oc_data(m, Poisson(θ), IC = IC, delta=delta, tau=tau)
    Ntrue = Nrep
    # FAR = get_far(ARL0)
    rej = 0
    trial = 0
    j = 0
    if exact
        thetaHat = θ
    else
        thetaHat = mean(y0)
    end
    while j < Nrep
        # Observe data
        y = rand(Poisson(θ + delta_sh))

        # ----- Bootstrap control limits
        updated_chart = update_series(ch, chart_statistic(y, thetaHat))
        # @inbounds Rvec[i] = R

        # Check for alarm
        if check_OC(updated_chart)
            rej += 1
            j += 1
            trial = 0
        else
            j += 1
            trial = 0
        end
    end
    return rej / Ntrue
end

config = Dict(
    :t => 1,
    :θ => 4.0,
    :ms => 50,
    :tau => [collect(0:10:100)],
    :delta => [0.5],
    :L => 0.5,
    :chart => ["EWMA"],
    :l => [0.2]
)

config_list = dict_list(config)

using GLM, LaTeXStrings, RCall, Random, KernelDensity
using Plots.PlotMeasures
include(scriptsdir("plot.jl"))

Random.seed!(2022-04-20)
out = DataFrame(a=Float64[], m=Int64[], bias=Float64[], v=Float64[], lambda=Float64[])

for i in 1:length(config_list)
    pars = config_list[i]
    file = plotsdir("FAR", savename(pars, "df.jld2"))
    if !isfile(file)
        @unpack t,θ,ms,L,chart,l,tau,delta = pars
        ch = eval(Symbol(chart))(l = l, L=L)
        println("Calculating true alpha...") 
        alpha_true = simulate_FAR(ch, θ, 10, Nrep = Int(1e08), exact=true) 

        println("Done.")


        println("Simulating chart alphas...")
        len_tau = length(tau)
        # ToDo: from here 
        nsim = 5000
        nrep = 40000
        alpha_shared = SharedArray{Float64}(len_tau, nsim)
        @sync @distributed for j in 1:len_tau
            println("j = " * string(j) * "/" * string(len_tau))
            alpha_sim = [simulate_FAR(ch, θ, ms + tau[j], Nrep = nrep, tau=ms, exact=false, delta=delta) for _ in 1:nsim]
            alpha_shared[j, :] = alpha_sim
        end
        
        alpha = convert(Array, alpha_shared)
        println("Done.")

        df = DataFrame(hcat(tau, alpha), :auto)
        sv = tostringdict(Dict(pairs(eachcol(df))))
        save(file, sv)
    else
        println("File already exists.")
        tmp = load(file)
        df = DataFrame(tmp)
    end
        # cur_colors = theme_palette(:auto)
        cur_colors = get_color_palette(:auto, 17)[1]
        pl = plot(xlab=L"\pi", ylab=L"n_{\textrm{oc}}", camera = (25, 20), right_margin=6mm, bottom_margin=4mm, dpi = 300)
        for j in 1:nrow(df)
            dens = kde(Vector(df[j, 2:end]))
            plot!(pl, dens.x,  fill(df[j, 1], length(dens.x)), dens.density, label="", colour=cur_colors)   
        end
        plot!(pl, fill(alpha_true, nrow(df)), df[:, 1], fill(0.0, nrow(df)), line = (:black, :dot, 2.0), label=L"\alpha_1")
        pl
        safesave(plotsdir("FAR", savename(pars,"power-convergence.png")), pl)
        
        pl = plot(df[:, 1], [mean(row[2:end]) for row in eachrow(df)], xlab=L"n_{\textrm{oc}}", ylab=L"\mathbb{E}[\pi]", dpi = 300, label="")
        hline!(pl, [alpha_true], line = (:black, :dot, 2.0), label=L"\alpha_1")
        safesave(plotsdir("FAR", savename(pars,"expected-power-convergence.png")), pl)
end