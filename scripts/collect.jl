using DrWatson
@quickactivate "CautiousLearning"
include(srcdir("cfg.jl"))

function sims_to_dataframe(cfg)
    sims = collect_results(datadir("sims", folder, string(cfg.um)))
    df = deepcopy(sims.out[1])
    for i in 2:nrow(sims)
        append!(df, deepcopy(sims.out[i]))
    end
    sort!(df, [:sim])

    #TODO: Append the configuration relative to each single simulation
    for n in fieldnames(typeof(cfg))
        if !(string(n) in ["simulation", "D", "IC", "ncond", "verbose", "delta", "tau"])
            df[!, n] = [string(getfield(cfg, n)) for _ in 1:nrow(df)]
        end
    end
    return df
end

output = vcat([sims_to_dataframe(cfg) for cfg in config]...)
safesave(datadir("sims", folder, "output", "output.jld2"), @strdict output)

using StatsPlots
# Safesave plots image
DrWatson._wsave(s, fig::T) where T <: Plots.Plot{Plots.GRBackend} = savefig(fig, s)

p1 = boxplot()
# p2 = boxplot()
hline!([parse(Int, output.Arl0[1])], style=:dash, label=false)
for um in unique(output.um)
    df = output[output.um .== um, :]
    CARL = [mean(df.rl[i]) for i in 1:nrow(df) if df.delta[i] == 0.0 && df.tau[i] == 0.0]
    boxplot!(p1, CARL, label=um, outliers=false)

    # CPerf = mean(df.rl[i] .>= df.Arl0[i] for i in 1:nrow(df) if df.delta[i] == 0.0 && df.tau[i] == 0.0)
    # boxplot!(p2, CPerf, label=um, outliers=false)
end
safesave(datadir("sims", folder, "output", "CARL0.png"), p1)

delta_performance = 500 * 0.02
for d in unique(output.delta)[2:end]
    for t in unique(output.tau)[2:end]
        p = boxplot()
        hline!([t], style=:dash, label=false)
        for um in unique(output.um)
            df = output[output.um .== um, :]
            CARL = [mean(df.rl[i]) for i in 1:nrow(df) if df.delta[i] == d && df.tau[i] == t]
            boxplot!(p, CARL, label=um, outliers=false)
        end
        name = "CARL_tau=" * string(t) * "_delta=" * string(d) * ".png"
        safesave(datadir("sims", folder, "output", name), p)
    end
end

# for cfg in config
#     CARL = mean(df.rl[i] for i in 1:nrow(df) if df.delta[i] == 0.0 && df.tau[i] == 0.0)
#     open(string(cfg.um),"a") do io
#         println(io, "d: 0.0")
#         println(io, "ave(CARL): ", mean(CARL))
#         println(io, "std(CARL): ", std(CARL))
#     end

#     for d in unique(df.delta)[2:end]
#         CARL = mean(df.rl[i] for i in 1:nrow(df) if df.delta[i] == d && df.tau[i] == 1)
#         open(string(cfg.um),"a") do io
#             println(io, "d: ", d)
#             println(io, "ave(CARL): ", mean(CARL))
#             println(io, "std(CARL): ", std(CARL))
#         end
#     end

    
#     safesave(datadir("sims", folder, "output", "output.jld2"), @strdict df)
# end

# CARL0 = mean(df.rl[i] for i in 1:nrow(df) if df.delta[i] == 0.0)
# Arl0 = 200
# mean(CARL0 .<= Arl0)

# CARL1 = mean(mean(df.rl[i] for i in 1:nrow(df) if df.delta[i] == 1.0 && df.tau[i] == 1))
