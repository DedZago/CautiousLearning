using DrWatson
@quickactivate "CautiousLearning"
include(srcdir("cfg.jl"))


function sims_to_dataframe(folder)
    df = DataFrame()
    for folder_sims in filter(isdir, readdir(folder, join=true))
        sims = collect_results(folder_sims)
        vectorized_df = [deepcopy(sims.out[i]) for i in 1:nrow(sims)]
        for n in fieldnames(typeof(sims.config[1]))
            if !(string(n) in ["simulation", "D", "IC", "ncond", "verbose", "delta", "tau"])
                for i in eachindex(vectorized_df)
                    vectorized_df[i][!, n] = [string(getfield(sims.config[i], n)) for _ in 1:nrow(vectorized_df[i])]
                end
            end
        end
        merged_df = vcat(vectorized_df...)
        append!(df, merged_df)
    end
    sort!(df, [:sim])

    #TODO: Append the configuration relative to each single simulation
    return df
end

overwrite = false
for folder in filter(isdir, readdir(datadir("sims"), join=true))
    if isdir(datadir("sims", folder, "output")) && !overwrite
        continue
    else
        output = sims_to_dataframe(folder)
        # output = vcat([sims_to_dataframe(folder) for cfg in config[1:3]]...)
        output.ARL = [mean(output.rl[i]) for i in 1:nrow(output)]
        output.SDRL = [std(output.rl[i]) for i in 1:nrow(output)]
        safesave(datadir("sims", folder, "output", "output.jld2"), @strdict output)

        using CSV
        output_R = output[:, Not(:rl)]
        CSV.write(datadir("sims", folder, "output", "output_R.csv"), output_R)

        using StatsPlots
        # Safesave plots image
        DrWatson._wsave(s, fig::T) where T <: Plots.Plot{Plots.GRBackend} = savefig(fig, s)

        p1 = boxplot()
        hline!(p1, [parse(Int, output.Arl0[1])], style=:dash, label=false)
        for um in unique(output.um)
            df = output[output.um .== um .&& output.delta .== 0.0 .&& output.tau .== 0.0, :]
            # CARL = [mean(df.rl[i]) for i in 1:nrow(df) if df.delta[i] == 0.0 && df.tau[i] == 0.0]
            boxplot!(p1, df.ARL, label=um, outliers=false)
            println("mean(df.ARL .<= df.ARL0: ", mean(df.ARL .<= parse.(Int, df.Arl0)))
        end
        safesave(datadir("sims", folder, "output", "CARL0.png"), p1)

        for d in unique(output.delta)[2:end]
            for t in unique(output.tau)[2:end]
                p = boxplot()
                hline!([t], style=:dash, label=false)
                for um in unique(output.um)
                    df = output[output.um .== um .&& output.delta .== d .&& output.tau .== t, :]
                    boxplot!(p, df.ARL, label=um, outliers=false)
                end
                name = "CARL_tau=" * string(t) * "_delta=" * string(d) * ".png"
                safesave(datadir("sims", folder, "output", name), p)
            end
        end
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
