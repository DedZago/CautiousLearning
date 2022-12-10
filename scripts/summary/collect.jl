using DrWatson
@quickactivate "CautiousLearning"
include(srcdir("cfg.jl"))

function collect_limits(folder)
    df = DataFrame()
    for folder_sims in filter(isdir, readdir(folder, join=true))
        if basename(folder_sims) == "output"
            continue
        end
        sims = collect_results(folder_sims * "/limits")
        ums = [string(sims.cfg[i].um) for i in 1:nrow(sims)]
        simulations = [sims.cfg[i].simulation for i in 1:nrow(sims)]
        thetas = [sims.res[i].thetaHat for i in 1:nrow(sims)]
        lims = [sims.res[i].L for i in 1:nrow(sims)]
        append!(df, DataFrame(hcat(ums, simulations, thetas, lims), [:um, :sim,:thetaHat, :L]))
    end
    sort!(df, [:sim])
    return df
end

function sims_to_dataframe(folder)
    df = DataFrame()
    for folder_sims in filter(isdir, readdir(folder, join=true))
        if basename(folder_sims) == "output" || !isdir(folder_sims * "/arls")
            continue
        end
        sims = collect_results(folder_sims * "/arls")
        if nrow(sims) == 0
            continue
        end
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

overwrite = true
using CSV

for folder in filter(isdir, readdir(datadir("sims"), join=true))
    if isdir(datadir("sims", folder, "output")) && !overwrite
        continue
    else
        output = collect_limits(folder)
        safesave(datadir("sims", folder, "output", "limits.jld2"), @strdict output)
        CSV.write(datadir("sims", folder, "output", "limits_R.csv"), output)
    end
end

for folder in filter(isdir, readdir(datadir("sims"), join=true))
    output = sims_to_dataframe(folder)
    # output = vcat([sims_to_dataframe(folder) for cfg in config[1:3]]...)
    output.ARL = [mean(output.rl[i]) for i in 1:nrow(output)]
    output.SDRL = [std(output.rl[i]) for i in 1:nrow(output)]
    safesave(datadir("sims", folder, "output", "output.jld2"), @strdict output)

    output_R = output[:, Not(:rl)]
    CSV.write(datadir("sims", folder, "output", "output_R.csv"), output_R)
end
