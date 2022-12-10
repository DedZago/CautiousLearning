println("Activating environment...")
using DrWatson
@quickactivate "CautiousLearning"

include(scriptsdir("simulate", "load_cfg.jl"))
verbose = false
println("Starting simulations...")
for cfg in config
    folder = "theta="*string(params(cfg.D)[1])*"_"*string(cfg.ch)
    println("Starting ", cfg.um)
    start_index = 1
    if add_simulations
        # If add_simulations, add them to already-present ones
        current_folder = datadir("sims", folder, string(cfg.um), "arls")
        start_index = length(readdir(current_folder)) + 1
    end
    @showprogress for i in start_index:(start_index + nsim_each)
        @show i
        sleep(0.5)
        cfg_thread = SimulationSettings(cfg, simulation = i + (current_thread-1) * nsim_each)
        loadname = datadir("sims", folder, string(cfg_thread.um), "limits", savename(cfg_thread, "jld2"))
        if isfile(loadname)
            GICPoutput = load(loadname)
            svname = datadir("sims", folder, string(cfg_thread.um), "arls", savename(cfg_thread, "jld2"))
            if isfile(svname)
                continue
            else
                out = runConditionalSimulations(GICPoutput, verbose=false)
                sv = @strdict config = cfg_thread out
                safesave(svname, sv)
                GC.gc()
            end#if
                
        else
            @warn "File " * loadname * " does not exist. Skipping."
            continue
        end#if
    end
end
