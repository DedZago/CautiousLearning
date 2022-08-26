println("Activating environment...")
using DrWatson
@quickactivate "CautiousLearning"

include(scriptsdir("simulate", "load_cfg.jl"))
verbose = false
println("Starting simulations...")
for cfg in config
    folder = "theta="*string(params(cfg.D)[1])*"_"*string(cfg.ch)
    println("Starting ", cfg.um)
    @showprogress for i in 1:nsim_each
        sleep(0.5)
        cfg_thread = SimulationSettings(cfg, simulation = i + (current_thread-1) * nsim_each)
        svname = datadir("sims", folder, string(cfg_thread.um), "limits", savename(cfg_thread, "jld2"))
        if isfile(svname)
            continue
        else
            sv = runGICP(cfg_thread, verbose=false)
            safesave(svname, sv)
            GC.gc()
        end
    end
end
