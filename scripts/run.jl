using DrWatson
@quickactivate "CautiousLearning"


include(srcdir("cfg.jl"))

nsim_each = parse(Int, ARGS[1])
current_thread = parse(Int, ARGS[2])
umVec = [FixedParameter()]
# umVec = [FixedParameter(), SelfStarting(), CautiousLearning(ATS=3)]
config = [SimulationSettings(um = um, seed = 2022-08-09) for um in umVec]
for cfg in config
    for i in 1:nsim_each
        cfg_thread = SimulationSettings(cfg, simulation = i + (current_thread-1) * nsim_each)
        svname = datadir("sims", "test", savename(cfg_thread, "jld2"))
        if isfile(svname)
            continue
        else
            out = runExperiment(cfg_thread)
            sv = @strdict config=cfg_thread out
            safesave(datadir("sims", "test", savename(cfg, "jld2")), sv)
        end
            
    end
end