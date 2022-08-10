using DrWatson
@quickactivate "CautiousLearning"

using ArgParse
using ProgressMeter
# Add command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--nsimulations_each", "-n"
            help = "Number of simulations for each process"
            arg_type = Int
            default = 1
            required = true
        "--current_thread", "-t"
            help = "Current index of the thread to launch."
            arg_type = Int
            default = 1
            required = true
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
nsim_each = parsed_args["nsimulations_each"]
current_thread = parsed_args["current_thread"]

include(srcdir("cfg.jl"))

umVec = [SelfStarting()]
# umVec = [FixedParameter(), SelfStarting(), CautiousLearning(ATS=3)]
config = [SimulationSettings(um = um, seed = 2022-08-10) for um in umVec]
for cfg in config
    @showprogress for i in 1:nsim_each
        cfg_thread = SimulationSettings(cfg, simulation = i + (current_thread-1) * nsim_each)
        svname = datadir("sims", "test", savename(cfg_thread, "jld2"))
        if isfile(svname)
            continue
        else
            out = runExperiment(cfg_thread)

            sv = @strdict config = cfg_thread out
            safesave(datadir("sims", "test", savename(cfg, "jld2")), sv)
            GC.gc()
        end
    end
end