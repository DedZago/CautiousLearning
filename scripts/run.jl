println("Activating environment...")
using DrWatson
@quickactivate "CautiousLearning"

println("Loading args...")
using ArgParse, ProgressMeter
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


println("Sourcing configuration...")
include(srcdir("cfg.jl"))
verbose = false

println("Starting simulations...")
for cfg in config
    folder = "theta="*string(theta)*"_"*string(ch)
    println("Starting ", cfg.um)
    @showprogress for i in 1:nsim_each
        sleep(0.5)
        cfg_thread = SimulationSettings(cfg, simulation = i + (current_thread-1) * nsim_each)
        svname = datadir("sims", folder, string(cfg.um), savename(cfg_thread, "jld2"))
        if isfile(svname)
            continue
        else
            out = runExperiment(cfg_thread, verbose=verbose)

            sv = @strdict config = cfg_thread out
            safesave(svname, sv)
            GC.gc()
        end
    end
end
