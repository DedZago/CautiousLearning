println("Activating environment...")
using DrWatson
@quickactivate "CautiousLearning"

#? Add some chosen values of delta to already-simulated outputs

println("Loading args...")
using ArgParse, ProgressMeter
# Add command line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--current_thread", "-t"
            help = "Current index of the thread to launch."
            arg_type = Int
            default = 1
            required = true
        "--nsim_each", "-n"
            help = "Number of simulations for each thread"
            arg_type = Int
            default = 1
            required = true
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
current_thread = parsed_args["current_thread"]
nsim_each = parsed_args["nsim_each"]

println("Sourcing configuration...")
include(srcdir("cfg.jl"))
verbose = false

delta_add = [0.25]

println("Starting simulations...")
for dir in readdir((datadir("sims")))
    for simdir in readdir(datadir("sims", dir))
        if basename(simdir) == "output"
            continue
        else
            toDo = readdir(datadir("sims", dir, simdir))
            num_toDo = length(toDo)
            for i in 1:nsim_each
                idx = i + (current_thread-1) * nsim_each
                if idx > num_toDo
                    continue
                else
                    file = toDo[idx]
                    res = load(datadir("sims", dir, simdir, file))
                    res["out"] = runExperiment(res, delta_add; verbose=true)
                    safesave(datadir("sims_up", dir, simdir, file), res)
                end
            end
        end
    end
end
