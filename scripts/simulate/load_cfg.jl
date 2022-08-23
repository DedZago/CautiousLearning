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