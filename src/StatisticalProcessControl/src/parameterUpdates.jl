using DrWatson
@quickactivate "cautiousBootstrap"

using Random, Distributions

include(projectdir("src", "simulate_data.jl"))      # generate Poisson data
include(projectdir("src", "spc-functions.jl"))      # control charts
include(projectdir("src", "cautious.jl"))           # cautious learning functions

function parameterUpdates(θhat, m, dist, UPR; IC = true, τ = 0.0, δ = 0.0, maxiter = 5000)
    θhatVec = zeros(maxiter)
    θhatCautVec = zeros(maxiter)
    dtVec = zeros(Int64, maxiter)
    dt = 1
    i = 1
    θhatVec[1] = θhat
    θhatCautVec[1] = θhat
    θhatCaut = θhat
    dtVec[1] = dt
    yVec = gen_ic_oc_data(maxiter, dist, IC=IC, τ=τ, δ=δ)

    isBoot = isa(UPR, TwoSidedCautiousBootstrap)
    if isBoot
        Nboot = 2000
        Rstar = fill(eps(), Nboot)                      # bootstrap control chart values
        # ch = AEWMA(λ=0.15, k=3.0)
        ch = AEWMA(λ=0.15, k = 3.0)
        R = eps()               # initial chart value
    end
    distType = (typeof(dist).name).wrapper          # process distribution name
    auxpars = []
    if distType == Binomial
        push!(auxpars, params(dist)[1])
    end

    UPR_cp = deepcopy(UPR)
    predDist = distType([auxpars; [θhatCaut]]...)

    while i < maxiter
        θhatCaut = θhatVec[i-dt+1]

        predDist = distType([auxpars; [θhatCaut]]...)
        # θhat_t = θhat(y1, ..., y_t-1)
        # Observe data
        y = yVec[i]

        # Calculate likelihood on y_{i-di+1}:y_i to get index
        t = get_t(UPR)
        # Use dynamic cautious learning

        # Update θ_i+1 using y_1, ..., y_i
        θhat = update_parameter(θhat, y, t, predDist)

        if isBoot 
            # ----- Bootstrap control limits
            ystar = Float64.(rand(predDist, Nboot))
            # @show ystar
            for b in 1:Nboot
                @inbounds Rstar[b] = update_chart(ch, Rstar[b], chart_statistic(ystar[b], predDist))
            end
            R = update_chart(ch, R, chart_statistic(y, predDist))
            # @show Rstar
            update_rule!(UPR_cp, y, Rstar, R)
        else
            update_rule!(UPR_cp, y)
        end
        dt = get_dt(UPR_cp)
        i += 1
        @inbounds dtVec[i] = dt
        # True sequence of θhat
        @inbounds θhatVec[i] = θhat
        @inbounds θhatCautVec[i] = get_cautious_μ(UPR_cp)  
    end
    return dtVec, θhatCautVec, Rstar
end


# Random.seed!(2)
# θ = 5.0
# # nb = 10
# # dist = Bernoulli(θ)
# D = Poisson
# dist = D(θ)
# # dist = Binomial(10, 0.5)
# m = 50
# # θhat = mean(rand(dist, m))/nb
# θhat = estimate_parameter(rand(dist, m), dist) 
# IC = false
# τ = 100
# δ = standardize_shift(0.5, dist)
# dyn = false
# maxiter = 200
# ATH = 3.0
# UPR = TwoSidedCautiousBootstrap(D(θhat), ATH, θhat, 1, m)
# res = parameterUpdates(θhat, m, dist, UPR, IC = IC, τ = τ, δ = δ, maxiter = maxiter)

# using StatsPlots, LaTeXStrings
# include(scriptsdir("plot.jl"))      # generate Poisson data
# t = 1:length(res[1])
# p = plot(t, t .- res[1], dpi=600, xlab=L"t", ylab=L"t - d_t", label="",legend=:outerright)
# if δ > 0.0
#     vline!(p, [τ], style=:dash, label=L"\tau")
# end
# safesave(plotsdir("CL-bootstrap-single", "t-minus-dt.png"), p)


# p = plot(res[2], xlab = L"t", ylab=L"\widehat{\vartheta}_t", label="", legend=:outerright, dpi=600)
# if δ > 0.0
#     vline!(p, [τ], style=:dash, label=L"\tau")
# end
# safesave(plotsdir("CL-bootstrap-single", "thetahat.png"), p)