function chart_statistic(x, thetaHat)
    return (x - thetaHat)/sqrt(thetaHat)
end

function runSimulation(ch, um, thetaHat, D, m; IC=true, tau=1, delta=0.0, maxrl=Int(1e04), seed = 123)
    # Random.seed!(seed)
    maxrl_i = Int(round(maxrl))
    thetaHatVec = zeros(maxrl_i)
    thetaHatVec[1] = thetaHat
    di = 1
    diVec = Array{Int}(undef, maxrl_i)
    diVec[1] = di
    thetaHatCaut = thetaHat
    thetaHatCautVec = zeros(maxrl_i)
    thetaHatCautVec[1] = thetaHatCaut

    valueVec = zeros(maxrl_i)
    valueVec[1] = get_value(ch)
    i = 1
    while i < maxrl_i
        thetaHatCaut = thetaHatVec[i - di + 1]
        y = gen_data_seq(D, i, IC=IC, tau=tau, delta=delta)
        ch = update_series(ch, chart_statistic(y, thetaHatCaut))
        thetaHat = update_parameter(thetaHat, y, i + m)
        i += 1
        valueVec[i] = get_value(ch)
        thetaHatVec[i] = thetaHat
        if check_update(ch, um)
            di = 1
        else
            di += 1
        end
        diVec[i] = di
        thetaHatCautVec[i] = thetaHatCaut
        if check_OC(ch)
            break
        end

    end

    # return @NamedTuple{t_alarm::Int64, chart_values::Vector{Float64}, limit_alarm::Float64, limit_cautious::Float64, parameter_updates::Vector{Float64}, di::Vector{Int64}}((t_alarm = i, chart_values = valueVec[1:i], limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um), parameter_updates = thetaHatCautVec[1:i], di = diVec[1:i]))
    return @NamedTuple{t_alarm::Int64}((t_alarm = i, ))

end


"""
	saControlLimits(ch, um, rlsim::Function, Arl0, thetaHat, Dist, m; kw...)

Computes the control limit for a control chart `ch` with update mechanism `um` such that it satisfies E[RL] = Arl0.

Setting `adjusted = false` ensures that E[RL] >= Arl0.
"""
function saControlLimits(ch, um, rlsim::Function, Arl0, thetaHat, Dist, m; Nfixed=500, Afixed=0.1, Amin=0.1, Amax=100, delta=0.1, q=0.55, gamma=0.02, Nmin=1000, z = 3.0, Cmrl=10.0, maxiter = 4e05, verbose=true, eps = 1e-04, adjusted=false, seed = Int(rand(1:1e06)))
    v = (z/gamma)^2
    sm = 0.0
    sp = 0.0
    u = 0.0

    h = get_limits(ch)
    rm = rp = hm = score = s2 = D = 0.0

    # Skip adaptation
    D = Amin
    i = 0

    if adjusted
        # Correction so that the resulting Arl0 estimate overshoots the true Arl0 with high probability
        # instead of being centered around the true Arl0 (Equation 12 of Capizzi & Masarotto (2016)).
        Arl0 = Arl0 / (1.0 - gamma)
    end
    
    if verbose println("Running SA ...") end

    while i < maxiter
        if verbose && i % floor(maxiter / 20) == 0
            println("i: ", i, "/", Integer(maxiter))
        end
        i += 1
        ch_temp = typeof(ch)(ch, L = h)                # Copy chart ch with hstart as control limit
        rl = rlsim(ch_temp, um, thetaHat, Dist, m, IC=true, maxrl = Int(round(Cmrl * Arl0 * sqrt(i + Nfixed))), seed=seed)[:t_alarm]
        score = (rl - Arl0)/Arl0

        h = max(eps, h - D * score / (i^q))
        hm = hm + (h - hm) / float(i)
        s2 = s2 + (score * score - s2) / float(i)

        if (i > Nmin) && (i > v * s2)
            if verbose println("Convergence!") end
            break
        end
    end
    
    if verbose println("saCL done.") end

    return (h=hm, iter=i)
end

function GICP(ch, um, yinit, thetaHat, runSimulation, m, Arl0; beta = 0.2, gamma = 0.03, verbose=true, seed = Int(rand(1:1e06)))
    #! Does not use seed because of conflicts
    #! Cannot optimize on the same value of seed, need to find a better mechanism
    CI = get_confint(yinit, conf = 1.0 - beta)

    if verbose println("Calculating GICP for lower extremum...") end
    Llow = saControlLimits(ch, um, runSimulation, Arl0, thetaHat, Poisson(CI[1]), m,
                            verbose=false, Amin=0.1, maxiter=1e05, gamma=gamma, q=0.55, adjusted=true)

    if verbose println("Calculating GICP for upper extremum...") end
    Lup = saControlLimits(ch, um, runSimulation, Arl0, thetaHat, Poisson(CI[2]), m,
                            verbose=false, Amin=0.1, maxiter=1e05, gamma=gamma, q=0.55, adjusted=true)
    if verbose println("GICP done.") end

    L = max(Llow[:h], Lup[:h])
    return L
end

function adjust_chart_gicp(ch::C, um, yinit, thetaHat, runSimulation, m, Arl0; beta = 0.2, gamma=0.03, verbose=true, seed = Int(rand(1:1e06))) where C <: UnivariateSeries
    L = GICP(ch, um, yinit, thetaHat, runSimulation, m, Arl0, beta = beta, gamma=gamma, verbose = verbose, seed=seed)
    return typeof(ch)(ch, L = L)
end

function extract_statistics(conditionalRun)
    return [cr[:t_alarm] for cr in conditionalRun]
end

function runConditionalSimulations(yinit, ncond, ch, um, D, m, Arl0; beta::Union{Bool, Float64} = 0.2, IC=true, tau=1, delta=0.0, maxrl=1e04, verbose=true, seed=Int(rand(1:1e06)))
    thetaHat = mean(yinit)
    rlIC = zeros(ncond)
    rlOC = [[zeros(ncond) for t in tau] for d in delta if false in IC]
    if isa(um, CautiousLearning)
        # Estimate cautious learning limit
        if verbose println("Calculating limits for target ATS...") end
        Ats0 = get_ATS(um)
        sa_ats = saControlLimits(ch, SelfStarting(), runSimulation, Ats0, thetaHat, Poisson(thetaHat),
                                m, verbose=false, Amin=0.1, maxiter=1e05,
                                gamma=0.015, adjusted=true, seed=seed)
        um = CautiousLearning(L = sa_ats[:h], ATS = Ats0)
    end

    if beta == false
        chart = deepcopy(ch)
    else
        chart = adjust_chart_gicp(ch, um, yinit, thetaHat, runSimulation, m, Arl0, beta=beta, verbose=verbose, seed=seed + 1)
    end
    
    if verbose println("Simulating run lengths...") end
    if true in IC
        seed_offset = 2
        icRun = [runSimulation(chart, um, thetaHat, D, m, IC=true, tau=1, delta=0.0, maxrl=maxrl, seed=seed+seed_offset+s) for s in 1:ncond]
        output = extract_statistics(icRun)
        rlIC = output
    end

    if false in IC
        for t in eachindex(tau)
            for d in eachindex(delta)
                seed_offset = 2 + ncond + t*length(delta) + d
                ocRun = [runSimulation(chart, um, thetaHat, D, m, IC=false, tau=tau[t], delta=delta[d], maxrl=maxrl, seed=seed+seed_offset+s) for s in 1:ncond]
                output = extract_statistics(ocRun)
                rlOC[d][t] = output
            end#for
        end#for
    end#if
    return (IC = rlIC, OC = rlOC)
end

function condSimToDf(out, tau, delta)
    # Save everything as a dataframe
    colnames = ["tau", "delta", "rl"]
    dfOutput = DataFrame([name => [] for name in colnames])
    rlIC = out[:IC]
    push!(dfOutput, [0.0, 0.0, rlIC])
    for t in eachindex(tau)
        for d in eachindex(delta)
            rlOC = out[:OC][d][t]
            push!(dfOutput, [tau[t], delta[d], rlOC])
        end#for
    end#for

    return dfOutput
end

# function runNestedSimulations(ncond, yinit, ch, um, D, m, Arl0; beta::Union{Bool, Float64} = 0.2, IC=true, tau=1, delta=0.0, maxrl=1e04, verbose=true, seed=Int(rand(1:1e06)))
#     nsim = length(yinit)

#     out = pmap((i) -> runConditionalSimulation(yinit[i], ncond, ch, um ,D, m, Arl0, beta=beta, IC=IC, tau=tau, delta=delta, maxrl = maxrl, verbose=verbose, seed=seed+i), 1:length(yinit))

#     colnames = ["tau", "delta", "rl"]
#     dfOutput = DataFrame([name => [] for name in colnames])
#     rlIC = permutedims(hcat([v[:IC] for v in out]...))
#     push!(dfOutput, [0.0, 0.0, rlIC])
#     for t in eachindex(tau)
#         for d in eachindex(delta)
#             rlOC = permutedims(hcat([v[:OC][d][t] for v in out]...))
#             push!(dfOutput, [tau[t], delta[d], rlOC])
#         end#for
#     end#for

#     return dfOutput
# end

function runExperiment(config::SimulationSettings)
    ncond   = config.ncond
    ch      = config.ch
    um      = config.um
    D       = config.D
    m       = config.m
    Arl0    = config.Arl0
    beta    = config.beta
    IC      = config.IC
    tau     = config.tau
    delta   = config.delta
    maxrl   = config.maxrl
    verbose = config.verbose
    simulation = config.simulation
    
    seed    = config.seed + simulation
    Random.seed!(seed)
    yinit = rand(D, m)
    out = runConditionalSimulations(yinit, ncond, ch, um, D, m, Arl0, beta=beta, IC=IC, delta=delta, tau=tau, verbose=verbose, maxrl=maxrl, seed=seed)
    if verbose println("Done.") end
    df = condSimToDf(out, tau, delta)
    sim = [config.simulation for _ in 1:nrow(df)]
    df = hcat(sim, df)
    rename!(df, :x1 => :sim)
    return df
end
