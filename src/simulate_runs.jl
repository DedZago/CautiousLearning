function chart_statistic(x, thetaHat)
    return (x - thetaHat)/sqrt(thetaHat)
end

function chart_statistic(y, D::Distribution)
    (y - mean(D)) / std(D)
end

function runSimulation(ch, um, thetaHat, D, m; IC=true, tau=1, delta=0.0, maxrl=Int(1e04))
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

function runSimulation(ch, um::CautiousLearningCM, thetaHat, D, m; IC=true, tau=1, delta=0.0, maxrl=Int(1e04))
    A = um.A
    B = um.B
    qi = 0
    qtilde = 0

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

        # C&M (2020)
        qtilde = qi + chart_statistic(y, thetaHatCaut)^(2.0)
        if qtilde < A * di - B
            di = 1
            qi = 0
        else
            di += 1
            qi = qtilde
        end
        diVec[i] = di
        thetaHatCautVec[i] = thetaHatCaut
        if check_OC(ch)
            break
        end

    end

    return @NamedTuple{t_alarm::Int64, chart_values::Vector{Float64}, limit_alarm::Float64, parameter_updates::Vector{Float64}, di::Vector{Int64}}((t_alarm = i, chart_values = valueVec[1:i], limit_alarm = get_limits(ch), parameter_updates = thetaHatCautVec[1:i], di = diVec[1:i]))
    # return @NamedTuple{t_alarm::Int64}((t_alarm = i, ))

end


"""
	saControlLimits(ch, um, rlsim::Function, Arl0, thetaHat, Dist, m; kw...)

Computes the control limit for a control chart `ch` with update mechanism `um` such that it satisfies E[RL] = Arl0.

Setting `adjusted = false` ensures that E[RL] >= Arl0.
"""
function saControlLimits(ch, um, rlsim::Function, Arl0, thetaHat, Dist, m; Nfixed=500, Afixed=0.1, Amin=0.1, Amax=100, delta=0.1, q=0.55, gamma=0.02, Nmin=1000, z = 3.0, Cmrl=10.0, maxiter = 4e05, verbose=true, eps = 1e-04, adjusted=false)
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
        rl = rlsim(ch_temp, um, thetaHat, Dist, m, IC=true, maxrl = Int(round(Cmrl * Arl0 * sqrt(i + Nfixed))))[:t_alarm]
        score = (rl - Arl0)/Arl0

        h = max(eps, h - D * score / (i^q))
        hm = hm + (h - hm) / float(i)
        s2 = s2 + (score * score - s2) / float(i)

        if (i > Nmin) && (i > v * s2)
            if verbose println("Convergence!") end
            break
        end
    end
    
    if verbose println("SA done.") end

    return (h=hm, iter=i)
end

function GICP(ch, um, yinit, thetaHat, runSimulation, m, Arl0; beta = 0.2, gamma = 0.03, verbose=true)
    CI = get_confint(yinit, ch, conf = 1.0 - beta)

    
    Lvec = []
    CI_valid = [th for th in CI if th != 0.0 && th != Inf]
    for i in eachindex(CI_valid)
        if verbose println("Calculating GICP for extremum " * string(i) * "/" * string(length(CI_valid)) * "...") end
        L_th = saControlLimits(ch, um, runSimulation, Arl0, thetaHat, Poisson(CI_valid[i]), m,
                                verbose=false, Amin=0.1, maxiter=1e05, gamma=gamma, q=0.55, adjusted=true)
        push!(Lvec, L_th[:h])
    end
    if verbose println("GICP done.") end

    #! Use abs.(Lvec) to consider both positive and negative control limits
    L = maximum(abs.(Lvec))
    return L
end

function adjust_chart_gicp(ch::C, um, yinit, thetaHat, runSimulation, m, Arl0; beta = 0.2, gamma=0.03, verbose=true) where C <: UnivariateSeries
    L = GICP(ch, um, yinit, thetaHat, runSimulation, m, Arl0, beta = beta, gamma=gamma, verbose = verbose)
    return typeof(ch)(ch, L = L)
end

function extract_statistics(conditionalRun)
    return [cr[:t_alarm] for cr in conditionalRun]
end

function computeGICP(yinit, ch, um, D, m, Arl0; beta::Union{Bool, Float64} = 0.2, verbose=true)
    thetaHat = mean(yinit)
    if isa(um, CautiousLearning)
        Ats0 = get_ATS(um)
        if Ats0 != 0
            # Calculate limit if Ats0 != 0, otherwise use zero-restarting chart
            if verbose println("Calculating limits for target ATS...") end
            sa_ats = saControlLimits(ch, AdaptiveEstimator(), runSimulation, Ats0, thetaHat, Poisson(thetaHat),
                                    m, verbose=false, Amin=0.1, maxiter=1e05,
                                    gamma=0.015, adjusted=true)
            um = CautiousLearning(L = sa_ats[:h], ATS = Ats0)
        else
            if verbose println("ATS = 0, skipping limit calculation.") end
        end
        # Estimate cautious learning limit
    end

    if beta == false
        chart = deepcopy(ch)
    else
        chart = adjust_chart_gicp(ch, um, yinit, thetaHat, runSimulation, m, Arl0, beta=beta, verbose=verbose)
    end
    
    return (thetaHat = thetaHat, m = m, L = chart.L, um = um)
end

function runGICP(cfg::SimulationSettings; verbose=true)
    ch      = cfg.ch
    um      = cfg.um
    D       = cfg.D
    m       = cfg.m
    Arl0    = cfg.Arl0
    beta    = cfg.beta
    simulation = cfg.simulation
        
    # Set simulation seed
    Random.seed!(cfg.seed + simulation)
    yinit = rand(D, m)
    res = computeGICP(yinit, ch, um, D, m, Arl0; beta=beta, verbose=verbose)
    if verbose println("Done.") end
    return @strdict cfg res yinit
end

function runConditionalSimulations(GICPoutput; verbose=true)
    cfg = deepcopy(GICPoutput["cfg"])
    res = deepcopy(GICPoutput["res"])
    ncond   = cfg.ncond
    ch      = cfg.ch
    D       = cfg.D
    m       = cfg.m
    IC      = cfg.IC
    tau     = cfg.tau
    delta   = cfg.delta
    maxrl   = cfg.maxrl
    simulation = cfg.simulation

    L::Float64 = res.L
    thetaHat::Float64 = res.thetaHat
    m::Int = res.m
    um = res.um
    chart = typeof(ch)(ch, L = L)

    rlIC = zeros(ncond)
    rlOC = [[zeros(ncond) for t in tau] for d in delta if false in IC]

    Random.seed!(cfg.seed + simulation + 1)
    
    if verbose println("Simulating run lengths...") end
    if true in IC
        icRun = [runSimulation(chart, um, thetaHat, D, m, IC=true, tau=1, delta=0.0, maxrl=maxrl) for _ in 1:ncond]
        output = extract_statistics(icRun)
        rlIC = output
    end

    if false in IC
        for t in eachindex(tau)
            for d in eachindex(delta)
                ocRun = [runSimulation(chart, um, thetaHat, D, m, IC=false, tau=tau[t], delta=delta[d], maxrl=maxrl) for _ in 1:ncond]
                output = extract_statistics(ocRun)
                rlOC[d][t] = output
            end#for
        end#for
    end#if
    out = (IC = rlIC, OC = rlOC, L = L, thetaHat = thetaHat, simulation = simulation, m=m)
    df = condSimToDf(out, tau, delta, simulation)
    return df
end

function condSimToDf(out, tau, delta, simulation)
    # Save everything as a dataframe
    colnames = ["sim", "thetaHat", "tau", "delta", "rl", "L", "m"]
    simulation = out[:simulation]
    thetaHat = out[:thetaHat]
    L = out[:L]
    m = out[:m]
    dfOutput = DataFrame([name => [] for name in colnames])
    rlIC = out[:IC]
    push!(dfOutput, [simulation, thetaHat, 0.0, 0.0, rlIC, L, m])
    for t in eachindex(tau)
        for d in eachindex(delta)
            rlOC = out[:OC][d][t]
            push!(dfOutput, [simulation, thetaHat, tau[t], delta[d], rlOC, L, m])
        end#for
    end#for

    return dfOutput
end
