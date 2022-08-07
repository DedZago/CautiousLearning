
function chart_statistic(x, thetaHat)
    return (x - thetaHat)/sqrt(thetaHat)
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

    return @NamedTuple{t_alarm::Int64, chart_values::Vector{Float64}, limit_alarm::Float64, limit_cautious::Float64, parameter_updates::Vector{Float64}, di::Vector{Int64}}((t_alarm = i, chart_values = valueVec[1:i], limit_alarm = get_limits(ch), limit_cautious = get_warning_limit(um), parameter_updates = thetaHatCautVec[1:i], di = diVec[1:i]))

end

function runExperiment(ch, um, D, m; IC=true, tau=1, delta=0.0, maxrl=1e04)
    yinit = rand(D, m)
    thetaHat = mean(yinit)

    return runSimulation(ch, um, thetaHat, D, m, IC=IC, tau=tau, delta=delta, maxrl=maxrl)
end


#! TODO: adapt to work with control charts
"""
	stochasticApproximationPCC(hstart::Real, rlsim::Function, Arl0::Real; kw...)

Computes the control limit for a control chart such that it satisfies E[RL] = Arl0.
"""
#? Set hstart to CH
#? 
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
    
    if verbose println("Done.") end

    return (h=hm, iter=i)
end