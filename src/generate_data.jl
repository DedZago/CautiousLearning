using Distributions

function standardize_shift(delta, D::Poisson)
    return std(D) * delta
end

function get_shifted_distribution(D::Poisson, delta)
    p = params(D)
    p_sh = (p[1] + standardize_shift(delta, D), )
    return Poisson(p_sh...)
end

# function get_shifted_distribution(D::Bernoulli, δ)
#     p = params(D)
#     p_sh = (p[1] + δ, )
#     return Bernoulli(p_sh...)
# end

# function get_shifted_distribution(D::Binomial, δ)
#     p = params(D)
#     p_sh = (p[1], p[2] + δ)
#     return Binomial(p_sh...)
# end

function gen_ic_oc_data(n, D::Distribution; IC = true, tau::Int = 1, delta = 0.0)
    out = Vector{Float64}(undef, n)

    if IC && delta != 0.0
        @warn "IC set to true and delta set to " * string(delta) * ", generating IC data."
    end

    if IC
        # Generate all IC
        out = Float64.(rand(D, n))
    else
        # First τ are IC
        nic = 1:(tau-1)
        # Last n - τ are OC
        noc = tau:n
        # Apply shift to the distribution
        dist = get_shifted_distribution(D, delta)
        if length(nic) == 0
            out = Float64.(rand(dist, n))
        else
            out[nic] = Float64.(rand(D, length(nic)))
            out[noc] = Float64.(rand(dist, length(noc)))
        end
    end
    return out
end

function gen_data_seq(D::Distribution, i; IC=true, tau = 1, delta = 0.0)
    if IC && delta != 0.0
        @warn "IC set to true and delta set to " * string(delta) * ", generating IC data."
    end

    if IC || i < tau
        dist = D
    elseif i >= tau
        dist = get_shifted_distribution(D, delta)
    end

    out = rand(dist)
    return out
end


function get_confint(x; D = Poisson, conf=0.95, eps = 1e-08)
    thetaHat = mean(x)
    n = length(x)
    alpha = 1.0 - conf
    return (lower = max(eps, thetaHat + quantile(Normal(0,1), alpha/2.0)*sqrt(thetaHat/n)),
            upper = thetaHat + quantile(Normal(0,1), 1.0 - alpha/2.0)*sqrt(thetaHat/n)
            )
end