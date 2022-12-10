using Distributions
using StatisticalProcessControl

# Update Poisson/Binomial mean for each data point
update_parameter(thetaHat, y, t) = thetaHat + (y - thetaHat) / t

function update_parameter(thetaHat, y, t, D::Binomial)
    pars = params(D)
    thetaHat + (y - pars[1]*thetaHat) / (pars[1] * t)
end

function update_parameter(thetaHat, y, t, ::Poisson)
    thetaHat + (y - thetaHat) / t
end

function update_parameter(thetaHat, y, t, ::Bernoulli)
    thetaHat + (y - thetaHat) / t
end
