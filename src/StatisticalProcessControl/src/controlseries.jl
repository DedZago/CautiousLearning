using StatisticalProcessControl
using Parameters

function check_IC(ch::CC) where CC <: DoubleSidedUnivariateSeries
    return -ch.L <= ch.C <= ch.L
end

function check_IC(ch::CC) where CC <: OneSidedUnivariateSeries
    if ch.upw
        return ch.C <= ch.L
    else
        return ch.C >= ch.L
    end
end

function check_OC(ch::CC) where CC <: UnivariateSeries
    return !check_IC(ch)
end

@with_kw struct EWMA <: DoubleSidedUnivariateSeries
    l::Float64 = 0.1
    L::Float64 = 0.0
    C::Float64 = 0.0
end

get_params(ch::EWMA) = @NamedTuple{l::Float64, L::Float64}((ch.l, ch.L))
get_value(ch::EWMA) = ch.C

function update_value(ch::EWMA, x)
    return (1.0 - ch.l) * get_value(ch) + ch.l * x
end

update_series(ch::EWMA, x) = EWMA(ch; C = update_value(ch, x))

@with_kw struct signedEWMA <: OneSidedUnivariateSeries
    l::Float64 = 0.1
    upw::Bool  = true
    L::Float64 = 0.0
    C::Float64 = 0.0
end

get_params(ch::signedEWMA) = @NamedTuple{l::Float64, upw::Bool, L::Float64}((ch.l, ch.upw, ch.L))
get_value(ch::signedEWMA) = ch.C

function update_value(ch::signedEWMA, x)
    if ch.upw
        return max(0.0, (1.0 - ch.l) * get_value(ch) + ch.l * x)
    else
        return min(0.0, (1.0 - ch.l) * get_value(ch) + ch.l * x)
    end
end

update_series(ch::signedEWMA, x) = signedEWMA(ch; C = update_value(ch, x))

##########################################################
#                       AEWMA                            #
##########################################################

# Default use Huber loss function
function huber(e, l, k)
    if e < -k
        return e + (1.0 - l)*k
    elseif e > k
        return e - (1.0 - l)*k
    else
        return l*e
    end
end

@with_kw struct AEWMA <: DoubleSidedUnivariateSeries
    l::Float64 = 0.1
    k::Float64 = 3.0
    L::Float64 = 0.0
    C::Float64 = 0.0
end

get_params(ch::AEWMA) = @NamedTuple{l::Float64, k::Float64, L::Float64}((ch.l, ch.k, ch.L))

function update_value(ch::AEWMA, x)
    e = x - get_value(ch)
    return get_value(ch) + huber(e, ch.l, ch.k)
end

update_series(ch::AEWMA, x) = AEWMA(ch; C = update_value(ch, x))

@with_kw struct signedAEWMA <: OneSidedUnivariateSeries
    l::Float64 = 0.1
    k::Float64 = 3.0
    upw::Bool = true
    L::Float64 = 0.0
    C::Float64 = 0.0
end

get_params(ch::signedAEWMA) = @NamedTuple{l::Float64, k::Float64, upw::Bool, L::Float64}((ch.l, ch.k, ch.upw, ch.L))

function update_value(ch::signedAEWMA, x)
    e = x - get_value(ch)
    if ch.upw
        return max(0.0, get_value(ch) + huber(e, ch.l, ch.k))
    else
        return min(0.0, get_value(ch) + huber(e, ch.l, ch.k))
    end
end

update_series(ch::signedAEWMA, x) = signedAEWMA(ch; C = update_value(ch, x))

#######################################################################
#                           Application                               #
#######################################################################


"""
Apply a control series to a vector x of observations

Parameters
----------
ch : A Series object
x : A vector of values

Returns
-------
"""

function apply_series(ch, x::Vector)
    z = similar(x)
    
    ch = update_series(ch, x[1])
    @inbounds z[1] = get_value(ch)
    for i in 2:length(z)
        ch = update_series(ch, x[i])
        @inbounds z[i] = get_value(ch)
    end
    return z
end

# mutable struct OnlineProcess{Z, C <: AbstractSeries}
#     z::Z
#     ch::C
# end

# OnlineProcess(ch::AbstractSeries) = OnlineProcess(init(ch), ch)

# function update!(OP::OnlineProcess, x)
#     OP.z = update_series(OP.ch, OP.z, x)
# end

# function online_monitor(ch::AbstractSeries, x::AbstractVector)
#     OP = OnlineProcess(ch)
#     for i in eachindex(x)
#         update!(OP, x[i])
#     end
#     return OP
# end

# --- Bootstrap limits
#abstract type Limits end
#abstract type BootstrapLimits <: Limits end

#struct UnivariateBootstrapLimits{C} <: BootstrapLimits where OP <: OnlineProcess
#    Rstar::Vector{OP}
#    ARL0::Integer
#    upper::Float64
#    lower::Float64
#end

## Constructor with series + generator
#function UnivariateBootstrapLimits(ch::C, ARL0::Integer) where C <: AbstractSeries
#    UnivariateBootstrapLimits([OnlineProcess(ch) for _ in 1:n], ARL0, )
#end

#compute_far(ARL0) = 1.0 - 1.0/ARL0

#function update_limits(LI::L, x, gen::Function) where L <: BootstrapLimits
#    # lower = get_lower(LI)
#    upper = get_upper(LI)
#    Rstar = get_seriess(LI)
#    n = length(Rstar)
#    if n == 0
#        throw(DomainError(n, "Bootstrap replicates are empty"))
#    end

#    ARL0 = get_ARL0(Rstar[1])
#    # Sample values that are below previous control limit
#    Rstar = sample(Rstar[get_values(LI) .<= upper], n, replace=true)
#    # Generate data from IC model
#    xstar = gen(n)

#    # Update bootstrap seriess
#    for i in 1:n
#        @inbounds Rstar[i] = update_series(Rstar[i], xstar[i])
#    end

#    vals = [get_value(ch) for ch in Rstar]
    
#    # Get new limits
#    upper = quantile(vals, 1.0 - (1.0 / ARL0))
#    return UnivariateBootstrapLimits(Rstar, upper)
#end

#function alarm(ch::C, lim::L) where {C <: AbstractSeries, L <: Limits}
#    if get_value(ch) < get_upper(lim)
#        return false
#    else
#        return true
#    end
#end

##----- Example usage
#ARL0 = 200
#n = 25*ARL0
#sim(n) = randn(n)
#ch = AEWMA(0.1, 3.0)
#chs = [AEWMA(0.1, 3.0) for _ in 1:n]
#lim = UnivariateBootstrapLimits(ch, n)

#x = 0.1
#@code_warntype update_limits(lim, x, sim)