using StatisticalProcessControl
using Parameters

struct AdaptiveEstimator <: AbstractUpdate end
function check_update(::C, ::AdaptiveEstimator) where C <: UnivariateSeries
    return true
end

get_warning_limit(um::AdaptiveEstimator) = Inf
get_ATS(um::AdaptiveEstimator) = Inf

struct FixedParameter <: AbstractUpdate end
function check_update(::C, ::FixedParameter) where C <: UnivariateSeries
    return false
end

get_warning_limit(um::FixedParameter) = 0.0
get_ATS(um::FixedParameter) = 1


@with_kw struct CautiousLearning <: AbstractUpdate
    L = 0.5
    ATS::Int = 3.0
end

get_warning_limit(um::CautiousLearning) = um.L
get_ATS(um::CautiousLearning) = um.ATS

"""
	check_update(ch, um)

Returns true if the chart is inside the update region, false if the chart is inside the warning region.
"""
function check_update(ch::C, um::CautiousLearning) where C <: DoubleSidedUnivariateSeries
    new_ch = typeof(ch)(ch, L = get_warning_limit(um))
    return check_IC(new_ch)
end

function check_update(ch::C, um::CautiousLearning) where C <: OneSidedUnivariateSeries
    ATS = get_ATS(um)
    if ATS > 0
        new_ch = typeof(ch)(ch, L = get_warning_limit(um))
        return check_IC(new_ch)
    else
        return get_value(ch) == 0.0
    end
end