using StatisticalProcessControl
using Parameters

struct SelfStarting <: AbstractUpdate end
function check_update(::C, ::SelfStarting) where C <: UnivariateSeries
    return true
end

struct FixedParameter <: AbstractUpdate end
function check_update(::C, ::FixedParameter) where C <: UnivariateSeries
    return false
end


@with_kw struct CautiousLearning <: AbstractUpdate
    L = 0.0
end

get_warning_limit(um::AbstractUpdate) = um.L

"""
	check_update(ch, um)

Returns true if the chart is inside the update region, false if the chart is inside the warning region.
"""
function check_update(ch::C, um::CautiousLearning) where C <: UnivariateSeries
    new_ch = typeof(ch)(ch, L = get_warning_limit(um))
    return check_IC(new_ch)
end



# function check_update(um::SelfStarting, ch::C) where C <: UnivariateSeries
#     return true
# end