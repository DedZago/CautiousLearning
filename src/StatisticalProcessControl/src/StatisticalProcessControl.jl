module StatisticalProcessControl

using Parameters

export EWMA,
        AEWMA,
        get_params,
        get_limits,
        get_value,
        update_value,
        update_series,
        check_IC,
        check_OC,
        apply_series,
        CautiousLearning,
        FixedParameter,
        SelfStarting,
        check_update,
        get_warning_limit,
        ControlChart,
        DataSimulation

abstract type AbstractSeries end

abstract type UnivariateSeries <: AbstractSeries end

abstract type SingleSidedUnivariateSeries <: UnivariateSeries end
abstract type DoubleSidedUnivariateSeries <: UnivariateSeries end
# Write your package code here.

get_params(ch::AbstractSeries) = @NamedTuple{}
get_limits(ch::AbstractSeries) = ch.L
get_value(ch::AbstractSeries) = ch.C
update_series(ch::AbstractSeries) = ch

abstract type AbstractUpdate end

include("controlseries.jl")
include("updatemechanism.jl")


end
