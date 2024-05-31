abstract type AbstractJulESInput end
abstract type AbstractJulESOutput end
abstract type AbstractScenario end
abstract type AbstractSubsystem end
abstract type AbstractScenarioModellingMethod end

"""
Interface:
    update(m::AbstractInflowModel, t::ProbTime)

    estimate_u0(m::AbstractInflowModel) -> u0::Vector{Float64}

    predict(m::AbstractInflowModel, u0::::Vector{Float64}) -> (Q, u) 
        where Q::Vector{Float64} and u::Vector{Vector{Float64}}
"""
abstract type AbstractInflowModel end
const ABSTRACT_INFLOW_MODEL = "AbstractInflowModel"
