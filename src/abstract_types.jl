abstract type AbstractJulESInput end
abstract type AbstractJulESOutput end
abstract type AbstractScenario end
abstract type AbstractSubsystem end
abstract type AbstractScenarioModellingMethod end

# TODO: Update code to return and store S for output reporting
"""
Interface: 
estimate_S0(m::AbstractInflowModel, t::ProbTime) -> S0::Vector{Float64}
predict(m::AbstractInflowModel, S0::::Vector{Float64}, t::ProbTime) -> (Q, S) where Q::Vector{Float64} and S::Vector{Vector{Float64}}
"""
abstract type AbstractInflowModel end
const ABSTRACT_INFLOW_MODEL = "AbstractInflowModel"
