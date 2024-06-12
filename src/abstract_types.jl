abstract type AbstractJulESInput end
abstract type AbstractJulESOutput end
abstract type AbstractScenario end
abstract type AbstractSubsystem end
abstract type AbstractScenarioModellingMethod end

"""
Interface:
    TuLiPa.assemble!(m::AbstractInflowModel)
    TuLiPa.getid(m::AbstractInflowModel) -> TuLiPa.Id
    estimate_u0(m::AbstractInflowModel, t::ProbTime) -> u0::Vector{Float64}
    predict(m::AbstractInflowModel, u0::::Vector{Float64}, t::ProbTime) -> Q::Vector{Float64}
"""
abstract type AbstractInflowModel end
const ABSTRACT_INFLOW_MODEL = "AbstractInflowModel"

# TODO: Interfaces
abstract type AbstractTwoStateIfmDataUpdater end
abstract type AbstractTwoStateIfmPredictor end

