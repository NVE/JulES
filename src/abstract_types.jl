# TODO: Interfaces
abstract type AbstractJulESInput end
abstract type AbstractJulESOutput end
abstract type AbstractScenario end

"""
In parts of the JulES algorithm the full problem is divided into subsystems that can be solved in parallel. 
AbstractSubsystem holds information about commodities, (priceareas) and (dataelements) in the subsystem.
It also decides what storage valuation models should be run for the subsystem, with what horizon and with what end value method.

Interface:
    get_commodities(s::AbstractSubsystem) -> commodities::Vector{CommodityName}
    is_subsystem_evp(s::AbstractSubsystem) -> Bool
    is_subsystem_stoch(s::AbstractSubsystem) -> Bool
    get_skipmed_impact(s::AbstractSubsystem) -> Bool
    get_endvaluemethod_xx(subsystem::EVPSubsystem) -> String

Optional:
    get_priceareas(subsystem::EVPSubsystem) -> Vector{String}
    get_dataelements(subsystem::EVPSubsystem) -> Vector{Int}
    get_horizonterm_xx(subsystem::EVPSubsystem) -> TermName
    get_duration_xx(subsystem::EVPSubsystem) -> Millisecond
"""
abstract type AbstractSubsystem end

"""
We use scenario modelling to consider uncertainty from all scenarios with only a few.
Scenarios can be reduced, altered and weighted with different methods

Interface:
    get_scenarios(::AbstractScenarioModellingMethod) -> Vector{AbstractScenario}
    choose_scenarios!(::AbstractScenarioModellingMethod, ::AbstractScenarioModellingMethod, ::TuLiPa.ProbTime, ::AbstractJulESInput) -> nothing
    perform_scenmod!(::AbstractScenarioModellingMethod, ::ScenarioIndex, ::Vector{Any}) -> nothing
    renumber_scenmodmethod!(scenmod::AbstractScenarioModellingMethod) -> nothing
"""
abstract type AbstractScenarioModellingMethod end

"""
Interface:
    TuLiPa.assemble!(m::AbstractInflowModel)
    TuLiPa.getid(m::AbstractInflowModel) -> TuLiPa.Id
    estimate_u0(m::AbstractInflowModel, t::ProbTime) -> u0::Vector{Float64}
    predict(m::AbstractInflowModel, u0::::Vector{Float64}, t::ProbTime) -> Q::Vector{Float64}
    get_numstates(m::AbstractInflowModel) -> Int  ( == length(u0))
    get_statename(::AbstractInflowModel, i::Int) -> String
    get_basin_area_m2(::AbstractInflowModel) -> Float64
"""
abstract type AbstractInflowModel end
const ABSTRACT_INFLOW_MODEL = "AbstractInflowModel"

# TODO: Interfaces
abstract type AbstractTwoStateIfmDataUpdater end
abstract type AbstractTwoStateIfmPredictor end

