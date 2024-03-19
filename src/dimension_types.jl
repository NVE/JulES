"""
Here we define some id/key-like types, that will make it easier to understand
the dimensions of the code, eg. understanding vectors and dicts in LocalDB
"""

mutable struct WeatherScenario <: AbstractScenario
    weather::Millisecond # as an offset from the simulationtime
    p_weather::Float64
    parentscenario::Int # index of parent scenario
end
get_probability(scen::WeatherScenario) = scen.p_weather

# Calculates end values with deterministic end value models
# TODO: Add name
struct EVPSubsystem <: AbstractSubsystem
    priceareas::Vector{String}
    dataelements::Vector{Int}
    evpduration::Millisecond
    stochduration::Millisecond
end
get_priceareas(subsystem::EVPSubsystem) = subsystem.priceareas
get_dataelements(subsystem::EVPSubsystem) = subsystem.dataelements
get_evpduration(subsystem::EVPSubsystem) = subsystem.evpduration
get_stochduration(subsystem::EVPSubsystem) = subsystem.stochduration
is_subsystem_evp(subsystem::EVPSubsystem) = true
is_subsystem_stoch(subsystem::EVPSubsystem) = true

# Collects end value from price prognosis models
struct StochSubsystem <: AbstractSubsystem
    priceareas::Vector{String}
    dataelements::Vector{Int}
    stochduration::Millisecond
end
get_priceareas(subsystem::StochSubsystem) = subsystem.priceareas
get_dataelements(subsystem::StochSubsystem) = subsystem.dataelements
get_stochduration(subsystem::StochSubsystem) = subsystem.stochduration
is_subsystem_evp(subsystem::StochSubsystem) = false
is_subsystem_stoch(subsystem::StochSubsystem) = true

# Only subsystem model (no ppp, ev or clearing)
struct ExogenSubsystem <: AbstractSubsystem end
is_subsystem_evp(subsystem::ExogenSubsystem) = false
is_subsystem_stoch(subsystem::ExogenSubsystem) = true

const ScenarioIx = Int
const SubsystemIx = Int
const TermName = String
const CommodityName = String
const CoreId = Int