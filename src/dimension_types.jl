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
    commodities::Vector{CommodityName}
    priceareas::Vector{String}
    dataelements::Vector{Int}
    duration_evp::Millisecond
    duration_stoch::Millisecond
end
get_commodities(subsystem::EVPSubsystem) = subsystem.commodities
get_priceareas(subsystem::EVPSubsystem) = subsystem.priceareas
get_dataelements(subsystem::EVPSubsystem) = subsystem.dataelements
get_duration_evp(subsystem::EVPSubsystem) = subsystem.duration_evp
get_duration_stoch(subsystem::EVPSubsystem) = subsystem.duration_stoch
is_subsystem_evp(subsystem::EVPSubsystem) = true
is_subsystem_stoch(subsystem::EVPSubsystem) = true

# Collects end value from price prognosis models
struct StochSubsystem <: AbstractSubsystem
    commodities::Vector{CommodityName}
    priceareas::Vector{String}
    dataelements::Vector{Int}
    duration_stoch::Millisecond
end
get_commodities(subsystem::StochSubsystem) = subsystem.commodities
get_priceareas(subsystem::StochSubsystem) = subsystem.priceareas
get_dataelements(subsystem::StochSubsystem) = subsystem.dataelements
get_duration_stoch(subsystem::StochSubsystem) = subsystem.duration_stoch
is_subsystem_evp(subsystem::StochSubsystem) = false
is_subsystem_stoch(subsystem::StochSubsystem) = true

# Only subsystem model (no ppp, ev or clearing)
struct ExogenSubsystem <: AbstractSubsystem 
    commodities::Vector{CommodityName}
end
get_commodities(subsystem::ExogenSubsystem) = subsystem.commodities
is_subsystem_evp(subsystem::ExogenSubsystem) = false
is_subsystem_stoch(subsystem::ExogenSubsystem) = true

const ScenarioIx = Int
const SubsystemIx = Int
const TermName = String
const CommodityName = String
const CoreId = Int