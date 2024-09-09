"""
Here we define some id/key-like types, that will make it easier to understand
the dimensions of the code, eg. understanding vectors and dicts in LocalDB
"""

# TODO: Also have type in alias
const ScenarioIx = Int
const SubsystemIx = Int
const TermName = String
const ShortTermName = "short"
const MedTermName = "med"
const LongTermName = "long"
const ClearingTermName = "clearing"
const MasterTermName = "master"
const SubTermName = "sub"
const CommodityName = String
const CoreId = Int
const MainTiming = "maintiming"
const StorageValues = "storagevalues"

mutable struct WeatherScenario <: AbstractScenario
    weatheroffset::Millisecond # as an offset from the simulationtime
    p_weather::Float64
    parentscenario::Int # index of parent scenario
end
get_probability(scen::WeatherScenario) = scen.p_weather
set_probability!(scen::WeatherScenario, value::Float64) = scen.p_weather = value

# Calculates end values with deterministic end value models
# TODO: Add name
struct EVPSubsystem <: AbstractSubsystem
    commodities::Vector{CommodityName}
    priceareas::Vector{String}
    dataelements::Vector{Int}
    horizonterm_evp::TermName
    duration_evp::Millisecond
    horizonterm_stoch::TermName
    duration_stoch::Millisecond
    endvaluemethod_evp::String
    skipmed_impact::Bool
end
get_commodities(subsystem::EVPSubsystem) = subsystem.commodities
get_priceareas(subsystem::EVPSubsystem) = subsystem.priceareas
get_dataelements(subsystem::EVPSubsystem) = subsystem.dataelements
get_horizonterm_evp(subsystem::EVPSubsystem) = subsystem.horizonterm_evp
get_duration_evp(subsystem::EVPSubsystem) = subsystem.duration_evp
get_horizonterm_stoch(subsystem::EVPSubsystem) = subsystem.horizonterm_stoch
get_duration_stoch(subsystem::EVPSubsystem) = subsystem.duration_stoch
is_subsystem_evp(subsystem::EVPSubsystem) = true
is_subsystem_stoch(subsystem::EVPSubsystem) = true
get_endvaluemethod_evp(subsystem::EVPSubsystem) = subsystem.endvaluemethod_evp
get_endvaluemethod_sp(subsystem::EVPSubsystem) = "evp"
get_skipmed_impact(subsystem::EVPSubsystem) = subsystem.skipmed_impact

# Collects end value from price prognosis models
struct StochSubsystem <: AbstractSubsystem
    commodities::Vector{CommodityName}
    priceareas::Vector{String}
    dataelements::Vector{Int}
    horizonterm_stoch::TermName
    duration_stoch::Millisecond
    endvaluemethod_sp::String
    skipmed_impact::Bool
end
get_commodities(subsystem::StochSubsystem) = subsystem.commodities
get_priceareas(subsystem::StochSubsystem) = subsystem.priceareas
get_dataelements(subsystem::StochSubsystem) = subsystem.dataelements
get_horizonterm_stoch(subsystem::StochSubsystem) = subsystem.horizonterm_stoch
get_duration_stoch(subsystem::StochSubsystem) = subsystem.duration_stoch
is_subsystem_evp(subsystem::StochSubsystem) = false
is_subsystem_stoch(subsystem::StochSubsystem) = true
get_endvaluemethod_sp(subsystem::StochSubsystem) = subsystem.endvaluemethod_sp
get_skipmed_impact(subsystem::StochSubsystem) = subsystem.skipmed_impact

# Only subsystem model (no ppp, evp or cp)
struct ExogenSubsystem <: AbstractSubsystem 
    commodities::Vector{CommodityName}
    endvaluemethod_sp::String
end
get_commodities(subsystem::ExogenSubsystem) = subsystem.commodities
is_subsystem_evp(subsystem::ExogenSubsystem) = false
is_subsystem_stoch(subsystem::ExogenSubsystem) = true
get_endvaluemethod_sp(subsystem::ExogenSubsystem) = subsystem.endvaluemethod_sp
get_skipmed_impact(subsystem::ExogenSubsystem) = false