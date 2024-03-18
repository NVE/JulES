"""
Here we define some id/key-like types, that will make it easier to understand
the dimensions of the code, eg. understanding vectors and dicts in LocalDB
"""

mutable struct WeatherScenario <: AbstractScenario
    weather::Millisecond # as an offset from the simulationtime
    p_weather::Float64
    parentscenario::Int # index of parent scenario
end
getprobability(scen::WeatherScenario) = scen.p_weather

# Calculates end values with deterministic end value models
# TODO: Add name
struct EVSubsystem <: AbstractSubsystem
    priceareas::Vector{String}
    dataelements::Vector{Int}
    evduration::Millisecond
    stochduration::Millisecond
end
isevsubsystem(subsystem::EVSubsystem) = true
isspsubsystem(subsystem::EVSubsystem) = true

# Collects end value from price prognosis models
struct SPSubsystem <: AbstractSubsystem
    priceareas::Vector{String}
    dataelements::Vector{Int}
    stochduration::Millisecond
end
isevsubsystem(subsystem::SPSubsystem) = false
isspsubsystem(subsystem::SPSubsystem) = true

# Only subsystem model (no ppp, ev or clearing)
struct ExogenSubsystem <: AbstractSubsystem end

const ScenarioIx = Int
const SubsystemIx = Int
const TermName = String
const CommodityName = String
const CoreId = Int