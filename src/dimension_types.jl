"""
Here we define some id/key-like types, that will make it easier to understand
the dimensions of the code, eg. understanding vectors and dicts in LocalDB
"""

mutable struct WeatherScenario <: AbstractScenario
    weather::Millisecond # as an offset from the simulationtime
    p_weather::Float64
    parentscenario::Int # index of parent scenario
end

const ScenarioIx = Int
const SubsystemIx = Int
const TermName = String
const CommodityName = String
const CoreId = Int
