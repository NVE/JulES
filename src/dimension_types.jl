"""
Here we define some id/key-like types, that will make it easier to understand
the dimensions of the code, eg. understanding vectors and dicts in LocalDB
"""

# Alias types to make code more readable
const ScenarioIx = Int
const SubsystemIx = Int
const TermName = String
const CommodityName = String
const CoreId = Int

"""
Main source of uncertainty in JulES is weather-related, as this is very importaint in
energy markets. Usually we will only use this uncertainty-dimension. However,
in order to support extentions of JulES, we allow for a second uncertainty-dimension,
which could e.g. be high/medium/low fuel prices, or high/medium/low availability of
capacity.
"""
mutable struct Scenario
    weather::Millisecond # as an offset from the simulationtime
    p_weather::Float64
    other::Dict{String, Tuple{Any,Float64}} # add possibility for other dimensions
    parentscenario::Int # index of parent scenario
end
function getprobability(s::Scenario)
    probability = s.p_weather
    for (key, value) in other
        probability *= last(value)
    end
    return probability
end
