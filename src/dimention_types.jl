"""
Here we define some id/key-like types, that will make it easier to understand
the dimentions of the code, eg. understanding vectors and dicts in LocalDB
"""

"""
Main source of uncertainty in JulES is weather-related, as this is very importaint in
energy markets. Usually we will only use this uncertainty-dimention. However,
in order to support extentions of JulES, we allow for a second uncertainty-dimention,
which could e.g. be high/medium/low fuel prices, or high/medium/low availability of
capacity.
"""
mutable struct Scenario
    weather::Int
    p_weather::Float64
    other::Int
    p_other::Float64
end
probability(s::Scenario) = s.p_other * s.p_weather

struct Subsystem
    name::String
end

const Core = Int
const Term = String
const Commodity = String

struct ScenarioTermCommodity
    scenario::Scenario
    term::Term
    commodity::Commodity    # TODO: In conflict with TuLiPa.Commodity?
end

struct TermCommodity
    term::Term
    commodity::Commodity    # TODO: In conflict with TuLiPa.Commodity?
end

struct ScenarioSubsystem
    scenario::Scenario
    subsystem::Subsystem
end

# These are used in the [problem]_dist
# slots in the LocalDB (eg. pp_dist)
# Dynamic load balancer may modify the core slot
mutable struct ScenarioCore
    const scenario::Scenario
    core::Core
end

mutable struct ScenarioSubsystemCore
    const scenario::Scenario
    const subsystem::Subsystem
    core::Core
end

mutable struct SubsystemCore
    const subsystem::Subsystem
    core::Core
end