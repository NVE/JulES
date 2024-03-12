

struct PricePrognosisProblem
    short
    med
    long
    short_prices
    med_prices
    long_prices
    short_nonstoragestates
end

function setstartstates!(p::PricePrognosisProblem, startstates)
end

# TODO: Support different terms? (med or long)
struct EndValueProblem
    prob
end

# TODO: Support different terms? (med or long)
struct ScenarioProblem
    prob
end

struct MasterProblem
    prob
end

struct ClearingProblem
    prob
    endstates
end


function create_ppp(input::AbstractJulESInput, scenario)
end

function create_evp(input::AbstractJulESInput, scenario, subsystem)
end

function create_mp(input::AbstractJulESInput, subsystem) 
end

function create_sp(input::AbstractJulESInput, scenario, subsystem)
end

function create_cp(input::AbstractJulESInput)
end
