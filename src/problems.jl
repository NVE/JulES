

struct PricePrognosisProblem
    short
    med
    long
    short_prices
    med_prices
    long_prices
    short_nonstoragestates
end

setstartstates!(p::PricePrognosisProblem, startstates)

struct EndValueProblem
end

struct ScenarioProblem
end

struct MasterProblem
end

struct ClearingProblem
    prob
    endstates
end


create_ppp(input::AbstractJulESInput, scenario) = nothing

function create_evp(input::AbstractJulESInput, subsystem, scenario)
    
end

create_mp(input::AbstractJulESInput, subsystem) = nothing
create_sp(input::AbstractJulESInput, subsystem, scenario) = nothing
create_cp(input::AbstractJulESInput) = nothing
