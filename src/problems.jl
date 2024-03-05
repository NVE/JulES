

struct PP_Problem
    short
    med
    long
    short_prices
    med_prices
    long_prices
    short_nonstoragestates
    medendvaluesdict # TODO: stored here or in local db?
end

setstartstates!(p::PP_Problem, startstates)

struct EV_Problem
end

struct CP_Problem
    prob
    endstates
end


create_pp_problem(input::AbstractJulESInput, scenario) = nothing

function create_ev_problem(input::AbstractJulESInput, subsystem, scenario)
    
end

create_mp_problem(input::AbstractJulESInput, subsystem) = nothing
create_sp_problem(input::AbstractJulESInput, subsystem, scenario) = nothing
create_cp_problem(input::AbstractJulESInput) = nothing
