

struct PP_Problem
end

struct EV_Problem
end

create_pp_problem(input::AbstractJulESInput, scenario) = nothing

function create_ev_problem(input::AbstractJulESInput, subsystem, scenario)
    
end

create_mp_problem(input::AbstractJulESInput, subsystem) = nothing
create_sp_problem(input::AbstractJulESInput, subsystem, scenario) = nothing
create_cp_problem(input::AbstractJulESInput) = nothing
