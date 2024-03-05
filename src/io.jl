


struct DefaultJulESInput <: AbstractJulESInput
    pp_dist::Dict
    ev_dist::Dict
    sp_dist::Dict
    # TODO: complete

    function DefaultJulESInput(dataset, config)
        # time

        # scenariogeneration? bør bo på 

        horizons = gethorizons(config)

        subsystems = getsubsystems(dataset, config, horizons)::Dict{String, Vector{DataElement}}

        # TODO: use throw-away-dummy-modelobjs for data proc
        pp_dist = distribute_pp(dataset, config)
        ev_dist = distribute_ev(dataset, config)

        return new(pp_dist, ev_dist)
    end
end

struct DefaultJulESOutput <: AbstractJulESOutput
    # TODO: complete
end

get_cores(input) = nothing
get_horizons(input) = nothing
get_simulation_period(input) = nothing
get_startstates_pp(input) = nothing


"""
How price prognosis (pp) problems are distributed on cores initially
"""
function get_pp_dist(input::DefaultJulESInput) end

function get_ev_dist(input::DefaultJulESInput) end

function get_sp_dist(input::DefaultJulESInput) end

function get_mp_dist(input::DefaultJulESInput) end

function get_cp_core(input::DefaultJulESInput) end 

