"""
Definition of default input and output types
"""


struct DefaultJulESInput <: AbstractJulESInput
    # TODO: complete

    function DefaultJulESInput(dataset, config)
        # time

        # scenariogeneration? bør bo på 

        horizons = gethorizons(config)

        subsystems = getsubsystems(dataset, config, horizons)::Dict{String, Vector{DataElement}}

        # TODO: use throw-away-dummy-modelobjs for data proc

        return new()
    end
end

struct DefaultJulESOutput <: AbstractJulESOutput
    # TODO: complete
end

# TODO: complete
get_cores(input) = nothing
get_horizons(input) = nothing
get_simulation_period(input) = nothing
get_startstates_ppp(input) = nothing



# Should live here and not in slot in PricePrognosisProblem?
function get_medendvaluesdict(input)
end

