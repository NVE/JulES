"""
Generic fallbacks for AbstractJulESInput and AbstractJulESOutput
"""

"""
How price prognosis problems (ppp) are distributed on cores initially
"""
function get_ppp_dist(input::AbstractJulESInput) 
    cores = get_cores(input)
    scenarios = get_scenarios(input)

    N = length(cores)
    
    dist = Vector{ScenarioCore}(undef, length(scenarios))
    
    for (i, s) in enumerate(scenarios)
        j = (i - 1) % N + 1
        dist[i] = ScenarioCore(s, cores[j])
    end
    
    return dist
end

"""
How end value problems (evp) are distributed on cores initially
"""
function get_evp_dist(input::AbstractJulESInput) 
    cores = get_cores(input)
    scenarios = get_scenarios(input)
    subsystems = get_subsystems(input)

    N = length(cores)
    out = Vector{ScenarioSubsystemCore}(undef, length(scenarios) * length(subsystems))

    if N >= (length(scenarios) * length(subsystems))
        k = 0
        for sub in subsystems
            for scen in scenarios
                k += 1
                out[k] = ScenarioSubsystemCore(scen, sub, cores[k])
            end
        end
        return out
    end

    # TODO: Do better when length(scenarios) < N < length(scenarios) * length(subsystems)

    k = 0
    for sub in subsystems
        for (i, s) in enumerate(scenarios)
            k += 1
            j = (i - 1) % N + 1
            out[k] = ScenarioSubsystemCore(s, sub, cores[j])
        end
    end
    return out
end

function get_cp_core(input::AbstractJulESInput) 
    return first(get_cores(input))
end 

# TODO: complete these two (must be seen together)
function get_sp_dist(input::AbstractJulESInput) end
function get_mp_dist(input::AbstractJulESInput) end
