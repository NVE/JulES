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
    
    dist = Vector{Tuple{ScenarioIx, CoreId}}(undef, length(scenarios))
    
    for (i, s) in enumerate(scenarios)
        j = (i - 1) % N + 1
        dist[i] = (s, cores[j])
    end
    
    return dist
end

"""
How end value problems (evp) are distributed on cores initially
"""
function get_evp_dist(input::AbstractJulESInput) 
    cores = get_cores(input)
    scenarios = get_scenarios_evp(input)
    subsystems = get_subsystems(input)

    N = length(cores)
    S = length(scenarios)
    Y = length(subsystems)

    out = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, S*Y)

    if N >= S*Y
        k = 0
        for sub in subsystems
            for scen in scenarios
                k += 1
                out[k] = (scen, sub, cores[k])
            end
        end
        return out
    end

    # TODO: Do better when S < N < S*Y

    k = 0
    for sub in subsystems
        for (i, s) in enumerate(scenarios)
            k += 1
            j = (i - 1) % N + 1
            out[k] = (s, sub, cores[j])
        end
    end
    return out
end

function get_cp_core(input::AbstractJulESInput) 
    return first(get_cores(input))
end 

# TODO: Is it a better default to minimize total object count on each core? (Can lead to one big subsystem on one core and several subsystems on others). More complicated.
"""
Will try to distribute small and large subsystems evenly on cores.

E.g. if cores = [1, 2, 3] and subsystems = [3, 6, 1, 4, 5, 2]
where the list of subsystems is sorted from biggest subsystem to smallest
where size is defined as number of objects in the subsystem. 

Then, distribution on cores should be: 
- Core 1: Index 1 and 6 (that is, subsystem 3 and 2) (i.e. biggest and smallest)
- Core 2: Index 2 and 5 (that is, subsystem 6 and 5) (i.e. next biggest and next smallest)
- Core 3: Index 3 and 4 (that is, subsystem 1 and 4) (i.e. next next biggest and next next smallest)

This way, a core that get a relatively big system in one pass, 
will get a relatively small system in the next pass. 
This should ensure to balance load better than e.g. distributing 
subsystems on cores by random choice.

Scenario problems (sp) will be put on the same core as master problems (mp).
"""
function get_mp_sp_dist(input::AbstractJulESInput)
    cores = get_cores(input)
    subsystems = get_subsystem_ids_by_decending_size(input)

    mp_dist = _distribute_subsystems_by_size!(subsystems, cores)
    
    N = get_num_sp_scenarios(input)
    sp_dist = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, N*length(mp_dist))
    i = 0
    for scen in 1:N
        for (sub, core) in mp_dist
            i += 1
            sp_dist[i] = (scen, sub, core)
        end
    end

    return (mp_dist, sp_dist)
end

function _distribute_subsystems_by_size!(subsystems::Vector{SubsystemIx}, cores::Vector{CoreId})
    N = length(cores)

    dist = Tuple{SubsystemIx, CoreId}[]

    loop_forward = true
    while length(subsystems)
        K = length(subsystems)
        M = min(K, N)

        if loop_forward
            for i in 1:M
                push!(dist, (subsystems[i], core[i]))
            end
            if K > M
                subsystems = subsystems[(M+1):end]
            end
            loop_forward = false                        
        else
            for i in 1:M
                j = K - (i-1)
                push!(dist, (subsystems[j], core[i]))
            end
            if K > M
                subsystems = subsystems[1:(K-M)]
            end
            loop_forward = true
        end
    end

    return dist
end

function get_subsystem_ids_by_decending_size(input::AbstractJulESInput)
    subsystems = get_subsystems(input)
    subsystems = [(length(s), id) for (id, s) in subsystems]
    sort!(subsystems, rev=true)
    return [id for (n, id) in subsystems]
end
