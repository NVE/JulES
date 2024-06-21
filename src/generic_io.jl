"""
Generic fallbacks for AbstractJulESInput and AbstractJulESOutput
"""

"""
How price prognosis problems (ppp) are distributed on cores initially
"""
function get_dist_ppp(input::AbstractJulESInput)
    cores = get_cores(input)
    N = length(cores)
    S = get_numscen_sim(input)
    
    dist = Vector{Tuple{ScenarioIx, CoreId}}(undef, S)
    
    for s in 1:S
        j = (s - 1) % N + 1
        #println("j in get_dist_ppp is $j")
        dist[s] = (s, cores[j])
    end
    
    return dist
end

"""
How end value problems (evp) are distributed on cores initially
"""
function get_dist_evp(input::AbstractJulESInput, subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}) 
    cores = get_cores(input)

    N = length(cores) #number of cores
    S = get_numscen_sim(input) #number of scenarios
    Y = length(subsystems)

    # println("Y is $Y")
    # println("N is $N")
    # println("S is $S")


    out = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, S*Y)

    if N >= S*Y
        k = 0
        for (subix, sub) in subsystems
            for scenix in 1:S
                k += 1
                out[k] = (scenix, subix, cores[k])
                #println("scenix is $scenix , k is $k")
            end
        end
        return out
    end

    # TODO: Do better when S < N < S*Y

    k = 0
    for (subix, sub) in subsystems
        for scenix in 1:S
            k += 1
            j = (scenix - 1) % N + 1
            #println("j is $j")
            out[k] = (scenix, subix, cores[j])
            #println("N < S*Y, scenix is $scenix , k is $k")
        end
    end
    return out
end

function get_core_cp(input::AbstractJulESInput) 
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
function get_dist_stoch(input::AbstractJulESInput, subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
    cores = get_cores(input)
    subsystems_desc = get_subsystem_ids_by_decending_size(subsystems)
    
    distribution_method = get_distribution_method(input)
   # println(method)
    #distribution_method = config[method]["subsystems"]["distribution_function"]

    if distribution_method == "random"
        dist_mp = _distribute_subsystems_randomly!(subsystems_desc, cores)
    elseif distribution_method == "by_size"
        dist_mp = _distribute_subsystems_by_size!(subsystems_desc, cores)
    elseif distribution_method == "greedy"
        _distribute_subsystems_elements_smarter!(subsystems, cores)
    end
    
    N = get_numscen_stoch(input)
    dist_sp = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, N*length(dist_mp))
    i = 0
    for scen in 1:N
        for (sub, core) in dist_mp
            i += 1
            dist_sp[i] = (scen, sub, core)
        end
    end

    return (dist_mp, dist_sp)
end



        

function _distribute_subsystems_randomly!(subsystems::Vector{SubsystemIx}, cores::Vector{CoreId})
   
    #println("subsystems: $subsystems")
    dist = Tuple{SubsystemIx, CoreId}[]

    for sub in subsystems
        # Randomly select a core
        core = rand(cores)

        #println("sub is $sub, core is $core")
        # Assign the subsystem to the randomly selected core
        push!(dist, (sub, core))
    end

    return dist
end



function _distribute_subsystems_by_size!(subsystems::Vector{SubsystemIx}, cores::Vector{CoreId})
    N = length(cores)

    dist = Tuple{SubsystemIx, CoreId}[]

    loop_forward = true
    while length(subsystems) > 0
        K = length(subsystems)
        M = min(K, N)
        # println("M is $M")
        # println("K is $K")
        if loop_forward
            for i in 1:M
                #println("i is $i")

                push!(dist, (subsystems[i], cores[i]))
            end
            subsystems = subsystems[(M+1):end]
            loop_forward = false                        
        else
            for i in 1:M
                j = K - (i-1)
                push!(dist, (subsystems[j], cores[i]))
            end
            subsystems = subsystems[1:(K-M)]
            loop_forward = true
        end
    end

    return dist
end

 
# elements_all = get_elements(input)
# elements_ix_in_subsystem = get_dataelements(s)
# elements_in_subsystem = elements[elements_ix]
# for element in elements_in_subsystem
#     if element.conceptname == "Storage"

function get_subsystem_ids_by_decending_size(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
    subsystems = [(length(get_dataelements(s)), ix) for (ix, s) in subsystems]
    sort!(subsystems, rev=true)
    return [ix for (n, ix) in subsystems]
end

function _distribute_subsystems_elements_smarter!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
   

    num_cores = length(cores)

    dist = Tuple{SubsystemIx, CoreId}[]
    # Initialize an array to hold the load (number of data elements) for each core
    core_loads = fill(0, num_cores)
    println("1")
    sorted_elements = get_sorted_subsystem_number_of_elements(subsystems)
    println(typeof(sorted_elements))
    println(sorted_elements)
    for (subsystem_index, num_elements) in sorted_elements
        # Find the core with the least load
        min_load_core_index = argmin(core_loads)
        # Assign the subsystem to this core
        push!(dist,(subsystem_index,min_load_core_index))
        # Update the load for this core
        core_loads[min_load_core_index] += num_elements
    end
    return dist
end



#function to get a tuple(subsystem index, number of data elements) that is sorted from highest to lowest number of data elements in the subsystem
#function get_sorted_subsystem_number_of_elements(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
function get_sorted_subsystem_number_of_elements(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
   
    dataelements_in_each_subsystem = Tuple{SubsystemIx, Int}[]
 
    for (ix, s) in subsystems
        
        data_elements_in_sub= get_dataelements(s)
        number_data_elements = length(data_elements_in_sub)
 
        # Create a tuple (data_element_index, number_data_elements)
        data_tuple = (ix, number_data_elements)
   
        # Add the tuple to the dataelements array
        push!(dataelements_in_each_subsystem, data_tuple)
    end
   
   
    #want to sort the tuple by the number of data_elements in each subsystem:
    sorted_dataelements = sort(dataelements_in_each_subsystem, by = x -> x[2], rev=true)
 
    return sorted_dataelements
end


function distribute_subsystems_cores(sorted_subsystems::Vector{Tuple{Int, Int}}, num_cores::Int)
    # Initialize an array to hold the load (number of data elements) for each core
    core_loads = fill(0, num_cores)
    # Initialize an array to hold the subsystems assigned to each core
    core_assignments = [Vector{Int}() for _ in 1:num_cores]
 
    for (subsystem_index, num_elements) in sorted_subsystems
        # Find the core with the least load
        min_load_core = argmin(core_loads)
        # Assign the subsystem to this core
        push!(core_assignments[min_load_core], subsystem_index)
        # Update the load for this core
        core_loads[min_load_core] += num_elements
    end
 
    return core_assignments
end

"""
Find which subsystems should have evp and stoch problems
"""
function get_subsystems_evp(allsubsystems::Vector{AbstractSubsystem})
    subsystems = Tuple{SubsystemIx, AbstractSubsystem}[]
    for (i, subsystem) in enumerate(allsubsystems)
        if is_subsystem_evp(subsystem)
            push!(subsystems, (i, subsystem))
        end
    end
    return subsystems
end
function get_subsystems_stoch(allsubsystems::Vector{AbstractSubsystem})
    subsystems = Tuple{SubsystemIx, AbstractSubsystem}[]
    for (i, subsystem) in enumerate(allsubsystems)
        if is_subsystem_stoch(subsystem)
            push!(subsystems, (i, subsystem))
        end
    end
    return subsystems
end