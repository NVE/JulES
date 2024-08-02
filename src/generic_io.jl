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

function get_dist_stoch(input::AbstractJulESInput, subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
    
    cores = get_cores(input)
    #Distributing the master problems/subsystems
    
    
    
    distribution_method_mp = get_distribution_method_mp(input)
    default = "bysize"

    valid_methods_mp = ["randdumb", "random", "bysize", "greedy", "storage", "sizepairing", "advanced"]

    # Check if distribution_method is valid
    if !(distribution_method_mp in valid_methods_mp)
        println("distribution method $distribution_method_mp is not valid. Using $default")
        distribution_method_mp = default
    end

    
    if distribution_method_mp == "randdumb"
        dist_mp = _distribute_subsystems_randdumb!(subsystems, cores)
    elseif distribution_method_mp == "random"
        dist_mp = _distribute_subsystems_random!(subsystems, cores) 
    elseif distribution_method_mp == "bysize"
        dist_mp = _distribute_subsystems_by_size!(subsystems, cores)
    elseif distribution_method_mp == "greedy"
        dist_mp = _distribute_subsystems_elements_greedy!(subsystems, cores)
    elseif distribution_method_mp == "storage"
        dist_mp = _distribute_subsystems_storage_greedy!(input, subsystems, cores)
    elseif distribution_method_mp == "sizepairing"
        dist_mp = _distribute_subsystems_big_small!(subsystems, cores)
    elseif distribution_method_mp == "advanced"
        dist_mp = _distribute_subsystems_advanced!(subsystems, cores)
        dist_sp =  _distribute_subscenarios_advanced!(dist_mp, cores, input)
        #Returning here to avoid the rest of the code, since dist_mp is already decided
        return (dist_mp, dist_sp)
    end
    

    
    #Distributing scenarioproblems sp
    distribution_method_sp = get_distribution_method_sp(input)
    default = "withmp"
    valid_methods_sp = ["withmp", "greedy"]

     # Check if distribution_method_sp is valid
     if !(distribution_method_sp in valid_methods_sp)
        println("distribution method $distribution_method_sp is not valid. Using $default")
        distribution_method_sp = default
    end

    if distribution_method_sp == "withmp"
        dist_sp = _distribute_sp_with_mp!(input, dist_mp)
    elseif distribution_method_sp == "greedy"
        #Using _get_core_load!() to get the core loads
        core_loads = _get_core_load!(input, subsystems, dist_mp)
        #using _distribute_scenarios_greedy!() to find the distribution of scenario problems
        dist_sp = _distribute_scenarios_greedy!(input, subsystems, dist_mp, core_loads)
    end
    return (dist_mp, dist_sp)
end

"""Function to get the distribution of data elements on the different cores after distributing master problems"""
function _get_core_load!(input::AbstractJulESInput, subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, dist_mp::Vector{Tuple{SubsystemIx, CoreId}})
    cores = get_cores(input)

    # Initialize a dictionary (CoreId => load) where load is initially zero
    core_loads = Dict{CoreId, Int}()
    for core in cores
        core_loads[core] = 0
    end

    
    # Create a dictionary to map SubsystemIx to AbstractSubsystem
    subsystem_dict = Dict{SubsystemIx, AbstractSubsystem}()
    for (ix, s) in subsystems
        subsystem_dict[ix] = s
    end


    for (sub_ix, core) in dist_mp
        s = subsystem_dict[sub_ix]
        num_elements = length(get_dataelements(s))

        # Update the load for the core directly in the dictionary
        if haskey(core_loads, core)
            core_loads[core] += num_elements
            
        else
            println("Warning: Core $core not found in core_loads dictionary. Skipping...")
        end

        
    end
    return core_loads
end

#withmp
"""The original way to distribute senario problems. Each sp is on the same core as its master problem"""
function _distribute_sp_with_mp!(input::AbstractJulESInput, dist_mp::Vector{Tuple{SubsystemIx, CoreId}})
    N = get_numscen_stoch(input)
        dist_sp = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, N*length(dist_mp))
        i = 0
        for scen in 1:N
            for (sub, core) in dist_mp
                i += 1
                dist_sp[i] = (scen, sub, core)
            end
        end
    return dist_sp
end

#greedy
"""Function to distribute the scenarios on the core with least data elements"""
function _distribute_scenarios_greedy!(input::AbstractJulESInput, subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, dist_mp::Vector{Tuple{SubsystemIx, CoreId}}, core_loads::Dict{CoreId, Int})
    elements = get_elements(input)
    
    N = get_numscen_stoch(input)

    dist_sp = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, N*length(dist_mp))
    
    #vektor med tupler med subsystem id, subsystem og antall dataelementer
    subsystem_elements = Tuple{SubsystemIx, Int}[]

    for (ix, s) in subsystems
        push!(subsystem_elements, (ix, length(get_dataelements(s))))
    end
   
    #sorting subsystem_elements from largest  number of data elements to lowest
    sort!(subsystem_elements, by=x -> x[2], rev=true)

    #function to distribute subsystems on cores greedy
    function assign_scenarios!(subsystem_elements, core_loads, dist_sp)
        i=1  # Index to keep track of the position in dist_sp
        for (ix, elements) in subsystem_elements
            for scen in 1:N 
            
                # Find the core with the minimum load
            
                min_load_core = findmin(core_loads)[2]
                dist_sp[i] = (scen, ix, min_load_core)
                core_loads[min_load_core] += elements
                
                i += 1
            end
        end
    end

    assign_scenarios!(subsystem_elements, core_loads, dist_sp)
    

    return dist_sp
end

#advanced MP
"""Making the largest MP get the two first cores by assigning largest MP to core 1 and "reserving" core two for its SP"""
function _distribute_subsystems_advanced!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
    
    dist = Tuple{SubsystemIx, CoreId}[]
    sorted_subsystems = get_subsystem_ids_by_decending_size(subsystems) #største subsystems først
    
    # Make a copy of the cores array to avoid modifying the original
    cores_copy = copy(cores)
    subsystems_copy = copy(subsystems)
    
    #the first/largest masterproblem (subsystem) gets the two first cores
    # Check if there are at least two cores
        if length(cores) >= 2
        
        push!(dist, (sorted_subsystems[1], cores[1]))
        
        
        # Remove the first subsystem and the first two cores
        deleteat!(subsystems_copy, sorted_subsystems[1])
        popfirst!(cores_copy)
        popfirst!(cores_copy)
    end

    #the rest of the masterproblems are distributed by_size
    dist_rest = _distribute_subsystems_by_size!(subsystems_copy, cores_copy)

    # Append dist_rest to dist
    append!(dist, dist_rest)

    return dist
end

#advanced SP
"""Function to split the scenario problems of the largest MP on the first two cores, and then distribute the rest of SP as before"""
function _distribute_subscenarios_advanced!(dist_mp::Vector{Tuple{SubsystemIx, CoreId}}, cores::Vector{CoreId}, input::AbstractJulESInput)
    N = get_numscen_stoch(input)
    dist_sp = Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}(undef, N * length(dist_mp))
    i = 0

    # Get the first two cores from the cores input
    first_core1 = cores[1]
    first_core2 = cores[2]
    
    # Find the largest MP problem
    first_sub = dist_mp[1][1]

    # Iterate through the first half of SP
    for scen in 1:div(N, 2)
        i += 1
        dist_sp[i] = (scen, first_sub, first_core1)
    end

    # The second half of scenarios get the next core
    for scen in div(N, 2)+1:N
        i += 1
        dist_sp[i] = (scen, first_sub, first_core2)
    end

    # Distribute the rest of the scenarios as before
    for j in 2:length(dist_mp)
        (sub, core) = dist_mp[j]
        for scen in 1:N
            i += 1
            dist_sp[i] = (scen, sub, core)
        end
    end

    # Remove uninitialized items at the end
    dist_sp = dist_sp[1:i]

    return dist_sp
end

#sizepairing
"""Function to sort all subsystems from biggest to smallest in a list, distribute smallest and biggest subsystem on the same core, removing them from the list and repeating until all subsystems are distributed on cores"""
function _distribute_subsystems_big_small!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
    sorted_subsystems = get_subsystem_ids_by_decending_size(subsystems)
    
    N = length(cores)
    dist = Tuple{SubsystemIx, CoreId}[]

    while length(sorted_subsystems) > 0
        K = length(sorted_subsystems)
        M = min(K, N)

        for i in 1:M
            if length(sorted_subsystems) == 0
                break
            end
            
            # Assign the smallest subsystem to the current core
            push!(dist, (sorted_subsystems[1], cores[i]))
            deleteat!(sorted_subsystems, 1)
            
            if length(sorted_subsystems) == 0
                break
            end
            
            # Assign the largest subsystem to the same core
            push!(dist, (sorted_subsystems[end], cores[i]))
            deleteat!(sorted_subsystems, length(sorted_subsystems))
        end
    end

    return dist
end
     
#randdumb
""" Function to distribute the subsystems on completely random cores"""   
function _distribute_subsystems_randdumb!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
   
    
    dist = Tuple{SubsystemIx, CoreId}[]

    for (ix, s) in subsystems
        # Randomly select a core
        core = rand(cores)

        
        # Assign the subsystem to the randomly selected core
        push!(dist, (ix, core))
    end

    return dist
end

#random
"""Function to distribute the subsystems on random cores, but on different cores until all cores have subproblems assigned to them, then again until all cores have 2 subproblems assigned to them, and so on"""
function _distribute_subsystems_random!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
   
    
    dist = Tuple{SubsystemIx, CoreId}[]

    corelist = copy(cores)

    
    for (ix, s) in subsystems
        if isempty(corelist)
            corelist = copy(cores)
        end
        core = rand(corelist)
        filter!(x -> x != core, corelist)  # Remove the selected core from the list
        print(corelist)
        # Assign the subsystem to the randomly selected core
        push!(dist, (ix, core))
        
    end

    return dist
end

#bysize
"""
Function that distributes the subsystems by size. 
"""
function _distribute_subsystems_by_size!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
    N = length(cores)

    subsystems_desc = get_subsystem_ids_by_decending_size(subsystems)

    dist = Tuple{SubsystemIx, CoreId}[]

    loop_forward = true
    while length(subsystems_desc) > 0
        K = length(subsystems_desc)
        M = min(K, N)
        
        if loop_forward
            for i in 1:M
                

                push!(dist, (subsystems_desc[i], cores[i]))
            end
            subsystems_desc = subsystems_desc[(M+1):end]
            loop_forward = false                        
        else
            for i in 1:M
                j = K - (i-1)
                push!(dist, (subsystems_desc[j], cores[i]))
            end
            subsystems_desc = subsystems_desc[1:(K-M)]
            loop_forward = true
        end
    end

    return dist
end

#storage
"""Function that first distributes the subsystems with the highest number of storage elements on the cores with the lowest number of dataelements, then the remaining subsystems are distributed on the core with fewest dataelements one by one. """
function _distribute_subsystems_storage_greedy!(input::AbstractJulESInput, subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
    elements = get_elements(input)

    dist = Tuple{SubsystemIx, CoreId}[]

    num_cores = length(cores)

    has_storage = Tuple{SubsystemIx, AbstractSubsystem, Int}[]
    not_storage = Tuple{SubsystemIx, AbstractSubsystem, Int}[]

    # Initialize an array to hold the load (number of data elements) for each core
    core_loads = fill(0, num_cores)

    for (ix, s) in subsystems
        
        
        elements_ix_in_subsystem = get_dataelements(s)
        elements_in_subsystem = elements[elements_ix_in_subsystem]

        storage = false
        num_storage = 0

        #iterating through elements in the subsystem to check number of Storage elements. Adding subsystem and number of Storage elements into has_storage and not_storage
        for element in elements_in_subsystem
            
            if element.conceptname == "Storage"
                num_storage += 1
                storage = true
            end
        end
        if storage == true
            push!(has_storage, (ix, s, num_storage))
        else
            push!(not_storage, (ix, s, length(get_dataelements(s))))
        end
    end
    #sorting has_storage from largest number of Storage elements to lowest
    sort!(has_storage, by=x -> x[3], rev=true)
    #sorting not_storage from largest number of data elements to lowest
    sort!(not_storage, by=x -> x[3], rev=true)

    #function to distribute subsystems on cores greedy
    function assign_subsystems!(subsystems, core_loads, dist)
        for (ix, s, num) in subsystems
            min_load_core = argmin(core_loads)
            core_id = cores[min_load_core]  # Get the actual CoreId
            push!(dist, (ix, core_id))
            core_loads[min_load_core] += num
        end
    end

    assign_subsystems!(has_storage, core_loads, dist)
    assign_subsystems!(not_storage, core_loads, dist)

    return dist
end


"""Function to sort subsystems by decending number of dataelements """
function get_subsystem_ids_by_decending_size(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
    subsystems = [(length(get_dataelements(s)), ix) for (ix, s) in subsystems]
    sort!(subsystems, rev=true)
    return [ix for (n, ix) in subsystems]
end

#greedy
""" Function to distribute subsystems greedy by number of dataelements. First sorting the subsystems in decreasing order, then placing them from biggest to smallest on the core with fewst dataelements."""
function _distribute_subsystems_elements_greedy!(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}}, cores::Vector{CoreId})
   

    num_cores = length(cores)

    dist = Tuple{SubsystemIx, CoreId}[]
    # Initialize an array to hold the load (number of data elements) for each core
    core_loads = fill(0, num_cores)
    dataelements_in_each_subsystem = sort(get_subsystem_number_of_elements(subsystems), by=x -> x[2], rev=true)
    
    for (subsystem_index, num_elements) in dataelements_in_each_subsystem
        # Find the core with the least load
        min_load_core_index = argmin(core_loads)
        # Assign the subsystem to this core
        push!(dist,(subsystem_index,min_load_core_index))
        # Update the load for this core
        core_loads[min_load_core_index] += num_elements
    end
    return dist
end



"""function to get a tuple(subsystem index, number of data elements) that is sorted from highest to lowest number of data elements in the subsystem"""
function get_subsystem_number_of_elements(subsystems::Vector{Tuple{SubsystemIx, AbstractSubsystem}})
   
    dataelements_in_each_subsystem = Tuple{SubsystemIx, Int}[]
 
    for (ix, s) in subsystems
        
        data_elements_in_sub= get_dataelements(s)
        number_data_elements = length(data_elements_in_sub)
 
        # Create a tuple (data_element_index, number_data_elements)
        data_tuple = (ix, number_data_elements)
   
        # Add the tuple to the dataelements array
        push!(dataelements_in_each_subsystem, data_tuple)
    end
   
    return dataelements_in_each_subsystem
end


"""
Find which subsystems should have evp problems
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

""" Find which subsystems should have stoch problems"""
function get_subsystems_stoch(allsubsystems::Vector{AbstractSubsystem})
    subsystems = Tuple{SubsystemIx, AbstractSubsystem}[]
    for (i, subsystem) in enumerate(allsubsystems)
        if is_subsystem_stoch(subsystem)
            push!(subsystems, (i, subsystem))
        end
    end
    return subsystems
end