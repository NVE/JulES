"""
This code does the orchestration of all the main simulation components:
- Initialization:
    - Local databases are created on all available cores and populated with:
        - User input, subsystems, scenarios, horizons and problem distribution
    - Specific problems are built on their allocated cores
    - One of the local databases initializes results. It will have the task of collecting results and timing from all 
        problems on all cores. Timing can be used to do dynamic load balancing (not implemented). Results and main
        problem (clearing or master problem) are on the same core to reduce data transfer overhead
- Each simulation step consist of
    - Scenario modelling
    - Update and solve problems on their allocated cores
    - Collect results and timing
    - Possibly do load balancing
- Minimal data transfer
    - Some data synced and distributed to all cores (horizons, scenarios, problem distribution, simulation state)
    - Most data is collected on demand and locally cached. This (fetch-needed-data-only + cache) minimizes 
        communication between cores, which is important for performance.
"""
# JulES API and reference implementation
# TODO: setup docstrings for automatic documentation

"""
Initialise the problems, run through all the simulation steps, store results and cleanup the databases
"""
function run_serial(input::AbstractJulESInput)
    (t, steps, steplength, skipmed, skipmax) = init_jules(input)
    totaltime = @elapsed for stepnr in 1:steps
        step_jules(t, steplength, stepnr, skipmed)
        t += steplength
        skipmed += Millisecond(steplength)
        if skipmed > skipmax
            skipmed = Millisecond(0)
        end
    end
    println(string("\nThe simulation took: ", round(totaltime/60; digits=2), " minutes"))
    println(string("Time usage per simulation step: ", round(totaltime/steps; digits=2), " seconds\n"))

    output = get_output_final(steplength, skipmax)
    cleanup_jules(input)
    return output
end

function init_jules(input::AbstractJulESInput)
    (t, steps, steplength, skipmed, skipmax) = get_simperiod(input)

    init_extensions(input)

    init_databases(input)
    
    return (t, steps, steplength, skipmed, skipmax)
end

function init_extensions(input::AbstractJulESInput)
    cores = get_cores(input)
    @sync for core in cores
        @spawnat core add_local_extensions()
    end
    return
end

function add_local_extensions()
    TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(ABSTRACT_INFLOW_MODEL, "TwoStateBucketIfm")] = includeTwoStateBucketIfm!
    TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(ABSTRACT_INFLOW_MODEL, "TwoStateNeuralODEIfm")] = includeTwoStateNeuralODEIfm!
    TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(TuLiPa.PARAM_CONCEPT, "ModeledInflowParam")] = includeModeledInflowParam!
    return
end

"""
Free local databases and clean-up temporary stuff.
Call gc on each core before returning.
"""
function cleanup_jules(input::AbstractJulESInput)
    cores = get_cores(input)
    @sync for core in cores
        @spawnat core free_local_db()
    end
    @sync for core in cores
        @spawnat core GC.gc()
    end
    return
end

"""
Create local db object on each core. Fill the db
with common information wich will be kept in sync.
The common inforation is 
- the user input
- the distribution of problems on cores
- horizons for each scenario, commodity and term (i.e. the time resolution)
In addition, problems are created on cores in accordance 
with the initial distribution. Which problems residing on which cores, may 
be changed at run time by dynamic load balancer. The precence of common (synced)
information on each core, will make it easier for the dynamic load balancer,
as this makes it easy to kill a problem on one core, and re-build it on another core.
"""
function init_databases(input::AbstractJulESInput)
    cores = get_cores(input)
    firstcore = first(cores)

    println("Add local dbs")
    @time begin
        @sync for core in cores
            @spawnat core create_local_db()
        end
    end

    println("Add local cores")
    @time begin
        @sync for core in cores
            @spawnat core add_local_core(core)
        end
    end

    println("Add local input")
    @time begin
        @sync for core in cores
            @spawnat core add_local_input(input)
        end
    end

    println("Add local dummyobjects")
    @time begin
        @sync for core in cores
            @spawnat core add_local_dummyobjects()
        end
    end

    println("Add local subsystems")
    @time begin
        wait(@spawnat firstcore add_local_subsystems())
    end

    println("Add local scenmod")
    @time begin
        @sync for core in cores
            @spawnat core add_local_scenariomodelling()
        end
    end
    
    # will calculate distribution on core c and then 
    # transfer this data to all other cores
    println("Add local problem distribution")
    @time begin
        wait(@spawnat firstcore add_local_problem_distribution())
    end

    println("Add local horizons")
    @time begin
        @sync for core in cores
            @spawnat core add_local_horizons()
        end
    end

    println("Add local problems")
    @time begin
        @sync for core in cores
            @spawnat core add_local_problems()
        end
        if !get_onlysubsystemmodel(input)
            @sync for core in cores
                @spawnat core add_local_cp()
            end
        end
    end

    println("Add local output")
    @time begin
        add_local_output()
    end
    return
end

"""
Store input in local db.
Input is treated read-only.
Input is serialized and copied to cores,
except the core that owns the input-object.
"""
function add_local_input(input::AbstractJulESInput)
    db = get_local_db()
    db.input = input
    return
end

"""
Store output on same core as cp
Needs to have all info and kept updated for load balancing
Updated after each simulation step
"""
function add_local_output()
    db = get_local_db()
    wait(@spawnat db.core_main init_local_output())
    return
end

"""
Add local core
"""
function add_local_core(core::CoreId)
    db = get_local_db()
    db.core = core
    return
end

"""
Build dummyobjects on each core
For use in scenario modelling, validate elements and collect storages
"""
function add_local_dummyobjects()
    db = get_local_db()

    # Horizons are needed to build modelobjects, but not used in scenario modelling
    dummyperiods = 10
    dummyperiodduration = Millisecond(Hour(24))
    dummyhorizon = TuLiPa.SequentialHorizon(dummyperiods, dummyperiodduration)

    # Make dummy elements
    elements = copy(get_elements(db.input))
    commodities = get_settings(db)["horizons"]["commodities"]
    for commodity in commodities
        set_horizon!(elements, commodity, dummyhorizon)
        if commodity == "Power"
            set_horizon!(elements, "Battery", dummyhorizon)
        end
    end
    (dummyobjects, dummydeps) = TuLiPa.getmodelobjects(elements, validate=true, deps=true)
    aggzonedict = Dict()
    for (k,v) in get_aggzone(get_settings(db))
        aggzonedict[TuLiPa.Id(TuLiPa.BALANCE_CONCEPT,"PowerBalance_" * k)] = [dummyobjects[TuLiPa.Id(TuLiPa.BALANCE_CONCEPT,"PowerBalance_" * vv)] for vv in v]
    end
    TuLiPa.aggzone!(dummyobjects, aggzonedict)
    db.dummyobjects = (dummyobjects, dummydeps)

    # Make dummy prog elements
    if haskey(db.input.dataset, "elements_ppp")
        elements_ppp = copy(get_elements_ppp(db.input))
        for commodity in commodities
            set_horizon!(elements_ppp, commodity, dummyhorizon)
            if commodity == "Power"
                set_horizon!(elements_ppp, "Battery", dummyhorizon)
            end
        end
        (dummyobjects_ppp, dummydeps_ppp) = TuLiPa.getmodelobjects(elements_ppp, validate=true, deps=true)
        db.dummyobjects_ppp = (dummyobjects_ppp, dummydeps_ppp)
    else
        db.dummyobjects_ppp = db.dummyobjects
    end
    return
end # TODO: Only validate once

"""
Make vector of subsystems
"""
function add_local_subsystems()
    db = get_local_db()

    subsystems = create_subsystems(db)
    subsystems_evp = get_subsystems_evp(subsystems)
    subsystems_stoch = get_subsystems_stoch(subsystems)

    cores = get_cores(db.input)
    @sync for core in cores
        @spawnat core set_local_subsystems(subsystems, subsystems_evp, subsystems_stoch)
    end

    return
end

# TODO: Implement this function for different methods
function create_subsystems(db)
    elements = get_elements(db.input)
    subsystems = AbstractSubsystem[]
    modelobjects, dependencies = db.dummyobjects
    deep_dependencies = TuLiPa.get_deep_dependencies(elements, dependencies)
    # filtered_dependencies = get_filtered_dependencies(elements, dependencies)
    # deep_dependencies = get_deep_dependencies(elements, filtered_dependencies; concepts=[PARAM_CONCEPT, METADATA_CONCEPT])
    if get_onlysubsystemmodel(db.input)
        commodities = get_commodities_from_dataelements(get_elements(db.input))
        endvaluemethod_sp = get_settings(db.input)["problems"]["stochastic"]["endcondition"] # TODO: Parse to struct
        return push!(subsystems, ExogenSubsystem(commodities, endvaluemethod_sp))
    else
        settings = get_settings(db.input)
        method = settings["subsystems"]["function"]
        if method == "twostorageduration"
            storagesystems = TuLiPa.getstoragesystems(modelobjects)
            shorttermstoragesystems = TuLiPa.getshorttermstoragesystems(storagesystems, Hour(settings["subsystems"]["shorttermstoragecutoff_hours"]))
            for storagesystem in shorttermstoragesystems
                commodities = get_commodities_from_storagesystem(storagesystem)
                main = Set()
                all = Set()
                for obj in storagesystem
                    i, element = get_element_from_obj(elements, obj)
                    # println(getelkey(elements[i]))
                    push!(main, i)
                    push!(all, i)
                end

                for (_i, _element) in enumerate(elements)
                    _deps = deep_dependencies[_element]
                    _add = false
                    for _dep in _deps
                        if _dep in main
                            _add = true
                        end
                    end
                    if _add
                        for _dep in _deps
                            if !(_dep in all)
                                elkey = TuLiPa.getelkey(elements[_dep])
                                if elkey.conceptname != TuLiPa.BALANCE_CONCEPT # getstoragesystems have already picked the balances we want to include, ignores power balances
                                    # println(elkey)
                                    push!(all, _dep)
                                end
                            end
                        end
                    end
                end
                # println(length(all))

                shortstochduration = parse_duration(settings["subsystems"], "shortstochduration")
                horizonterm_stoch = get_term_ppp(get_horizons(db.input), commodities, shortstochduration)

                priceareas = get_priceareas(storagesystem)
                skipmed_impact = false
                subsystem = StochSubsystem(commodities, priceareas, collect(all), horizonterm_stoch, shortstochduration, "startequalstop", skipmed_impact)
                push!(subsystems, subsystem)
            end
            num_shortterm = length(subsystems)
            println("Number of shortterm storagesystems $num_shortterm")

            longtermstoragesystems = TuLiPa.getlongtermstoragesystems(storagesystems, Hour(settings["subsystems"]["shorttermstoragecutoff_hours"]))
            for storagesystem in longtermstoragesystems
                commodities = get_commodities_from_storagesystem(storagesystem)
                if length(commodities) == 1
                    continue # TODO: error and fix dataset linvasselv and vakkerjordvatn have two subsystems, one not connected to power market, send liste til Carl 
                end  

                # all = Set()
                # for obj in storagesystem
                #     i, element = get_element_from_obj(elements, obj)
                #     for dep in deep_dependencies[element]
                #         if !(dep in all)
                #             println(getelkey(elements[dep]))
                #             push!(all, dep)
                #         end
                #     end
                # end 

                main = Set()
                all = Set()
                for obj in storagesystem
                    i, element = get_element_from_obj(elements, obj)
                    # println(getelkey(elements[i]))
                    push!(main, i)
                    push!(all, i)
                end

                for (_i, _element) in enumerate(elements)
                    _deps = deep_dependencies[_element]
                    _add = false
                    for _dep in _deps
                        if _dep in main
                            _add = true
                        end
                    end
                    if _add
                        for _dep in _deps
                            if !(_dep in all)
                                elkey = TuLiPa.getelkey(elements[_dep])
                                if elkey.conceptname != TuLiPa.BALANCE_CONCEPT # getstoragesystems have already picked the balances we want to include, ignores power balances
                                    # println(elkey)
                                    push!(all, _dep)
                                end
                            end
                        end
                    end
                end
                # println(length(all))

                # completed = Set()
                # remaining = Set()
                # for obj in storagesystem
                #     i, element = get_element_from_obj(elements, obj)
                #     push!(remaining, i)
                # end
                # while length(remaining) > 0
                #     i = pop!(remaining)
                #     push!(completed, i)
                #     # # Deps over
                #     elkey = getelkey(elements[i])
                #     println(elkey)
                #     for dep in dependencies[elkey]
                #         if !(dep in completed) && (dep <= length(elements))
                #             if !(elkey.conceptname in ["Commodity", "Arrow"])
                #                 push!(remaining, dep)
                #             end
                #         end
                #     end
                #     # Deps under
                #     for (_i, _element) in enumerate(elements)
                #         _elkey = getelkey(_element)
                #         if (_elkey != elkey)
                #             _deps = dependencies[_elkey]
                #             for _dep in _deps
                #                 if _dep == i
                #                     if (_dep in completed) && !(_elkey.conceptname in ["TimeVector", "TimeValues"])
                #                         push!(remaining, _i)
                #                     end
                #                 end
                #             end
                #         end
                #     end
                # end
                # println(length(completed))

                longstochduration = parse_duration(settings["subsystems"], "longstochduration")
                horizonterm_stoch = get_term_ppp(get_horizons(db.input), commodities, longstochduration)

                priceareas = get_priceareas(storagesystem)
                skipmed_impact = true  
                if has_longevduration(settings)
                    longevduration = parse_duration(settings["subsystems"], "longevduration")
                    horizonterm_evp = get_term_ppp(get_horizons(db.input), commodities, longevduration)

                    subsystem = EVPSubsystem(commodities, priceareas, collect(all), horizonterm_evp, longevduration, horizonterm_stoch, longstochduration, "ppp", skipmed_impact)
                else
                    subsystem = StochSubsystem(commodities, priceareas, collect(all), horizonterm_stoch, longstochduration, "ppp", skipmed_impact)
                end
                push!(subsystems, subsystem)
            end
            println("Number of longterm storagesystems $(length(subsystems)-num_shortterm)")
            num_ignored = length(shorttermstoragesystems) + length(longtermstoragesystems) - length(subsystems)
            println("Number of ignored storagesystems not connected to power $num_ignored")
        else
            error("getsubsystem() not implemented for $(method)")
        end
    end
    return subsystems
end

function has_longevduration(settings)
    for key in keys(settings["subsystems"])
        if startswith(key, "longevduration")
            return true
        end
    end
    return false  
end

# Which time resolution (short, med, long) should we use horizons and prices from
# TODO: Should we use different terms for master and subproblems?
function get_term_ppp(horizons, commodities, duration)
    dummycommodity = commodities[1]
    horizon_short = horizons[(ShortTermName, dummycommodity)]
    if duration <= TuLiPa.getduration(horizon_short) # TODO: also account for slack in case of reuse of watervalues
        return ShortTermName
    end
    horizon_med = horizons[(MedTermName, dummycommodity)]
    if duration <= TuLiPa.getduration(horizon_med) # TODO: also account for slack in case of reuse of watervalues
        return MedTermName
    end
    horizon_long = horizons[(LongTermName, dummycommodity)]
    @assert duration <= TuLiPa.getduration(horizon_long) # TODO: also account for slack in case of reuse of watervalues
    return LongTermName   
end

function get_filtered_dependencies(elements, dependencies)
    filtered_dependencies = Dict{TuLiPa.ElementKey,Vector{Int}}()
    for element in elements # remove dependencies of elemements not in elements list
        filtered = [x for x in dependencies[TuLiPa.getelkey(element)] if x < length(elements)]
        filtered_dependencies[TuLiPa.getelkey(element)] = filtered
    end
    return filtered_dependencies
end

function get_elkey_from_element(dataelements::Vector{TuLiPa.DataElement}, element::TuLiPa.DataElement)
    for _element in dataelements
        if _element == element
            return TuLiPa.getelkey(element)
        end
    end
    error("element not in dataelements")
end

function get_element_from_obj(dataelements::Vector{TuLiPa.DataElement}, obj::Any)
    objid = TuLiPa.getid(obj)
    conceptname = objid.conceptname
    instancename = objid.instancename
    for (i, dataelement) in enumerate(dataelements)
        if (dataelement.conceptname == conceptname) && (dataelement.instancename == instancename)
            return (i, dataelement)
        end
    end
end
    
function get_commodities_from_dataelements(elements::Vector{TuLiPa.DataElement})
    commodities = CommodityName[]
    for element in elements
        if element.conceptname == TuLiPa.BALANCE_CONCEPT
            commodity = element.value[TuLiPa.COMMODITY_CONCEPT]
            if !(commodity in commodities)
                push!(commodities, commodity)
            end
        end
    end
    return commodities
end

function get_obj_from_id(objects::Vector, id::TuLiPa.Id)
    for obj in objects
        if TuLiPa.getid(obj) == id
            return obj
        end
    end
end

function get_commodities_from_storagesystem(storagesystem::Vector)
    commodities = Set{CommodityName}()
    for obj in storagesystem
        if obj isa TuLiPa.Flow
            for arrow in TuLiPa.getarrows(obj)
                commodity = TuLiPa.getcommodity(TuLiPa.getbalance(arrow))
                if !(commodity in commodities)
                    push!(commodities, TuLiPa.getinstancename(TuLiPa.getid(commodity)))
                end
            end
        end
    end
    return collect(commodities)
end

function set_local_subsystems(subsystems, subsystems_evp, subsystems_stoch)
    db = get_local_db()
    
    db.subsystems = subsystems
    db.subsystems_evp = subsystems_evp
    db.subsystems_stoch = subsystems_stoch
    return
end

function get_priceareas(objects)
    priceareas = []
    for obj in objects
        if obj isa TuLiPa.Flow
            for arrow in TuLiPa.getarrows(obj)
                balance = TuLiPa.getbalance(arrow)
                if (TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(balance))) == "Power") && !TuLiPa.isexogen(balance)
                    pricearea = TuLiPa.getinstancename(TuLiPa.getid(balance))
                    push!(priceareas, pricearea)
                end
            end
        end
    end
    return unique(priceareas)
end


"""
Initial scenario modelling for simulation, prognosis and stochastic
"""
function add_local_scenariomodelling()
    db = get_local_db()
    scenmod_data = get_scenmod_data(db)
    numscen_data = get_numscen_data(db)
    scenarios_data = get_scenarios(scenmod_data)
    settings = get_settings(db)

    # Simulation scenario modelling - choose scenarios for price prognosis and endvalue problems for the whole simulation
    # TODO: Add possibility to update simulation scenarios (for long simulations)
    numscen_sim = get_numscen_sim(db.input)
    @assert numscen_sim <= numscen_data
    if numscen_sim == numscen_data
        db.scenmod_sim = scenmod_data
    else
        db.scenmod_sim = get_scenmod(scenarios_data, settings["scenariogeneration"]["simulation"], numscen_sim, collect(values(first(get_dummyobjects_ppp(db)))))
    end

    # Stochastic scenario modelling - choose scenarios for the price stochastic models
    numscen_stoch = get_numscen_stoch(db.input)
    @assert numscen_stoch <= numscen_sim
    if numscen_stoch == numscen_sim
        db.scenmod_stoch = db.scenmod_sim
    else
        db.scenmod_stoch = get_scenmod(scenarios_data, settings["scenariogeneration"]["stochastic"], numscen_stoch, collect(values(first(get_dummyobjects(db)))))
    end

    return
end

"""
Initial allocation of problems to cores.
Done in order to minimize difference in expected work 
by different cores, based on static information (i.e. differences in problem size).
As static information will leave some room for improvement, the initial allocation 
may be improved by dynamic load balancer during simulation, in order to optimize 
further (i.e. move problem-core-location) based on collected timing data.
The initial distribution on cores is calculated on one core, and then transfered
to all other cores. 
"""
function add_local_problem_distribution()
    db = get_local_db()

    if !get_onlysubsystemmodel(db.input)
        db.dist_ppp = get_dist_ppp(db.input)
        println(db.dist_ppp)
        db.dist_evp = get_dist_evp(db.input, db.subsystems_evp)
        println(db.dist_evp)
    end
    db.dist_ifm = get_dist_ifm(db.input)
    println(db.dist_ifm)
    (dist_mp, dist_sp) = get_dist_stoch(db.input, db.subsystems_stoch)
    println(dist_mp)
    println(dist_sp)
    db.dist_sp = dist_sp
    db.dist_mp = dist_mp

    db.core_main = get_core_main(db.input)
    println(db.core_main)

    dists = (db.dist_ifm, db.dist_ppp, db.dist_evp, db.dist_sp, db.dist_mp, db.core_main)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != db.core
            @spawnat core set_local_dists(dists)
        end
    end

    return
end

function set_local_dists(dists)
    db = get_local_db()

    (dist_ifm, dist_ppp, dist_evp, dist_sp, dist_mp, core_main) = dists
    db.dist_ifm = dist_ifm
    db.dist_ppp = dist_ppp
    db.dist_evp = dist_evp
    db.dist_sp = dist_sp
    db.dist_mp = dist_mp
    db.core_main = core_main
    
    return
end

"""
Add set of horizons (for all scenarios, terms and commodities) on each core, with a promise to keep them in sync.
The "master" horizons come from the db.ppp (price prognosis problems). When we update and solve a ppp,
the horizons in this problem may (depending on type) be updated (e.g. we may change durations for some periods, 
or change which hours are mapped to load blocks). We want to store a synced set of all horizons on all cores 
(to enable dynamic load balancing), so if a master horizon changes, we need to transfer this information to all 
non-master copies of this horizon residing on other cores than the core holding the master. Since we for 
non-master horizons do update-by-transfer, we want to turn off update-by-solve behaviour for these horizons. 
Hence, we wrap them in ExternalHorizon, which specializes the update! method to do nothing. 
Which cores own which scenarios are defined in db.dist_ppp at any given time. 
"""
function add_local_horizons()
    db = get_local_db()
    horizons = get_horizons(db.input)

    d = Dict{Tuple{ScenarioIx, TermName, CommodityName}, TuLiPa.Horizon}()

    if !get_onlysubsystemmodel(db.input)
        for (scenarioix, ownercore) in db.dist_ppp
            for ((term, commodity), horizon) in horizons
                if (ownercore != db.core) && !get_onlysubsystemmodel(db.input)
                    horizon = TuLiPa.getlightweightself(horizon)
                    horizon = deepcopy(horizon)
                    externalhorizon = TuLiPa.ExternalHorizon(horizon)
                    d[(scenarioix, term, commodity)] = externalhorizon
                elseif ownercore == db.core
                    horizon = deepcopy(horizon) # TODO: Only deepcopy parts of horizon
                    d[(scenarioix, term, commodity)] = horizon
                end
            end
        end
    else
        commoditites = get_settings(db.input)["horizons"]["commodities"]
        for commodity in commoditites
            term = MasterTermName
            for (subix, ownercore) in db.dist_mp
                if ownercore == db.core
                    d[(1, term, commodity)] = horizons[term, commodity]
                end
            end
            term = SubTermName
            for (scenarioix, subix, ownercore) in db.dist_sp
                if ownercore == db.core
                    d[(scenarioix, term, commodity)] = deepcopy(horizons[term, commodity])
                end
            end
        end
    end
    db.horizons = d
    return
end

function add_local_cp()
    db = get_local_db()
    if db.core == db.core_main
        create_cp()
    end
    return
end

function add_local_problems()
    db = get_local_db()

    create_ifm()

    for (scenix, core) in db.dist_ppp
        if core == db.core
            create_ppp(scenix)
        end
    end

    for (scenix, subix, core) in db.dist_evp
        if core == db.core
            create_evp(scenix, subix)
        end
    end

    for (subix, core) in db.dist_mp
        if core == db.core
            create_mp(subix)
        end
    end

    for (scenix, subix, core) in db.dist_sp
        if core == db.core
            create_sp(scenix, subix)
        end
    end
    return
end

# TODO: Use or remove delta
function step_jules(t, steplength, stepnr, skipmed)
    db = get_local_db()
    cores = get_cores(db)
    firstcore = first(cores)

    println(t)
    println("Garbage collection")
    @time begin
        if mod(stepnr, 20) == 0
            @sync for core in cores
                @spawnat core GC.gc()
            end
        end
    end
    
    println("Startstates")
    @time begin
        @sync for core in cores
            @spawnat core update_startstates(stepnr, t)
        end
    end
    
    println("Scenario modelling")
    @time begin
        if stepnr == 1 
            wait(@spawnat firstcore update_scenmod_sim())
        end
        # TODO: Add option to do scenariomodelling per individual or group of subsystem (e.g per area, commodity ...)
        # TODO: Do scenario modelling based on ifm
        wait(@spawnat firstcore update_scenmod_stoch(t, skipmed))
    end

    println("Solve inflow models")
    @time begin
        @sync for core in cores
            @spawnat core solve_ifm(t, stepnr)
        end

        @sync for core in cores
            @spawnat core synchronize_ifm_output()
        end

        @sync for core in cores
            @spawnat core update_ifm_derived()
        end
    end

    println("Solve price prognosis problems")
    @time begin
        @sync for core in cores
            @spawnat core solve_ppp(t, steplength, stepnr, skipmed)
        end

        @sync for core in cores
            @spawnat core synchronize_horizons(skipmed)
        end
    end
	
    println("End value problems")
    @time begin
        @sync for core in cores
            @spawnat core solve_evp(t, stepnr, skipmed)
        end
    end

    println("Subsystem problems")
    @time begin
        @sync for core in cores
            @spawnat core solve_stoch(t, stepnr, skipmed)
        end
    end

    println("Clearing problem")
    @time begin
        if !get_onlysubsystemmodel(db.input)
            @sync for core in cores
                @spawnat core solve_cp(t, stepnr, skipmed)
            end
        end
    end

    println("Update output")
    @time begin
        wait(@spawnat db.core_main update_output(t, stepnr))
    end
	
    # do dynamic load balancing here
    return
end

# TODO: input parameters ok?

# TODO: consistent use of states and duals

function set_scenmodchanges_sim(changes)
    db = get_local_db()
    set_changes(db.scenmod_sim, changes)
    return
end
function set_scenmodchanges_stoch(changes)
    db = get_local_db()
    set_changes(db.scenmod_stoch, changes)
    return
end

function update_scenmod(scenmodmethod, scenmodmethodoptions, renumber, simtime, skipmed)
    db = get_local_db()

    if (length(get_scenarios(scenmodmethod)) != length(get_scenarios(scenmodmethodoptions))) && (skipmed.value == 0)
        choose_scenarios!(scenmodmethod, scenmodmethodoptions, simtime, db.input) # see JulES/scenariomodelling.jl
        renumber && renumber_scenmodmethod!(scenmodmethod)
    end
    return
end

function update_scenmod_sim()
    db = get_local_db()

    update_scenmod(db.scenmod_sim, get_scenmod_data(db), true, get_simstarttime(db.input), Millisecond(0))
    changes = get_changes(db.scenmod_sim)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != db.core
            @spawnat core set_scenmodchanges_sim(changes)
        end
    end
    return
end
function update_scenmod_stoch(simtime, skipmed)
    db = get_local_db()

    update_scenmod(db.scenmod_stoch, db.scenmod_sim, false, simtime, skipmed)
    changes = get_changes(db.scenmod_stoch)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != db.core
            @spawnat core set_scenmodchanges_stoch(changes)
        end
    end
    return
end