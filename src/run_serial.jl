"""
Possible re-design of JulES

Design goals
- small and clear code
- minimal communication
    - each local db stores the location of all problems
    - on-demand remote data collection with local cache
- time "everything"
- dynamic load balancer

"""
# JulES API and reference implementation

# TODO: setup docstrings for automatic documentation

function run_serial(output::AbstractJulESOutput, input::AbstractJulESInput)
    (t, N, delta, skipmed, skipmax) = init_jules(output, input)
    for stepnr in 1:N
        step_jules(output, t, delta, stepnr, skipmed)
        t += delta
        skipmed += Millisecond(delta)
        if skipmed > skipmax
            skipmed = Millisecond(0)
        end
    end
    cleanup_jules(output)
    return
end

function init_jules(output::AbstractJulESOutput, input::AbstractJulESInput)
    (t, N, delta, skipmed, skipmax) = get_simperiod(input)

    init_databases(input)

    # preallocate_output(output, input) TODO
    
    return (t, N, delta, skipmed, skipmax)
end

"""
Free local databases and clean-up temporary stuff in output-object
"""
function cleanup_jules(output::AbstractJulESOutput)
    @sync for core in get_cores(output)
        @spawnat core free_local_db()
    end
    cleanup_output(output)
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
            @spawnat core add_local_dummyobjects(core)
        end
    end

    println("Add local subsystems")
    @time begin
        c = first(cores)
        future = @spawnat c add_local_subsystems(c)
        wait(future)
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
        future = @spawnat c add_local_problem_distribution(c)
        wait(future)
    end

    println("Add local horizons")
    @time begin
        @sync for core in cores
            @spawnat core add_local_horizons(core)
        end
    end

    println("Add local problems")
    @time begin
        @sync for core in cores
            @spawnat core add_local_problems(core)
        end
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
# TODO: Only validate once
function add_local_dummyobjects(thiscore)
    db = get_local_db()

    # Horizons are needed to build modelobjects, but not used in scenario modelling
    dummyperiods = 10
    dummyperiodduration = Millisecond(Hour(24))
    dummyhorizon = SequentialHorizon(dummyperiods, dummyperiodduration)

    # Make dummy elements
    elements = copy(get_elements(db.input))
    commodities = get_settings(db)["horizons"]["commodities"]
    for commodity in commodities
        set_horizon!(elements, commodity, dummyhorizon)
        if commodity == "Power"
            set_horizon!(elements, "Battery", dummyhorizon)
        end
    end
    (dummyobjects, dummydeps) = getmodelobjects(elements, validate=true, deps=true)
    aggzonedict = Dict()
    for (k,v) in get_aggzone(get_settings(db))
        aggzonedict[Id(BALANCE_CONCEPT,"PowerBalance_" * k)] = [modelobjects[Id(BALANCE_CONCEPT,"PowerBalance_" * vv)] for vv in v]
    end
    aggzone!(dummyobjects, aggzonedict)
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
        (dummyobjects_ppp, dummydeps_ppp) = getmodelobjects(elements_ppp, validate=true, deps=true)
        db.dummyobjects_ppp = (dummyobjects_ppp, dummydeps_ppp)
    else
        db.dummyobjects_ppp = db.dummyobjects
    end
    return
end

"""
Make vector of subsystems
"""

function add_local_subsystems(thiscore)
    db = get_local_db()

    subsystems = get_subsystems(db)
    subsystems_evp = get_subsystems_evp(subsystems)
    subsystems_stoch = get_subsystems_stoch(subsystems)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core set_local_subsystems(subsystems, subsystems_evp, subsystems_stoch)
        end
    end

    return
end

# TODO: Implement this function for different methods
function get_subsystems(db)
    subsystems = []
    if get_onlysubsystemmodel(db.input)
        commodities = get_commoditites_from_dataelements(get_elements(db.input))
        endvaluemethod_sp = get_settings(db.input)["subsystems"]["endvaluemethod_sp"] # TODO: Parse to struct
        return push!(subsystems, ExogenSubsystem(commodities, endvaluemethod_sp))
    else
        settings = get_settings(db.input)
        method = settings["subsystems"]["function"]
        if method == "twostorageduration"
            modelobjects, dependencies = db.dummyobjects
            storagesystems = getstoragesystems(modelobjects)
            shorttermstoragesystems = getshorttermstoragesystems(storagesystems, Hour(settings["subsystems"]["shorttermstoragecutoff_hours"]))
            for storagesystem in shorttermstoragesystems
                subsystemdeps = Int[]
                for obj in storagesystem
                    objdeps = dependencies[getid(obj)]
                    for dep in objdeps
                        push!(subsystemdeps, dep)
                    end
                end
                commodities = get_commodities_from_storagesystem(storagesystem)
                priceareas = get_priceareas(storagesystem)
                subsystem = StochSubsystem(commodities, priceareas, unique(deps), Hour(settings["subsystems"]["shortstochduration_hours"]), "start_equal_stop")
                push!(subsystems, subsystem)
            end

            longtermstoragesystems = getlongtermstoragesystems(storagesystems, Hour(settings["subsystems"]["shorttermstoragecutoff_hours"]))
            for storagesystem in storagesystems
                subsystemdeps = Int[]
                for obj in storagesystem
                    objdeps = dependencies[getid(obj)]
                    for dep in objdeps
                        push!(subsystemdeps, dep)
                    end
                end
                commodities = get_commodities_from_storagesystem(storagesystem)
                priceareas = get_priceareas(storagesystem)
                subsystem = EVPSubsystem(commodities, priceareas, unique(deps), Day(settings["subsystems"]["longevduration_days"]), Day(settings["subsystems"]["longstochduration_days"]), "ppp")
                push!(subsystems, subsystem)
            end
        else
            error("getsubsystem() not implemented for $(method)")
        end
    end
    return subsystems
end
    
function get_commoditites_from_dataelements(elements::Vector{DataElement})
    commodities = CommodityName[]
    for element in elements
        if element.conceptname == BALANCE_CONCEPT
            commodity = element.value[COMMODITY_CONCEPT]
            if !(commodity in commodities)
                push!(commodities, commodity)
            end
        end
    end
    return commodities
end

function get_commodities_from_storagesystem(storagesystem::Vector)
    commodities = Set{CommodityName}()
    for obj in storagesystem
        if obj isa Flow
            for arrow in getarrows(obj)
                commodity = getcommodity(getbalance(arrow))
                if !(commodity in commodities)
                    push!(commodities, commodity)
                end
            end
        end
    end
    return commodities
end

function set_local_subsystems(subsystems, subsystems_evp, subsystems_stoch)
    db = get_local_db()
    
    db.subsystems = subsystems
    db.subsystems_evp = subsystems_evp
    db.subsystems_stoch = db.subsystems_stoch
    return
end

function get_priceareas(objects)
    priceareas = []
    for obj in objects
        if obj isa Flow
            for arrow in getarrows(obj)
                pricearea = getinstancename(getid(getbalance(arrow)))
                push!(priceareas, area)
            end
        end
    end
    return unique(areas)
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

    # Simulation scenario modelling - choose scenarios for the whole simulation
    numscen_sim = get_numscen_sim(db.input)
    @assert numscen_sim <= numscen_data
    if numscen_sim == numscen_data
        db.scenmod_sim = scenmod_data
    else
        db.scenmod_sim = get_scenmod(scenarios_data, settings["scenariogeneration"]["simulation"], numscen_sim, collect(values(first(get_dummyobjects_ppp(db)))))
    end

    # Prognosis scenario modelling - choose scenarios for the price prognosis models
    numscen_ppp = get_numscen_ppp(db.input)
    @assert numscen_ppp <= numscen_sim
    if numscen_ppp == numscen_sim
        db.scenmod_ppp = db.scenmod_sim
    else
        db.scenmod_ppp = get_scenmod(scenarios_data, settings["scenariogeneration"]["prognosis"], numscen_ppp, collect(values(first(get_dummyobjects_ppp(db)))))
    end

    # EVP scenario modelling - choose scenarios for the end values models
    numscen_evp = get_numscen_evp(db.input)
    @assert numscen_evp <= numscen_ppp
    if numscen_evp == numscen_ppp
        db.scenmod_evp = db.scenmod_ppp
    else
        db.scenmod_evp = get_scenmod(scenarios_data, settings["scenariogeneration"]["endvalue"], numscen_evp, collect(values(first(get_dummyobjects(db)))))
    end

    # Stochastic scenario modelling - choose scenarios for the price stochastic models
    numscen_stoch = get_numscen_stoch(db.input)
    @assert numscen_stoch <= numscen_evp
    if numscen_stoch == numscen_evp
        db.scenmod_stoch = db.scenmod_evp
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
function add_local_problem_distribution(thiscore)
    db = get_local_db()

    dist_ppp = get_dist_ppp(db.input)
    dist_evp = get_dist_evp(db.input, db.subsystems_evp)
    (dist_mp, dist_sp) = get_dist_stoch(db.input, db.subsystems_stoch)
    core_cp = get_core_cp(db.input)

    db.dist_ppp = dist_ppp
    db.dist_evp = dist_evp
    db.dist_sp = dist_sp
    db.dist_mp = dist_mp
    db.core_cp = core_cp

    dists = (dist_ppp, dist_evp, dist_sp, dist_mp, core_cp)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core set_local_dists(dists)
        end
    end

    return
end

function set_local_dists(dists)
    (dist_ppp, dist_evp, dist_sp, dist_mp, core_cp) = dists

    db = get_local_db()
    
    db.dist_ppp = dist_ppp
    db.dist_evp = dist_evp
    db.dist_sp = dist_sp
    db.dist_mp = dist_mp
    db.core_cp = core_cp
    
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
function add_local_horizons(thiscore)
    db = get_local_db()
    horizons = get_horizons(db.input)
    d = Dict{Tuple{ScenarioIx, TermName, CommodityName}, Horizon}()
    for (scenarioix, ownercore) in db.dist_ppp
        for ((term, commodity), horizon) in horizons
            if ownercore != thiscore
                horizon = getlightweightself(horizon)
                horizon = deepcopy(horizon)
                externalhorizon = ExternalHorizon(horizon)
                d[(scenarioix, term, commodity)] = externalhorizon
            else
                d[(scenarioix, term, commodity)] = horizon
            end
        end
    end
    db.horizons = d
    return
end

function add_local_problems(thiscore)
    db = get_local_db()

    for (scenix, core) in db.dist_ppp
        if core == thiscore
            create_ppp(db, scenix)
        end
    end

    for (scenix, subix, core) in db.dist_evp
        if core == thiscore
            create_evp(db, scenix, subix)
        end
    end

    for (subix, core) in db.dist_mp
        if core == thiscore
            create_mp(db, subsystem)
        end
    end

    for (scenix, subix, core) in db.dist_sp
        if core == thiscore
            create_sp(db, scenix, subix)
        end
    end

    if thiscore == db.core_cp
        create_cp(db)
    end
    return
end

# TODO: Use or remove delta
function step_jules(output::AbstractJulESOutput, t, delta, stepnr, skipmed)
    cores = get_cores(output)
    T = typeof(output) # So we can dispatch on output-type (to add extensibility)

    println(t)

    println("Price prognosis problems")
    @time begin
        c = first(cores)
        if stepnr == 1 
            f = @spawnat c update_scenmod_sim()
            wait(f)
        end
        f = @spawnat c update_scenmod_ppp(t, skipmed)
        wait(f)

        @sync for core in cores
            @spawnat core solve_ppp(t, delta, stepnr, skipmed)
        end
    end

    println("End value problems")
    @time begin
        # TODO: Add option to do scenariomodelling per individual or group of subsystem (e.g per area, commodity ...)
        f = @spawnat c update_scenmod_evp(t, skipmed)
        wait(f)

        @sync for core in cores
            @spawnat core solve_evp(t, delta, stepnr)
        end
    end

    println("Subsystem problems")
    @time begin
        # TODO: Add option to do scenariomodelling per individual or group of subsystem (e.g per area, commodity ...)
        f = @spawnat c update_scenmod_stoch(t, skipmed)
        wait(f)

        @sync for core in cores
            @spawnat core solve_mp(t, delta, stepnr)
        end
    end

    println("Clearing problem")
    @time begin
        @sync for core in cores
            @spawnat core solve_cp(t, delta, stepnr)
        end
    end

    # update_output(output, t, delta, stepnr) # TODO

    # do dynamic load balancing here
    return
end

# Principle for problem-solving: 
#    1. Collect information from other problems. 
#       If not cached locally, collect remote data and cache result locally. 
#       This (fetch-needed-data-only + cache) minimizes communication between cores, 
#       which is important for performance.
#    2. Update and solve problems
#    3. Possibly do syncronization
#
#    (each step should also store timing info for the load balancer and output-report)

# TODO: input parameters ok?

# TODO: consistent use of states and duals

function set_scenmodchanges_sim(changes)
    db = get_local_db()
    set_changes(db.scenmod_sim, changes)
end
function set_scenmodchanges_ppp(changes)
    db = get_local_db()
    set_changes(db.scenmod_ppp, changes)
end
function set_scenmodchanges_evp(changes)
    db = get_local_db()
    set_changes(db.scenmod_evp, changes)
end
function set_scenmodchanges_stoch(changes)
    db = get_local_db()
    set_changes(db.scenmod_stoch, changes)
end

function update_scenmod(scenmodmethod, scenmodmethodoptions, renumber, simtime, skipmed)
    db = get_local_db()

    if (length(get_scenarios(scenmodmethod)) != length(get_scenarios(scenmodmethodoptions))) && (skipmed.value == 0)
        choose_scenarios!(scenmodmethod, scenmodmethodoptions, simtime, db.input) # see JulES/scenariomodelling.jl
        renumber && renumber_scenmodmethod!(scenmodmethod)
    end
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
end
function update_scenmod_ppp(simtime, skipmed)
    db = get_local_db()

    update_scenmod(db.scenmod_ppp, db.scenmod_sim, true, simtime, skipmed)
    changes = get_changes(db.scenmod_ppp)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != db.core
            @spawnat core set_scenmodchanges_ppp(changes)
        end
    end
end
function update_scenmod_evp(simtime, skipmed)
    db = get_local_db()

    update_scenmod(db.scenmod_evp, db.scenmod_ppp, false, simtime, skipmed)
    changes = get_changes(db.scenmod_evp)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != db.core
            @spawnat core set_scenmodchanges_evp(changes)
        end
    end
end
function update_scenmod_stoch(simtime, skipmed)
    db = get_local_db()

    update_scenmod(db.scenmod_stoch, db.scenmod_evp, false, simtime, skipmed)
    changes = get_changes(db.scenmod_stoch)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != db.core
            @spawnat core set_scenmodchanges_stoch(changes)
        end
    end
end