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
    (t, N, delta, skipmed, skipmax) = get_simulation_period(input)

    init_databases(input)

    preallocate_output(output, input)
    
    return (t, N, delta)
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

    @sync for core in cores
        @spawnat core create_local_db()
    end

    @sync for core in cores
        @spawnat core add_local_input(input)
    end

    @sync for core in cores
        @spawnat core add_local_dummyobjects(input)
    end

    @sync for core in cores
        @spawnat core add_local_scenmodmethods(input)
    end

    c = first(cores)
    # will calculate distribution on core c and then 
    # transfer this data to all other cores
    future = @spawnat c add_local_problem_distribution(c)
    wait(future)

    @sync for core in cores
        @spawnat core add_local_horizons(core)
    end

    @sync for core in cores
        @spawnat core add_local_problems(core)
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
Build dummyobjects on each core
For use in scenario modelling, validate elements and collect storages
"""
# TODO: Only validate once
function add_local_dummyobjects(thiscore)
    db = get_local_db()

    # Horizons are needed to build modelobjects, but not used in scenario modelling
    dummyperiods = 10
    dummyperiodduration = Millisecond(Hour(24))
    power_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)
    hydro_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)

    # Make dummy elements
    db.dummyobjects = make_obj(db.input.dataset["elements"], hydro_horizon, power_horizon, validate=true)

    # Make dummy prog elements
    if haskey(db.input.dataset, "progelements")
        db.dummyprogobjects = make_obj(db.input.dataset["progelements"], hydro_horizon, power_horizon, validate=true)
    else
        db.dummyprogobjects = db.dummyobjects
    end
    return
end

"""
Initial scenario modelling for simulation, prognosis and stochastic

"""
function add_local_scenmodemethods()
    db = get_local_db()
    datascenmodmethod = db.input.datascenmodmethod
    datanumscen = length(datascenmodmethod)
    settings = db.input.settings

    # Simulation scenario modelling - choose scenarios for the whole simulation
    simnumscen = getnumscen_sim(db.input)
    @assert simnumscen <= datanumscen
    if simnumscen == datanumscen
        db.simscenmodmethod = db.datascenmodmethod
    else
        db.simscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["simulation"], simnumscen, values(db.input.dummyprogobjects))
    end

    # Prognosis scenario modelling - choose scenarios for the price prognosis models
    prognumscen = getnumscen_ppp(db.input)
    @assert prognumscen <= simnumscen
    if prognumscen == simnumscen
        db.progscenmodmethod = db.simscenmodmethod
    else
        db.progscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["prognosis"], prognumscen, values(db.input.dummyprogobjects))
    end

    # EVP scenario modelling - choose scenarios for the end values models
    evnumscen = getnumscen_evp(db.input)
    @assert evnumscen <= prognumscen
    if evnumscen == prognumscen
        db.evscenmodmethod = db.progscenmodmethod
    else
        db.evscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["endvalue"], evnumscen, values(db.input.dummyobjects))
    end

    # Stochastic scenario modelling - choose scenarios for the price stochastic models
    stochnumscen = getnumscen_sp(db.input)
    @assert stochnumscen <= evnumscen
    if stochnumscen == evnumscen
        db.stochscenmodmethod = db.progscenmodmethod
    else
        db.stochscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["stochastic"], stochnumscen, values(db.input.dummyobjects))
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

    ppp_dist = get_ppp_dist(db.input)
    evp_dist = get_evp_dist(db.input)
    (mp_dist, sp_dist) = get_mp_sp_dist(db.input)
    cp_core = get_cp_core(db.input)

    db.ppp_dist = ppp_dist
    db.evp_dist = evp_dist
    db.sp_dist = sp_dist
    db.mp_dist = mp_dist
    db.cp_core = cp_core

    dists = (ppp_dist, evp_dist, sp_dist, mp_dist, cp_core)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core set_local_dists(dists)
        end
    end

    return
end

function set_local_dists(dists)
    (ppp_dist, evp_dist, sp_dist, mp_dist, cp_core) = dists

    db = get_local_db()
    
    db.ppp_dist = ppp_dist
    db.evp_dist = evp_dist
    db.sp_dist = sp_dist
    db.mp_dist = mp_dist
    db.cp_core = cp_core
    
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
Which cores own which scenarios are defined in db.ppp_dist at any given time. 
"""
function add_local_horizons(thiscore)
    db = get_local_db()
    horizons = get_horizons(db.input)
    d = Dict{Tuple{ScenarioIx, TermName, CommodityName}, Horizon}()
    for (scenarioix, ownercore) in db.ppp_dist
        for ((term, commodity), horizon) in horizons
            horizon = getlightweightself(horizon)
            horizon = deepcopy(horizon)
            if ownercore != thiscore
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

    d = Dict{ScenarioIx, PricePrognosisProblem}()
    for (scenarioix, core) in db.ppp_dist
        if core == thiscore
            d[scenarioix] = create_ppp(db, scenario)
        end
    end
    db.ppp = d

    d = Dict{Tuple{ScenarioIx, SubsystemIx}, EndValueProblem}()
    for (scenario, subsystem, core) in db.evp_dist
        if core == thiscore
            d[(scenarioix, subsystem)] = create_evp(db, scenario, subsystem)
        end
    end
    db.evp = d

    d = Dict{SubsystemIx, MasterProblem}()
    for (scenarioix, core) in db.mp_dist
        if core == thiscore
            d[subsystem] = create_mp(db, subsystem)
        end
    end
    db.mp = d

    d = Dict{Tuple{ScenarioIx, SubsystemIx}, ScenarioProblem}()
    for (scenario, subsystem, core) in db.sp_dist
        if core == thiscore
            d[(scenarioix, subsystem)] = create_sp(db, scenarioix, subsystem)
        end
    end
    db.sp = d

    if thiscore == db.cp_core
        db.cp = create_cp(db)
    end
    return
end

function step_jules(output::AbstractJulESOutput, t, delta, stepnr, skipmed)
    cores = get_cores(output)

    # do input models here

    # Do global scenario modelling
    c = first(cores)
    if stepnr == 1 
        f = @spawnat c update_simscenariomodelling!(c, skipmed)
        wait(f)
    end
    f = @spawnat c update_progscenariomodelling!(c, t, skipmed)
    wait(f)

    T = typeof(output) # So we can dispatch on output-type (to add extensibility)

    @sync for core in cores
        @spawnat core solve_ppp(T, t, delta, stepnr, skipmed, core)
    end

    # TODO: Add option to do scenariomodelling per individual or group of subsystem (e.g per area, commodity ...)
    f = @spawnat c update_evpscenariomodelling!(c, t, skipmed)
    wait(f)

    @sync for core in cores
        @spawnat core solve_evp(T, t, delta, stepnr, core)
    end

    # TODO: Add option to do scenariomodelling per individual or group of subsystem (e.g per area, commodity ...)
    f = @spawnat c update_stochscenariomodelling!(c, t, skipmed)
    wait(f)

    @sync for core in cores
        @spawnat core solve_mp(T, t, delta, stepnr, core)
    end

    @sync for core in cores
        @spawnat core solve_cp(T, t, delta, stepnr, core)
    end

    update_output(output, t, delta, stepnr)

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
function solve_evp(T, t, delta, stepnr, thiscore)
    update_startstates_evp(stepnr)
    update_endstates_evp()
    update_prices_evp()
    solve_local_evp(t)
    return
end

function solve_mp(T, t, delta, stepnr, thiscore)
    update_startstates_mp(stepnr)
    update_endstates_sp(stepnr)
    scale_inflow_sp(stepnr)
    update_prices_mp()
    update_prices_sp()
    update_statedependent_mp()
    update_mp(t)
    update_sp(t)
    solve_benders(stepnr)
    final_solve_mp()
    return
end

function solve_cp(T, t, delta, stepnr, thiscore)
    db = get_local_db()
    if thiscore == db.cp_core        
        db.time_cp_startstates = @elapsed update_startstates_cp(stepnr)
        db.time_cp_endstates   = @elapsed update_endstates_cp()
        db.time_cp_cuts        = @elapsed update_cuts_cp()
        db.time_cp_update      = @elapsed update!(db.cp, t)
        db.time_cp_solve       = @elapsed solve!(db.cp)
    end
    return
end

function update_startstates_evp(stepnr)

end

function setchanges_simscenmodmethod(changes)
    db = get_local_db()
    setchanges(db.simscenmodmethod, changes)
end
function setchanges_progscenmodmethod(changes)
    db = get_local_db()
    setchanges(db.progscenmodmethod, changes)
end
function setchanges_evscenmodmethod(changes)
    db = get_local_db()
    setchanges(db.evscenmodmethod, changes)
end
function setchanges_stochscenmodmethod(changes)
    db = get_local_db()
    setchanges(db.stochscenmodmethod, changes)
end

function update_scenariomodelling!(thiscore, scenmodmethod, scenmodmethodoptions, renumber, simtime, skipmed)
    db = get_local_db()

    if (length(scenmodmethod) != length(scenmodmethodoptions)) && (skipmed.value == 0)
        choose_scenarios!(scenmodmethod, scenmodmethodoptions, simtime, db.input) # see JulES/scenariomodelling.jl
        renumber && renumber_scenmodmethod!(scenmodmethod)
    end
end

function update_simscenariomodelling!(thiscore)
    db = get_local_db()

    update_scenariomodelling!(thiscore, db.simscenmodmethod, db.input.datascenmodmethod, true, db.input.simstarttime, Millisecond(0))
    changes = getchanges(db.simscenmodmethod)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core setchanges_simscenmodmethod(changes)
        end
    end
end
function update_progscenariomodelling!(thiscore, simtime, skipmed)
    db = get_local_db()

    update_scenariomodelling!(thiscore, db.progscenmodmethod, db.simscenmodmethod, true, simtime, skipmed)
    changes = getchanges(db.progscenmodmethod)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core setchanges_progscenmodmethod(changes)
        end
    end
end
function update_evscenariomodelling!(thiscore, simtime, skipmed)
    db = get_local_db()

    update_scenariomodelling!(thiscore, db.evscenmodmethod, db.progscenmodmethod, false, simtime, skipmed)
    changes = getchanges(db.evscenmodmethod)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core setchanges_evscenmodmethod(changes)
        end
    end
end
function update_stochscenariomodelling!(thiscore, simtime, skipmed)
    db = get_local_db()

    update_scenariomodelling!(thiscore, db.stochscenmodmethod, db.evscenmodmethod, false, simtime, skipmed)
    changes = getchanges(db.stochscenmodmethod)

    cores = get_cores(db.input)
    @sync for core in cores
        if core != thiscore
            @spawnat core setchanges_stochscenmodmethod(changes)
        end
    end
end