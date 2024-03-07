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
    (t, N, delta) = init_jules(output, input)
    for stepnr in 1:N
        step_jules(output, t, delta, stepnr)
        t += delta
    end
    cleanup_jules(output)
    return
end

function init_jules(output::AbstractJulESOutput, input::AbstractJulESInput)
    (t, N, delta) = get_simulation_period(input)

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
        @spawnat core add_local_problem_distribution()
    end

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
Initial allocation of problems to cores.
Done in order to minimize difference in expected work 
by different cores, based on static information (i.e. differences in problem size).
As static information will leave some room for improvement, the initial allocation 
may be improved by dynamic load balancer during simulation, in order to optimize further 
(i.e. move problem-core-location) based on collected timing data.
"""
function add_local_problem_distribution()
    db = get_local_db()
    db.pp_dist = get_pp_dist(db.input)   # must give same result on each core (deterministic)
    db.ev_dist = get_ev_dist(db.input)   # must give same result on each core (deterministic)
    db.sp_dist = get_sp_dist(db.input)   # must give same result on each core (deterministic)
    db.mp_dist = get_mp_dist(db.input)   # must give same result on each core (deterministic)
    db.cp_core = get_cp_core(db.input)   # must give same result on each core (deterministic)
    return
end

"""
Add set of horizons (for all scenarios, terms and commodities) on each core, with a promise to keep them in sync.
The "master" horizons come from the db.pp (price prognosis) problems. When we update and solve a pp problem,
the horizons in this problem may (depending on type) be updated (e.g. we may change durations for some periods, 
or change which hours are mapped to load blocks). We want to store a synced set of all horizons on all cores 
(to-enable-dynamic-load-balancing), so if a master horizon changes, we need to transfer this information to all 
non-master copies of this horizon residing on other cores than the core holding the master. Since we update 
non-master-horizon by-transfer, we want to turn off update-by-solve behaviour for these horizons. Hence, we wrap 
them in ExternalHorizon, which specializes the update! method to do nothing. 
Which cores own which scenarios are defined in db.pp_dist at any given time. 
"""
# TODO: No-data-version-of-non-master-horizon
function add_local_horizons(thiscore)
    db = get_local_db()
    horizons = get_horizons(db.input)
    d = Dict()
    for (scenario, ownercore) in db.pp_dist
        d[scenario] = Dict()
        for ((term, commodity), horizon) in horizons
            horizon = deepcopy(horizon)
            if ownercore != thiscore
                horizon = ExternalHorizon(horizon)    
            end
            d[scenario][(term, commodity)] = horizon
        end
    end
    db.horizon = d
    return
end

# TODO: Ensure N scenarioproblems that may change scenario (maybe use Vector instead of Dict)
function add_local_problems(thiscore)
    db = get_local_db()

    d = Dict()
    for (scenario, core) in db.pp_dist
        if core == thiscore
            d[scenario] = create_pp_problem(db.input, scenario)
        end
    end
    db.pp = d

    d = Dict()
    for subsystem in keys(db.ev_dist)
        for (scenario, core) in db.ev_dist[subsystem]
            if core == thiscore
                d[(subsystem, scenario)] = create_ev_problem(db.input, subsystem, scenario)
            end
        end
    end
    db.ev = d

    d = Dict()
    for (subsystem, core) in db.mp_dist
        if core == thiscore
            d[subsystem] = create_mp_problem(db.input, subsystem)
        end
    end
    db.mp = d

    d = Dict()
    for subsystem in keys(db.sp_dist)
        for (scenario, core) in db.sp_dist[subsystem]
            if core == thiscore
                d[(subsystem, scenario)] = create_sp_problem(db.input, subsystem, scenario)
            end
        end
    end
    db.sp = d

    if thiscore == db.cp_core
        db.cp = create_cp_problem(db.input)
    end
    return
end

function step_jules(output::AbstractJulESOutput, t, delta, stepnr)
    cores = get_cores(output)

    # do input models here
    # generate_scenarios
    # choose_scenarios

    @sync for core in cores
        @spawnat core solve_pp_problems(t, delta, stepnr, core)
    end

    @sync for core in cores
        @spawnat core solve_ev_problems(t, delta, stepnr, core)
    end

    @sync for core in cores
        @spawnat core solve_mp_problems(t, delta, stepnr, core)
    end

    @sync for core in cores
        @spawnat core solve_cp_problem(t, delta, stepnr, core)
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
function solve_pp_problems(T, t, delta, stepnr, thiscore)
    update_startstates_pp(stepnr)
    update_endstates_pp(stepnr) # only long, rest happens in solve
    solve_local_pp_problems(t)
    syncronize_horizons(thiscore)
    return
end

function solve_ev_problems(T, t, delta, stepnr, thiscore)
    update_startstates_ev(stepnr)
    update_endstates_ev()
    update_prices_ev()
    solve_local_ev_problems(t)
    return
end

function solve_mp_problems(T, t, delta, stepnr, thiscore)
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

function solve_cp_problem(T, t, delta, stepnr, thiscore)
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

function update_startstates_pp(stepnr)
    db = get_local_db()
    problems = db.pp
    if stepnr == 1
        startstates = get_startstates_pp(db.input)
    else
        if has_startstates_pp(db)
            startstates = get_startstates_pp(db)
        else
            startstates = get_startstates_pp_from_cp(db)
        end
    end
    for p in problems
        setstartstates!(p, startstates)
    end
end

function solve_local_pp_problems(t)
    db = get_local_db()
    for p in db.pp
        update!(p.long, t)
        solve!(p.long)
        # transfer states from long to med
        update!(p.med, t)
        solve!(p.med)
        # transfer states from med to short
        update!(p.short, t)
        solve!(p.short)
    end
end

function syncronize_horizons(thiscore)
    db = get_local_db()
    owner_scenarios = [s for (s, c) in db.pp_dist if c == thiscore]
    for ((scenario, term, commodity), horizon) in db.horizons
        if !(scenario in owner_scenarios)
            continue
        end
        changes = getchanges(horizon)
        if length(changes)
            @sync for (s, c) in db.pp_dist if !(s in owner_scenarios)
                @spawnat c transfer_horizon_changes(s, term, commodity, changes)
            end        
        end
    end
end

function transfer_horizon_changes(s, term, commodity, changes)
    db = get_local_db()
    h = db.horizons[(s, term, commodity)]
    setchanges!(h, changes)
    return
end

function update_startstates_ev(stepnr)

end
