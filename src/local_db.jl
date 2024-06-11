"""
In-memory process-local database for JulES.

We utilize some of Julia's reflection (names, getfield) 
and metaprogramming (eval, Expr) capabilities to implement 
a database object in the global scope of a Julia process.
We learned about this technique from this package: 
https://github.com/ChrisRackauckas/ParallelDataTransfer.jl

The database is just a struct with fields for all the different
things we want to store on cores while running JulES. 
    
The input field holds read-only data that is identical on all cores

The horizons field holds horizon objects, where some of which are stateful 
and managed on a particular core, while others are syncronized copies of
horizons managed on other cores

Some fields hold optimization problems. A problem resides
at exactly one core at a time, but may be moved to antother core 
between time steps. See fields ppp, evp, mp, sp and cp.

Some fields hold info about on which cores problems are stored. 
This also enables transfer of data between e.g. optimizion problems
residing on different cores. See fields dist_ppp, dist_evp, dist_mp, 
dist_sp and core_cp.

Many of the fields contain timing data. This is useful both for results
and to inform dynamic load balancer.

The div field holds a Dict object. Possible extentions of JulES may use
this field to store data.
"""

const _LOCAL_DB_NAME = :_local_db

# TODO: Complete this (add more timings and maybe other stuff)
mutable struct LocalDB
    core::CoreId
    
    input::Union{Nothing, AbstractJulESInput}
    output::Union{Nothing, AbstractJulESOutput}
    horizons::Dict{Tuple{ScenarioIx, TermName, CommodityName}, TuLiPa.Horizon}

    dummyobjects::Tuple
    dummyobjects_ppp::Tuple # TODO: Move dummyobjects, scenariogeneration and subsystems to io?

    startstates::Dict{String, Float64}

    subsystems::Vector{AbstractSubsystem}
    subsystems_evp::Vector{Tuple{SubsystemIx, AbstractSubsystem}}
    subsystems_stoch::Vector{Tuple{SubsystemIx, AbstractSubsystem}}

    scenmod_sim::AbstractScenarioModellingMethod
    scenmod_ppp::AbstractScenarioModellingMethod
    scenmod_evp::AbstractScenarioModellingMethod
    scenmod_stoch::AbstractScenarioModellingMethod

    ifm::Dict{String, Union{AbstractInflowModel, Nothing}}
    ppp::Dict{ScenarioIx, PricePrognosisProblem}
    prices_ppp::Dict{Tuple{ScenarioIx, TermName, TuLiPa.Id}, Tuple{Int, Vector{Float64}}}
    evp::Dict{Tuple{ScenarioIx, SubsystemIx}, EndValueProblem}
    mp::Dict{SubsystemIx, MasterProblem}
    sp::Dict{Tuple{ScenarioIx, SubsystemIx}, ScenarioProblem}
    cp::Union{Nothing, ClearingProblem}

    dist_ifm::Vector{Tuple{String, CoreId}}
    dist_ppp::Vector{Tuple{ScenarioIx, CoreId}}
    dist_evp::Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}
    dist_mp::Vector{Tuple{ScenarioIx, CoreId}}
    dist_sp::Vector{Tuple{ScenarioIx, SubsystemIx, CoreId}}
    core_cp::CoreId

    # time_cp::Float64

    div::Dict

    ifm_output::Dict{String, Dict{ScenarioIx, Tuple{Vector{DateTime}, Vector{Float64}}}}
    ifm_derived::Dict{String, Dict{ScenarioIx, Tuple{Vector{DateTime}, Vector{Float64}}}}

    function LocalDB()
        return new(
            -1,
        
            nothing,   # input
            nothing,   # output
            Dict{Tuple{ScenarioIx, TermName, CommodityName}, TuLiPa.Horizon}(),   # horizons

            (),   # dummyobjects
            (),   # dummyobjects_ppp

            Dict{String, Float64}(),    # startstates

            AbstractSubsystem[],       # subsystems
            AbstractSubsystem[],       # subsystems_evp
            AbstractSubsystem[],       # subsystems_stoch

            NothingScenarioModellingMethod(), # scenmod_sim
            NothingScenarioModellingMethod(), # scenmod_ppp
            NothingScenarioModellingMethod(), # scenmod_evp
            NothingScenarioModellingMethod(), # scenmod_stoch

            Dict(), # ifm
            Dict{ScenarioIx, PricePrognosisProblem}(),                               # ppp
            Dict{Tuple{ScenarioIx, TermName, TuLiPa.Id}, Tuple{Int, Vector{Float64}}}(), # prices_ppp
            Dict{Tuple{ScenarioIx, SubsystemIx}, EndValueProblem}(),                 # evp
            Dict{SubsystemIx, MasterProblem}(),                                      # mp
            Dict{Tuple{ScenarioIx, SubsystemIx}, ScenarioProblem}(),                 # sp
            nothing,                                                                 # cp

            [], # dist_ifm
            Tuple{ScenarioIx, CoreId}[],                # dist_ppp
            Tuple{ScenarioIx, SubsystemIx, CoreId}[],   # dist_evp
            Tuple{SubsystemIx, CoreId}[],               # dist_mp
            Tuple{ScenarioIx, SubsystemIx, CoreId}[],   # dist_sp
            -1,   # core_cp

            Dict(),   # div

            Dict(),   # ifm_output
            Dict(),   # ifm_weighted

        )
    end
end

get_core(db::LocalDB) = db.core
get_input(db::LocalDB) = db.input
get_output(db::LocalDB) = db.output
get_horizons(db::LocalDB) = db.horizons
get_dummyobjects(db::LocalDB) = db.dummyobjects
get_dummyobjects_ppp(db::LocalDB) = db.dummyobjects_ppp
get_startstates(db::LocalDB) = db.startstates
get_subsystems(db::LocalDB) = db.subsystems
get_subsystems_evp(db::LocalDB) = db.subsystems_evp
get_subsystems_stoch(db::LocalDB) = db.subsystems_stoch
get_scenmod_sim(db::LocalDB) = db.scenmod_sim
get_scenmod_ppp(db::LocalDB) = db.scenmod_ppp
get_scenmod_evp(db::LocalDB) = db.scenmod_evp
get_scenmod_stoch(db::LocalDB) = db.scenmod_stoch
get_ppp(db::LocalDB) = db.ppp
get_evp(db::LocalDB) = db.evp
get_mp(db::LocalDB) = db.mp
get_sp(db::LocalDB) = db.sp
get_cp(db::LocalDB) = db.cp
get_dist_ppp(db::LocalDB) = db.dist_ppp
get_dist_evp(db::LocalDB) = db.dist_evp
get_dist_mp(db::LocalDB) = db.dist_mp
get_dist_sp(db::LocalDB) = db.dist_sp
get_core_cp(db::LocalDB) = db.core_cp
get_div(db::LocalDB) = db.div

get_ifm(db::LocalDB) = db.ifm
get_dist_ifm(db::LocalDB) = db.dist_ifm
get_ifm_output(db::LocalDB) = db.ifm_output
get_ifm_derived(db::LocalDB) = db.ifm_derived

get_ifm_weights(db::LocalDB) = get_ifm_weights(get_input(db))
get_ifm_normfactors(db::LocalDB) = get_ifm_normfactors(get_input(db))
get_ifm_elements(db::LocalDB) = get_ifm_elements(get_input(db))

get_cores(db::LocalDB) = get_cores(get_input(db))
get_dataset(db::LocalDB) = get_dataset(get_input(db))
get_mainconfig(db::LocalDB) = get_mainconfig(get_input(db))
get_settings(db::LocalDB) = get_settings(get_input(db))
get_datayear(db::LocalDB) = get_datayear(get_input(db))
get_weatheryear(db::LocalDB) = get_weatheryear(get_input(db))
get_onlysubsystemmodel(db::LocalDB) = get_onlysubsystemmodel(get_input(db))
get_steps(db::LocalDB) = get_steps(get_input(db))
get_steplength(db::LocalDB) = get_steplength(get_input(db))
get_simstarttime(db::LocalDB) = get_simstarttime(get_input(db))
get_scenmod_data(db::LocalDB) = get_scenmod_data(get_input(db))
get_numscen_data(db::LocalDB) = get_numscen_data(get_input(db))
get_tnormaltype(db::LocalDB) = get_tnormaltype(get_input(db))
get_tphaseintype(db::LocalDB) = get_tphaseintype(get_input(db))
get_phaseinoffset(db::LocalDB) = get_phaseinoffset(get_input(db))
get_phaseindelta(db::LocalDB) = get_phaseindelta(get_input(db))
get_phaseinsteps(db::LocalDB) = get_phaseinsteps(get_input(db))

function create_local_db()
    if (_LOCAL_DB_NAME in names(Main)) 
        db = getfield(Main, _LOCAL_DB_NAME)
        isnothing(db) || error("$_LOCAL_DB_NAME already exists")
    end
    Core.eval(Main, Expr(:(=), _LOCAL_DB_NAME, LocalDB()))
    return
end

function get_local_db()
    (_LOCAL_DB_NAME in names(Main)) || error("$_LOCAL_DB_NAME has not been created")
    db = getfield(Main, _LOCAL_DB_NAME)
    isnothing(db) && error("$_LOCAL_DB_NAME has been freed")
    return db::LocalDB
end

function free_local_db()
    Core.eval(Main, Expr(:(=), _LOCAL_DB_NAME, nothing))
    return
end