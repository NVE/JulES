"""
In-memory-process-local-database for JulES.

We utilize some of Julia's reflection (fieldnames, getfield) 
and metaprogramming (eval, Expr) capabilities to implement 
a database object in the global scope of a Julia process.
We leared this technique from this package: 
https://github.com/ChrisRackauckas/ParallelDataTransfer.jl

The database is just a struct with slots for all the different
things we want to store on cores while running JulES. 
    
Some slots hold read-only data that is identical on all cores (input)

Some slots hold stateful objects that are managed on a particular core
while other cores hold synced copies (horizons)

Some slots hold optimization problems where a problem resides
at exactly one core at a time, but may be moved to antother core 
over time (ppp, evp, mp, sp, cp)

Some slots hold info about on which cores problems are stored. This
enables JulES to locate and use any problem from any core, which is
very useful e.g. when solving stochastic optimization problems where
the master problem and the different subproblems reside on different
cores, or e.g. when we need to get results from a problem residing on
one core in order to use it as input to another problem residing on a 
different core (ppp_dist, evp_dist, mp_dist, sp_dist, cp_core)

Many of the slots contain timing data.
"""

const _LOCAL_DB_NAME = :_local_db

# TODO: Complete this (add more timings and maybe other stuff)
mutable struct LocalDB
    input::AbstractJulESInput
    horizons::Dict{ScenarioTermCommodity, Horizon}

    ppp::Dict{Scenario, PricePrognosisProblem}
    evp::Dict{ScenarioSubsystem, EndValueProblem}
    sp::Dict{ScenarioSubsystem, ScenarioProblem}
    mp::Dict{Subsystem, MasterProblem}
    cp::Union{Nothing, ClearingProblem}

    ppp_dist::Vector{ScenarioCore}
    evp_dist::Vector{ScenarioSubsystemCore}
    mp_dist::Vector{SubsystemCore}
    sp_dist::Vector{ScenarioSubsystemCore}
    cp_core::Int

    cp_time_solve::Float64
    cp_time_update::Float64
    cp_time_cuts::Float64
    cp_time_startstates::Float64
    cp_time_endstates::Float64

    # TODO: Fix by handle types, will fail now
    function LocalDB()
        n = length(fieldnames(typeof(LocalDB)))
        return new([nothing for _ in 1:n]...)
    end
end

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