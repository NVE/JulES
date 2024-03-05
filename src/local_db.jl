
# Define simple local-database API 

const _LOCAL_DB_NAME = :_local_db

# TODO: Complete this
mutable struct LocalDB
    input
    horizon

    pp
    ev
    mp
    sp
    cp

    pp_dist
    ev_dist
    mp_dist
    sp_dist
    cp_core

    cp_time_solve
    cp_time_update
    cp_time_cuts
    cp_time_startstates
    cp_time_endstates

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
    Core.eval(Main, Expr(:(=), _LOCAL_DB_NAME, LocalDB())
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