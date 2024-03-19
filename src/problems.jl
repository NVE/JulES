# TODO: Support different terms? (med or long)
struct EndValueProblem
    prob
end

# TODO: Support different terms? (med or long)
struct ScenarioProblem
    prob
end

struct MasterProblem
    prob
end

function create_evp(db::LocalDB, scenarioix, subsystem)
end

function create_mp(db::LocalDB, subsystem) 
end

function create_sp(db::LocalDB, scenarioix, subsystem)
end

