

struct PricePrognosisProblem
    long
    med
    short
    long_prices
    med_prices
    short_prices
    short_nonstoragestates
end

function setstartstates!(p::PricePrognosisProblem, startstates)
end

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

struct ClearingProblem
    prob
    endstates
end

function make_obj(elements::Vector{DataElement}, hydro_horizon::Horizon, power_horizon::Horizon; validate::Bool=false)
    elements1 = copy(elements)

    set_horizon!(elements1, "Power", power_horizon)
    set_horizon!(elements1, "Battery", power_horizon)
    set_horizon!(elements1, "Hydro", hydro_horizon)

    modelobjects = getmodelobjects(elements1; validate=validate)

    return modelobjects
end

function create_ppp(db::LocalDB, scenarioix)
    probmethods
    (aggzone, aggsuplyn, removestoragehours, residualarealist) = simplifyinputs

    lhh = db.horizons[(scenarioix, "Long", "Hydro")]
    lph = db.horizons[(scenarioix, "Long", "Power")]
    longobjects = make_obj(db.input.progelements, lhh, lph)
    longprobmethod = 

    simplify!(longobjects; aggzone=aggzone, aggsupplyn=aggsuplyn, removestoragehours=removestoragehours, residualarealist=residualarealist)
    addPowerUpperSlack!(longobjects)
    longprob = buildprob(probmethods[1], longobjects)


end

function create_evp(db::LocalDB, scenarioix, subsystem)
end

function create_mp(db::LocalDB, subsystem) 
end

function create_sp(db::LocalDB, scenarioix, subsystem)
end

function create_cp(db::LocalDB)
end
