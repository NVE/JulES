get_prob(cp::ClearingProblem) = cp.prob
get_endstates(cp::ClearingProblem) = cp.endstates

get_startstates_from_cp() = get_endstates(get_local_db().cp)

function create_cp()
    db = get_local_db()
    settings = get_settings(db)

    elements = copy(get_elements(db.input))
    horizons = get_horizons(db.input)

    for ((term, commodity), horizon) in horizons
        if term == ClearingTermName
            set_horizon!(elements, commodity, horizon)
            if commodity == "Power"
                set_horizon!(elements, "Battery", horizon)
            end
        end
    end
    modelobjects = getmodelobjects(elements)
    add_PowerUpperSlack!(modelobjects)

    for (subix, core) in db.dist_mp # or get list of cuts from each core?
        future = @spawnat core get_lightcuts(subix)
        cuts = deepcopy(fetch(future))

        # Change statevars so that they represents clearing version of objects
        for i in 1:length(cuts.statevars)
            varin = getvarin(cuts.statevars[i])
            (varid, varix) = getvarout(cuts.statevars[i])

            newt = getnumperiods(gethorizon(modelobjects[varid]))
            cuts.statevars[i] = StateVariableInfo(varin, (varid, newt))
        end
        cutid = getid(cuts)
        modelobjects[cutid] = cuts
    end

    probmethod = parse_methods(settings["problems"]["clearing"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.cp = ClearingProblem(prob, Dict{String, Float64}(), Dict())

    return
end

function solve_cp(t, stepnr, skipmed)
    db = get_local_db()

    if db.core_cp == db.core
        timing = db.output.timing_cp
        timing[stepnr, 3] = @elapsed begin
            update_startstates_cp(db, stepnr, t)
            update_cuts(db, skipmed)
            update_nonstoragestates_cp(db)
            update_statedependent_cp(db, stepnr, t)
            timing[stepnr, 1] = @elapsed update!(db.cp.prob, t)
            set_minstoragevalue!(db.cp.prob, minstoragevaluerule)
            timing[stepnr, 2] = @elapsed solve!(db.cp.prob)
            get_startstates!(db.cp.prob, db.input.dataset["detailedrescopl"], db.input.dataset["enekvglobaldict"], db.cp.endstates)
        end
    end
end

# Util functions for solve_cp ----------------------------------------------------------------------------------
function minstoragevaluerule(storage::Storage)
    commodity = getinstancename(getid(getcommodity(getbalance(storage))))
    if commodity == "Hydro"
        return 0.001
    end
    return 0.0
 end

function set_minstoragevalue!(problem::Prob, costrule::Function)
    for modelobject in getobjects(problem)
        if modelobject isa Storage
            id = getid(modelobject)
            balance = getbalance(modelobject)
            horizon = gethorizon(balance)
            T = getnumperiods(horizon)
            coeff = getobjcoeff(problem, id, T)
            cost = costrule(modelobject)
            newcoeff = min(-cost, coeff)
            if !(coeff ≈ newcoeff)
                setobjcoeff!(problem, id, T, newcoeff)
            end
        end
    end
    return
end

function update_statedependent_cp(db, stepnr, t)
    settings = get_settings(db)

    # Statedependent prod and pumping
    init = false
    if stepnr == 1
        init = true
    end

    get_statedependentprod(settings["problems"]["clearing"]) && statedependentprod!(db.cp.prob, db.startstates, init=init)
    get_statedependentpump(settings["problems"]["clearing"]) && statedependentpump!(db.cp.prob, db.startstates)

    # Headlosscosts
    if get_headlosscost(settings["problems"]["clearing"])
        for (_subix, _core) in db.dist_mp
            future = @spawnat _core get_headlosscost_data_from_mp(_subix, t)
            headlosscost_data = fetch(future)

            for (resid, headlosscost, T) in headlosscost_data
                obj = find_obj_by_id(getobjects(db.cp.prob), resid)
                T = getnumperiods(gethorizon(obj))
                
                setobjcoeff!(db.cp.prob, resid, T, headlosscost)
            end
        end
    end
end

function get_headlosscost_data_from_mp(subix, t) # TODO: get method from config
    db = get_local_db()

    mp = db.mp[subix]

    return get_headlosscost_data(ReservoirCurveSlopeMethod(), mp.prob, t)
end

function update_nonstoragestates_cp(db)
    scenix = 1 # which scenario to use?

    for (_scenix, _core) in db.dist_ppp
        if scenix == _scenix
            future = @spawnat _core get_nonstoragestates_short(scenix)
            nonstoragestates_short = fetch(future)
            setoutgoingstates!(db.cp.prob, nonstoragestates_short)
        end
    end
end

function get_nonstoragestates_short(scenix)
    db = get_local_db()

    return db.ppp[scenix].nonstoragestates_short
end

function update_cuts(db, skipmed)
    for (_subix, _core) in db.dist_mp
        if skipmed_check(_subix, skipmed)
            future = @spawnat _core get_cutsdata(_subix)
            (cutid, constants, slopes) = fetch(future)

            cuts_cp = find_obj_by_id(getobjects(db.cp.prob), cutid)
            cuts_cp.constants = constants
            cuts_cp.slopes = slopes

            updatecuts!(db.cp.prob, cuts_cp)
        end
    end
end

function get_lightcuts(subix)
    db = get_local_db()
    cuts = db.mp[subix].cuts
    return getlightweightself(cuts)
end

function get_cutsdata(subix)
    db = get_local_db()

    cuts = db.mp[subix].cuts
    return (cuts.id, cuts.constants, cuts.slopes)
end

function update_startstates_cp(db, stepnr, t)
    if stepnr == 1 # Non storage startstates free in first step
        set_startstates!(db.cp.prob, getstorages(getobjects(db.cp.prob)), db.startstates)
    else
        set_startstates!(db.cp.prob, getobjects(db.cp.prob), db.startstates)
    end
end

# Util functions create_cp ------------------------------------------------------------------------------------


# TODO: Rename to update_startstates
function get_startstates!(clearing::Prob, detailedrescopl::Dict, enekvglobaldict::Dict, startstates::Dict{String, Float64})
    startstates_ = get_states(getobjects(clearing))
    getoutgoingstates!(clearing, startstates_)
    
    for var in keys(startstates_)
        value = round(startstates_[var], digits=10) # avoid approx 0 negative values, ignored by solvers so no problem?
        startstates[getinstancename(first(getvarout(var)))] = value
    end

    for area in Set(values(detailedrescopl))
        startstates["Reservoir_" * area * "_hydro_reservoir"] = 0.0
    end
    for res in keys(detailedrescopl)
        resname = "Reservoir_" * res
        areaname = "Reservoir_" * detailedrescopl[res] * "_hydro_reservoir"
        startstates[areaname] += startstates[resname] * enekvglobaldict[res]
    end

    # Avoid reservoirs being filled more than max, gives infeasible solution
    # - If aggregated reservoir capacity is lower than the sum capacities
    # - If reservoir is full in model, numerical tolerance can bring variable value slightly over cap
    # - TODO: Add warning/logging if this happens
    for resname in keys(startstates)
        resmax = resname * "_max"
        if haskey(startstates, resmax)
            if startstates[resname] > startstates[resmax]
                startstates[resname] = startstates[resmax]
            end
        end
    end
end