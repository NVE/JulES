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
    modelobjects = TuLiPa.getmodelobjects(elements)
    add_PowerUpperSlack!(modelobjects)

    for (subix, core) in db.dist_mp # or get list of cuts from each core?
        future = @spawnat core get_lightcuts(subix)
        ret = fetch(future)
        if ret isa RemoteException
            throw(ret)
        end
        cuts = deepcopy(ret)

        # Change statevars so that they represents clearing version of objects
        for i in 1:length(cuts.statevars)
            varin = TuLiPa.getvarin(cuts.statevars[i])
            (varid, varix) = TuLiPa.getvarout(cuts.statevars[i])

            newt = TuLiPa.getnumperiods(TuLiPa.gethorizon(modelobjects[varid]))
            cuts.statevars[i] = TuLiPa.StateVariableInfo(varin, (varid, newt))
        end
        cutid = TuLiPa.getid(cuts)
        modelobjects[cutid] = cuts
    end

    probmethod = parse_methods(settings["problems"]["clearing"]["solver"])
    prob = TuLiPa.buildprob(probmethod, modelobjects)

    db.cp = ClearingProblem(prob, Dict{String, Float64}(), Dict())

    return
end

function solve_cp(t, delta, stepnr, skipmed)
    db = get_local_db()

    if db.core_cp == db.core
        timing = db.output.timing_cp
        timing[stepnr, 3] = @elapsed begin
            update_startstates_cp(db, stepnr, t)
            update_cuts(db, skipmed)
            update_nonstoragestates_cp(db)
            update_statedependent_cp(db, stepnr, t)
            timing[stepnr, 1] = @elapsed TuLiPa.update!(db.cp.prob, t)
            set_minstoragevalue!(db.cp.prob, minstoragevaluerule)
            timing[stepnr, 2] = @elapsed TuLiPa.solve!(db.cp.prob)
            get_startstates!(db.cp.prob, db.input.dataset["detailedrescopl"], db.input.dataset["enekvglobaldict"], db.cp.endstates)
        end
    end
end

# Util functions for solve_cp ----------------------------------------------------------------------------------
function minstoragevaluerule(storage::TuLiPa.Storage)
    commodity = TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(TuLiPa.getbalance(storage))))
    if commodity == "Hydro"
        return 0.001
    end
    return 0.0
 end

function set_minstoragevalue!(problem::TuLiPa.Prob, costrule::Function)
    for modelobject in TuLiPa.getobjects(problem)
        if modelobject isa TuLiPa.Storage
            id = TuLiPa.getid(modelobject)
            balance = TuLiPa.getbalance(modelobject)
            horizon = TuLiPa.gethorizon(balance)
            T = TuLiPa.getnumperiods(horizon)
            coeff = TuLiPa.getobjcoeff(problem, id, T)
            cost = costrule(modelobject)
            newcoeff = min(-cost, coeff)
            if !(coeff ≈ newcoeff)
                TuLiPa.setobjcoeff!(problem, id, T, newcoeff)
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

    get_statedependentprod(settings["problems"]["clearing"]) && TuLiPa.statedependentprod!(db.cp.prob, db.startstates, init=init)
    get_statedependentpump(settings["problems"]["clearing"]) && TuLiPa.statedependentpump!(db.cp.prob, db.startstates)

    # Headlosscosts
    if get_headlosscost(settings["problems"]["clearing"])
        for (_subix, _core) in db.dist_mp
            future = @spawnat _core get_headlosscost_data_from_mp(_subix, t)

            ret = fetch(future)
            if ret isa RemoteException
                throw(ret)
            end
            headlosscost_data = ret

            for (resid, headlosscost, T) in headlosscost_data
                obj = find_obj_by_id(TuLiPa.getobjects(db.cp.prob), resid)
                T = TuLiPa.getnumperiods(TuLiPa.gethorizon(obj))
                
                TuLiPa.setobjcoeff!(db.cp.prob, resid, T, headlosscost)
            end
        end
    end
end

function get_headlosscost_data_from_mp(subix, t) # TODO: get method from config
    db = get_local_db()

    mp = db.mp[subix]

    return TuLiPa.get_headlosscost_data(TuLiPa.ReservoirCurveSlopeMethod(), mp.prob, t)
end

function update_nonstoragestates_cp(db)
    scenix = 1 # which scenario to use?

    for (_scenix, _core) in db.dist_ppp
        if scenix == _scenix
            future = @spawnat _core get_nonstoragestates_short(scenix)

            ret = fetch(future)
            if ret isa RemoteException
                throw(ret)
            end

            nonstoragestates_short = ret
            TuLiPa.setoutgoingstates!(db.cp.prob, nonstoragestates_short)
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

            ret = fetch(future)
            if ret isa RemoteException
                throw(ret)
            end

            (cutid, constants, slopes) = ret

            cuts_cp = find_obj_by_id(TuLiPa.getobjects(db.cp.prob), cutid)
            cuts_cp.constants = constants
            cuts_cp.slopes = slopes

            TuLiPa.updatecuts!(db.cp.prob, cuts_cp)
        end
    end
end

function get_lightcuts(subix)
    db = get_local_db()
    cuts = db.mp[subix].cuts
    return TuLiPa.getlightweightself(cuts)
end

function get_cutsdata(subix)
    db = get_local_db()

    cuts = db.mp[subix].cuts
    return (cuts.id, cuts.constants, cuts.slopes)
end

function update_startstates_cp(db, stepnr, t)
    if stepnr == 1 # Non storage startstates free in first step
        set_startstates!(db.cp.prob, TuLiPa.getstorages(TuLiPa.getobjects(db.cp.prob)), db.startstates)
    else
        set_startstates!(db.cp.prob, TuLiPa.getobjects(db.cp.prob), db.startstates)
    end
end

# Util functions create_cp ------------------------------------------------------------------------------------


# TODO: Rename to update_startstates
function get_startstates!(clearing::TuLiPa.Prob, detailedrescopl::Dict, enekvglobaldict::Dict, startstates::Dict{String, Float64})
    startstates_ = get_states(TuLiPa.getobjects(clearing))
    TuLiPa.getoutgoingstates!(clearing, startstates_)
    
    for var in keys(startstates_)
        value = round(startstates_[var], digits=10) # avoid approx 0 negative values, ignored by solvers so no problem?
        startstates[TuLiPa.getinstancename(first(TuLiPa.getvarout(var)))] = value
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