function create_mp(db::LocalDB, subix::SubsystemIx)
    scenix = 1 # TODO: Which scenario should be represented in the master problem? Not important due to phasein?
    settings = get_settings(db)

    # TODO: Not use ifm here as in clearnig?
    modelobjects, horizons = make_modelobjects_stochastic(db, scenix, subix, true)

    maxcuts = settings["problems"]["stochastic"]["maxcuts"] # preallocate fixed number of cuts, no cut selection
    lb = settings["problems"]["stochastic"]["lb"] # lower bound of the future value in the first iteration
    numscen_stoch = get_numscen_stoch(db.input)
    cutobjects = getcutobjects(modelobjects)
    cuts = initialize_cuts!(modelobjects, cutobjects, maxcuts, lb, numscen_stoch)
    states = get_states(cutobjects)
    cuts.statevars = [var for (var, value) in states] # TODO: Find better way of getting same order in states and cuts.statevars

    probmethod = parse_methods(settings["problems"]["stochastic"]["master"]["solver"])
    prob = TuLiPa.buildprob(probmethod, modelobjects)

    div = Dict()
    div[MainTiming] = zeros(5)
    if has_result_storagevalues(settings)
        if has_headlosscost(settings["problems"]["stochastic"]["master"])
            num_storagevalues = get_numscen_stoch(db.input)*2 + 2 # scenarios + master operative + master operative after headlosscost adjustment
        else
            num_storagevalues = get_numscen_stoch(db.input)*2 + 1 # scenarios + master operative 
        end
        div[StorageValues] = zeros(num_storagevalues, length(states))
    else
        div[StorageValues] = nothing
    end

    db.mp[subix] = MasterProblem(prob, cuts, states, div)

    return
end

function create_sp(db::LocalDB, scenix::ScenarioIx, subix::SubsystemIx)
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    modelobjects, horizons = make_modelobjects_stochastic(db, scenix, subix, false)
    states = get_states(modelobjects) # different order than mp.cuts.statevars, so only use length

    probmethod = parse_methods(settings["problems"]["stochastic"]["subs"]["solver"])
    prob = TuLiPa.buildprob(probmethod, modelobjects)

    div = Dict()
    div[MainTiming] = zeros(3)

    db.sp[(scenix, subix)] = ScenarioProblem(prob, horizons, zeros(length(states)), -1.0, div)

    return
end

function solve_stoch(t, stepnr, skipmed)
    db = get_local_db()
    settings = get_settings(db)

    for (subix, core) in db.dist_mp
        if core == db.core
            subsystem = db.subsystems[subix]
            if skipmed_check(subsystem, skipmed)
                mp = db.mp[subix]
                maintiming = mp.div[MainTiming]

                maintiming[4] = @elapsed begin
                    has_headlosscost(settings["problems"]["stochastic"]["master"]) && TuLiPa.resetheadlosscosts!(mp.prob)
                    update_probabilities(mp.cuts, db.scenmod_stoch) # TODO: Add possibility for scenario modelling per subsystem
                    set_startstates!(mp.prob, TuLiPa.getstorages(TuLiPa.getobjects(mp.prob)), db.startstates)
                    update_prices_mp(stepnr, subix, subsystem)
                    update_statedependent_mp(stepnr, subsystem, mp.prob, db.startstates, settings)
                end
                maintiming[1] = @elapsed TuLiPa.update!(mp.prob, t)

                @sync for _core in get_cores(db)
                    @spawnat _core update_sps(t, stepnr, subix)
                end
                
                solve_benders(stepnr, subix)
                maintiming[4] += @elapsed save_storagevalues(mp.prob, mp.cuts, mp.div[StorageValues])
                maintiming[3] = @elapsed final_solve_mp(t, mp.prob, mp.cuts, mp.div[StorageValues], settings)
            end
        end
    end

end

function update_sps(t, stepnr, subix)
    db = get_local_db()
    subsystem = db.subsystems[subix]

    for (_scenix, _subix, _core) in db.dist_sp
        if (_core == db.core) && (subix == _subix)
            maintiming = db.sp[(_scenix, subix)].div[MainTiming]
            maintiming[3] = @elapsed begin
                update_horizons_sp(_scenix, subix, subsystem)
                update_endconditions_sp(_scenix, subix, t)
                perform_scenmod_sp(_scenix, subix)
                update_prices_sp(stepnr, _scenix, subix, subsystem)
            end
            maintiming[1] = @elapsed update_sp(t, _scenix, subix)
        end
    end
end

# Util functions for solve_stoch ----------------------------------------------------------------------------------------------

function final_solve_mp(t::TuLiPa.ProbTime, prob, cuts, storagevalues, settings)
    if has_headlosscost(settings["problems"]["stochastic"]["master"])
        TuLiPa.updateheadlosscosts!(TuLiPa.ReservoirCurveSlopeMethod(), prob, t)
        TuLiPa.solve!(prob)
        !isnothing(storagevalues) && final_save_storagevalues(prob, cuts, storagevalues)
    end
end 

function final_save_storagevalues(prob, cuts, storagevalues)
    for (j, statevar) in enumerate(cuts.statevars) # master / operative water values after headlosscost
        obj = get_obj_from_id(cuts.objects, first(TuLiPa.getvarout(statevar))) # TODO: OK to assume objid = varoutid?
        balance = TuLiPa.getbalance(obj)

        storagevalues[length(cuts.probabilities)*2 + 2, j] = TuLiPa.getcondual(prob, TuLiPa.getid(balance), TuLiPa.getnumperiods(TuLiPa.gethorizon(balance)))
        if haskey(balance.metadata, TuLiPa.GLOBALENEQKEY)
            storagevalues[length(cuts.probabilities)*2 + 2, j] = storagevalues[length(cuts.probabilities)*2 + 2, j] / balance.metadata[TuLiPa.GLOBALENEQKEY]
        end
    end
    return
end


function save_storagevalues(prob, cuts, storagevalues)
    if !isnothing(storagevalues)
        for (j, statevar) in enumerate(cuts.statevars)
            obj = get_obj_from_id(cuts.objects, first(TuLiPa.getvarout(statevar))) # TODO: OK to assume objid = varoutid?
            balance = TuLiPa.getbalance(obj)
            for i in 1:length(cuts.probabilities) # scenario storage values
                minslope = 0
                maxslope = -1e9
                for k in 1:TuLiPa.getnumcuts(cuts)
                    val = cuts.scenslopes[i,k,j]
                    minslope = min(minslope, val)
                    maxslope = max(maxslope, val)
                end
                if haskey(balance.metadata, TuLiPa.GLOBALENEQKEY)
                    minslope = minslope / balance.metadata[TuLiPa.GLOBALENEQKEY]
                    maxslope = maxslope / balance.metadata[TuLiPa.GLOBALENEQKEY]
                end
                storagevalues[(i-1)*2+1, j] = minslope
                storagevalues[(i-1)*2+2, j] = maxslope
            end

            # master / operative water values
            storagevalues[length(cuts.probabilities)*2 + 1, j] = TuLiPa.getcondual(prob, TuLiPa.getid(balance), TuLiPa.getnumperiods(TuLiPa.gethorizon(balance)))
            if haskey(balance.metadata, TuLiPa.GLOBALENEQKEY)
                storagevalues[length(cuts.probabilities)*2 + 1, j] = storagevalues[length(cuts.probabilities)*2 + 1, j] / balance.metadata[TuLiPa.GLOBALENEQKEY]
            end
        end
    end
    return
end

"""
Two possibilities to warm start benders solve:
- Start with solving master problem if cuts from previous step are available. NB! we reuse cuts although scenarios can have changed
- Start with solving subproblems with start reservoirs from db.startstates if cuts from previous step are NOT available
"""
function solve_benders(stepnr, subix)
    db = get_local_db()
    settings = get_settings(db)

    mp = db.mp[subix]
    maintiming = mp.div[MainTiming]

    count = 0
    cutreuse = true
    if stepnr == 1
        cutreuse = false
    end
    ub = 0.0
    lb = mp.cuts.lower_bound
    reltol = settings["problems"]["stochastic"]["reltol"] # relative tolerance

    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1) && count < 15
        maintiming[4] += @elapsed begin
            if (count == 1 && cutreuse)
                TuLiPa.updatecuts!(mp.prob, mp.cuts)
            elseif count != 0
                TuLiPa.updatelastcut!(mp.prob, mp.cuts)
            end
        end

        maintiming[2] += @elapsed begin
            if TuLiPa.getnumcuts(mp.cuts) != 0
                count == 0 && TuLiPa.setwarmstart!(mp.prob, false)
                if cutreuse
                    try
                        TuLiPa.solve!(mp.prob)
                        count == 0 && TuLiPa.clearcuts!(mp.cuts)
                        count += 1
                    catch
                        count == 0 && println("Retrying first iteration without cuts from last time step")
                        count > 0 && println("Restarting iterations without cuts from last time step")
                        TuLiPa.clearcuts!(mp.prob, mp.cuts)
                        cutreuse = false
                        count = 0
                    end
                else
                    TuLiPa.solve!(mp.prob)
                    count += 1
                end
                count == 0 && TuLiPa.setwarmstart!(mp.prob, true)
                lb = TuLiPa.getvarvalue(mp.prob, TuLiPa.getfuturecostvarid(mp.cuts), 1)
                TuLiPa.getoutgoingstates!(mp.prob, mp.states)
            end
        end

        maintiming[4] += @elapsed begin
            if TuLiPa.getnumcuts(mp.cuts) != 0

            end
        end

        futures = []
        @sync for (_scenix, _subix, _core) in db.dist_sp
            if _subix == subix
                if TuLiPa.getnumcuts(mp.cuts) != 0
                    f = @spawnat _core solve_sp(_scenix, _subix, mp.states)
                else
                    f = @spawnat _core solve_sp_with_startreservoirs(_scenix, _subix, mp.states)
                end
                push!(futures, f)
            end
        end

        maintiming[4] += @elapsed begin
            ub = 0.0
            cutix = TuLiPa.getcutix(mp.cuts) + 1
            if cutix > TuLiPa.getmaxcuts(mp.cuts)
                cutix = 1
            end
            for future in futures
                scenix, objectivevalue, scenslopes, scenconstant = fetch(future)

                ub += objectivevalue*mp.cuts.probabilities[scenix]
                mp.cuts.scenslopes[scenix, cutix, :] .= scenslopes
                mp.cuts.scenconstants[scenix, cutix] = scenconstant
            end
        
            TuLiPa.updatecutparameters!(mp.prob, mp.cuts)
        end
    end
    maintiming[5] = count
    return
end

function solve_sp(scenix, subix, states)
    db = get_local_db()

    sp = db.sp[(scenix, subix)]
    maintiming = sp.div[MainTiming]

    maintiming[3] += @elapsed TuLiPa.setingoingstates!(sp.prob, states)
    maintiming[2] += @elapsed TuLiPa.solve!(sp.prob)
    maintiming[3] += @elapsed begin
        get_scencutparameters!(sp, states)

        objectivevalue = TuLiPa.getobjectivevalue(sp.prob)
        scenslopes = sp.scenslopes
        scenconstant = sp.scenconstant
    end
    return (scenix, objectivevalue, scenslopes, scenconstant)
end

function solve_sp_with_startreservoirs(scenix, subix, states)
    db = get_local_db()

    sp = db.sp[(scenix, subix)]
    maintiming = sp.div[MainTiming]

    maintiming[3] += @elapsed set_startstates!(sp.prob, TuLiPa.getstorages(TuLiPa.getobjects(sp.prob)), db.startstates)
    maintiming[2] += @elapsed TuLiPa.solve!(sp.prob)
    maintiming[3] += @elapsed begin
        get_scencutparameters!(sp, states) # TODO: Even better replacing states with db.startstates?

        objectivevalue = TuLiPa.getobjectivevalue(sp.prob)
        scenslopes = sp.scenslopes
        scenconstant = sp.scenconstant
    end
    return (scenix, objectivevalue, scenslopes, scenconstant)
end

# TODO: Simplify TuLiPa version of getscencutparameters?
function get_scencutparameters!(sp::ScenarioProblem, states::Dict{TuLiPa.StateVariableInfo, Float64})
    sp.scenconstant = TuLiPa.getobjectivevalue(sp.prob)

    for (i, (statevar, value)) in enumerate(states)
        (id, ix) = TuLiPa.getvarin(statevar)
        slope = TuLiPa.getfixvardual(sp.prob, id, ix)
        sp.scenconstant -= slope * value
        sp.scenslopes[i] = slope
    end

    return
end

function update_sp(t, scenix, subix)
    db = get_local_db()

    sp = db.sp[(scenix, subix)]
    scentime = get_scentphasein(t, get_scenarios(db.scenmod_stoch)[scenix], db.input)
    TuLiPa.update!(sp.prob, scentime)
    return
end

function update_statedependent_mp(stepnr, subsystem, prob, startstates, settings)
    has_statedependentprod(settings["problems"]["stochastic"]["master"]) && TuLiPa.statedependentprod!(prob, startstates, init=(stepnr==1))
    has_statedependentpump(settings["problems"]["stochastic"]["master"]) && TuLiPa.statedependentpump!(prob, startstates)
    return
end

update_prices_mp(stepnr, subix, subsystem::ExogenSubsystem) = nothing
function update_prices_mp(stepnr, subix, subsystem)
    db = get_local_db()
    mp = db.mp[subix]
    scenix = 1 # Which price to use for master problem?

    term_ppp = get_horizonterm_stoch(subsystem)
    for obj in TuLiPa.getobjects(mp.prob)
        update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    end
    return
end

update_prices_sp(stepnr, scenix, subix, subsystem::ExogenSubsystem) = nothing
function update_prices_sp(stepnr, scenix, subix, subsystem)
    db = get_local_db()
    sp = db.sp[(scenix, subix)]
              
    parentscenix = get_scenmod_stoch(db).scenarios[scenix].parentscenario
    term_ppp = get_horizonterm_stoch(subsystem)
    for obj in TuLiPa.getobjects(sp.prob)
        update_prices_obj(db, parentscenix, subix, stepnr, obj, term_ppp)
    end
    return
end

function update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    if obj isa TuLiPa.ExogenBalance
        periods = TuLiPa.getperiods(TuLiPa.gethorizon(obj))
        bid = TuLiPa.getid(obj)

        isupdated = isupdated_prices(db, scenix, term_ppp, bid, stepnr)
        !isupdated && update_local_price(db, scenix, term_ppp, bid, stepnr)

        updated, allvalues = db.prices_ppp[(scenix, term_ppp, bid)]
        obj.price.values .= allvalues[periods]
    end
end

function update_local_price(db, scenix, term_ppp, bid, stepnr)
    core_ppp = get_core_ppp(db, scenix)
    future = @spawnat core_ppp get_prices_from_core(scenix, term_ppp, bid)
    db.prices_ppp[(scenix, term_ppp, bid)] = (stepnr, fetch(future)) # TODO: Should we collect all prices or just relevant periods?
end

function get_core_ppp(db::LocalDB, scenix)
    for (_scenix, _core) in db.dist_ppp
        if _scenix == scenix
            return _core
        end
    end
end

function isupdated_prices(db, scenix, term_ppp, bid, stepnr)
    if haskey(db.prices_ppp, (scenix, term_ppp, bid))
        stored_stepnr, __ = db.prices_ppp[(scenix, term_ppp, bid)]
        return stored_stepnr == stepnr
    end
    return false
end

function find_obj_by_id(objects::Any, bid::TuLiPa.Id)
    for obj in objects
        if TuLiPa.getid(obj) == bid
            return obj
        end
    end
    return nothing
end

function get_prices_from_core(scenix, term_ppp, bid)
    db = get_local_db()

    ppp = get_ppp_term(db.ppp[scenix], term_ppp)

    obj = find_obj_by_id(ppp.objects, bid)
    horizon = TuLiPa.gethorizon(obj)

    return [-TuLiPa.getcondual(ppp, bid, t) for t in 1:TuLiPa.getnumperiods(horizon)]
end

function get_ppp_term(ppp, term::TermName)
    if term == LongTermName
        return ppp.longprob
    elseif term == MedTermName
        return ppp.medprob
    elseif term == ShortTermName
        return ppp.shortprob
    end
end

function perform_scenmod_sp(scenix, subix)
    db = get_local_db()
    sp = db.sp[(scenix, subix)]

    scenmod_stoch = get_scenmod_stoch(db)
    perform_scenmod!(scenmod_stoch, scenix, TuLiPa.getobjects(sp.prob))
    return
end

update_horizons_sp(scenix, subix, subsystem::ExogenSubsystem) = nothing
function update_horizons_sp(scenix, subix, subsystem)
    db = get_local_db()
    horizons = get_horizons(db)

    sp = db.sp[(scenix, subix)]

    parentscenix = get_scenmod_stoch(db).scenarios[scenix].parentscenario
    for commodity in collect(keys(sp.horizons))
        if commodity == "Battery"
            commodity_ppp = "Power"
        else
            commodity_ppp = commodity
        end
        if sp.horizons[commodity].subhorizon isa TuLiPa.ShortenedHorizon
            sp.horizons[commodity].subhorizon.subhorizon = horizons[parentscenix, get_horizonterm_stoch(subsystem), commodity_ppp]
        else
            @assert sp.horizons[commodity] isa TuLiPa.ShortenedHorizon
            sp.horizons[commodity].subhorizon = horizons[parentscenix, get_horizonterm_stoch(subsystem), commodity_ppp]
        end
    end
    return
end

function update_endconditions_sp(scenix, subix, t)
    db = get_local_db()

    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]

    endvaluemethod_sp = get_endvaluemethod_sp(subsystem)
    parentscenix = get_scenmod_stoch(db).scenarios[scenix].parentscenario

    storages = TuLiPa.getstorages(TuLiPa.getobjects(sp.prob))
    if endvaluemethod_sp == "monthly_price"
        exogenprice = find_firstprice(TuLiPa.getobjects(sp.prob))
        scentime = get_scentphasein(t, get_scenarios(db.scenmod_stoch)[scenix], db.input)
        scenprice = TuLiPa.getparamvalue(exogenprice, scentime + TuLiPa.getduration(TuLiPa.gethorizon(storages[1])), TuLiPa.MsTimeDelta(Week(4))) 

        for obj in storages
            enddual = scenprice * TuLiPa.getbalance(obj).metadata[TuLiPa.GLOBALENEQKEY]
            T = TuLiPa.getnumperiods(TuLiPa.gethorizon(TuLiPa.getbalance(obj)))
            TuLiPa.setobjcoeff!(sp.prob, TuLiPa.getid(obj), T, -enddual)
        end
    elseif endvaluemethod_sp == "startequalstop"
        set_endstates!(sp.prob, storages, db.startstates)
    elseif endvaluemethod_sp == "evp" # TODO: Store bid and period in sp (or subsystem?)
        core_evp = get_core_evp(db, parentscenix, subix)
        for obj in storages
            commodityname = TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(TuLiPa.getbalance(obj))))
            horizon_sp = TuLiPa.gethorizon(TuLiPa.getbalance(obj))
            duration_stoch = TuLiPa.getdurationtoend(horizon_sp)
            term_evp = get_horizonterm_evp(subsystem)
            horizon_evp = db.horizons[(parentscenix, term_evp, commodityname)]
            period_evp = TuLiPa.getendperiodfromduration(horizon_evp, duration_stoch)
            bid = TuLiPa.getid(TuLiPa.getbalance(obj))
            future = @spawnat core_evp get_balancedual_evp(parentscenix, subix, bid, period_evp)
            dual_evp = fetch(future)

            period_sp = TuLiPa.getnumperiods(horizon_sp)
            TuLiPa.setobjcoeff!(sp.prob, TuLiPa.getid(obj), period_sp, dual_evp)
        end
    elseif endvaluemethod_sp == "ppp"
        detailedrescopl = get_dataset(db)["detailedrescopl"]
        enekvglobaldict = get_dataset(db)["enekvglobaldict"]
        for obj in storages
            balance = TuLiPa.getbalance(obj)
            bid = TuLiPa.getid(balance)
            instancename = split(TuLiPa.getinstancename(bid), "Balance_")
            if haskey(detailedrescopl, instancename[2])
                balancename = detailedrescopl[instancename[2]]
                bid = TuLiPa.Id(bid.conceptname, instancename[1] * "Balance_" * balancename * "_hydro_reservoir") # TODO: This should be in the dataset
            end
            endperiod = TuLiPa.getlastperiod(TuLiPa.gethorizon(TuLiPa.getbalance(obj)))
            term_ppp = get_horizonterm_stoch(subsystem)
            core_ppp = get_core_ppp(db, parentscenix)
            future = @spawnat core_ppp get_balancedual_ppp(parentscenix, bid, endperiod, term_ppp)
            dual_ppp = fetch(future)
            if haskey(enekvglobaldict, instancename[2])
                dual_ppp *= enekvglobaldict[instancename[2]]
            end

            numperiods = TuLiPa.getnumperiods(TuLiPa.gethorizon(TuLiPa.getbalance(obj)))
            TuLiPa.setobjcoeff!(sp.prob, TuLiPa.getid(obj), numperiods, dual_ppp)
        end
    end

    return
end

function get_core_evp(db::LocalDB, scenix, subix)
    for (_scenix, _subix, _core) in db.dist_evp
        if (_scenix == scenix) && (_subix == subix)
            return _core
        end
    end
end

function get_core_sp(db::LocalDB, scenix, subix)
    for (_scenix, _subix, _core) in db.dist_sp
        if (_scenix == scenix) && (_subix == subix)
            return _core
        end
    end
end

function get_balancedual_evp(scenix, subix, bid, period)
    db = get_local_db()

    evp = db.evp[(scenix, subix)]
    dual_evp = TuLiPa.getcondual(evp.prob, bid, period)
    return dual_evp
end

function update_probabilities(cuts, scenmod)
    cuts.probabilities = [get_probability(scenario) for scenario in scenmod.scenarios]
end

# Util function under create_mp, create_sp -------------------------------------------------------------------------------------------------
function make_modelobjects_stochastic(db, scenix, subix, master)
    subsystem = get_subsystems(db)[subix]

    subelements, numperiods_powerhorizon, horizons = get_elements_with_horizons_stochastic(db, scenix, subsystem, master)

    if (get_numscen_sim(db.input) == get_numscen_stoch(db.input)) || master
        ifmscenix = scenix
    else
        ifmscenix = get_numscen_sim(db.input) + scenix
    end
    add_scenix_to_InflowParam(subelements, ifmscenix)

    aggzonecopl = get_aggzonecopl(get_aggzone(get_settings(db.input)))
    change_elements!(subelements, aggzonecopl=aggzonecopl)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = TuLiPa.getmodelobjects(subelements, validate=false)

    if !(subsystem isa ExogenSubsystem)
        if master && (get_horizonterm_stoch(subsystem) == ShortTermName) 
            # Removes spills from upper and lower storages in PHS, to avoid emptying reservoirs in master problem. TODO: Find better solution
            for id in keys(modelobjects)
                instance = TuLiPa.getinstancename(id)
                if occursin("Spill", instance) && occursin("_PHS_", instance)
                    pop!(modelobjects, id)
                end
            end
        end
    end

    return collect(values(modelobjects)), horizons
end

# Aggregate modelobjects and remove modelobjects not relevant for subsystems
function change_elements!(elements::Vector{TuLiPa.DataElement}; aggzonecopl::Dict=Dict()) # TODO: Replace with more user settings    
    delix = []
    powerbasebalances = []
    for (i,element) in enumerate(elements)
        
        # Power balance Arrows
        if element.conceptname == "Arrow"
            if element.value["Balance"] in keys(aggzonecopl)
                value = copy(element.value)
                value["Balance"] = aggzonecopl[element.value["Balance"]]
                elements[i] = TuLiPa.DataElement(element.conceptname, element.typename, element.instancename, value)
            end
        end

        # TODO: Get these from config
        if element.typename == "HydroRampingWithout"
            push!(delix,i)
        end
    end
    for deli in sort(delix; rev=true)
        popat!(elements, deli)
    end
    return elements
end

add_prices!(elements, subsystem::ExogenSubsystem, numperiods_powerhorizon, aggzonecopl) = nothing
function add_prices!(elements, subsystem, numperiods_powerhorizon, aggzonecopl)
    for pricearea in subsystem.priceareas
        if pricearea in keys(aggzonecopl)
            pricearea = aggzonecopl[element.value["Balance"]]
        end

        push!(elements, TuLiPa.getelement(TuLiPa.BALANCE_CONCEPT, "ExogenBalance", pricearea, 
        (TuLiPa.COMMODITY_CONCEPT, "Power"),
        (TuLiPa.PRICE_CONCEPT, "Price_" * pricearea)))
        push!(elements, TuLiPa.getelement(TuLiPa.PRICE_CONCEPT, "VectorPrice", "Price_" * pricearea,
        ("Vector", zeros(Float64, numperiods_powerhorizon))))
    end
    return 
end

function get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term, stochastic, ismaster)
    horizons = get_horizons(db)
    subelements = get_subelements(db, subsystem, ismaster)
    numscen_sim = get_numscen_sim(db.input)
    numscen_stoch = get_numscen_stoch(db.input)
    local numperiods_powerhorizon::Int

    probhorizons = Dict{CommodityName, TuLiPa.Horizon}()
    for commodity in get_commodities(subsystem)
        horizon = get_shortenedhorizon(horizons, scenix, term, commodity, startduration, endduration, stochastic, numscen_sim, numscen_stoch)
        probhorizons[commodity] = horizon
        set_horizon!(subelements, commodity, horizon)
        if commodity == "Power"
            numperiods_powerhorizon = TuLiPa.getnumperiods(horizon)
        end
    end

    return subelements, numperiods_powerhorizon, probhorizons
end

function get_subelements(db, ::ExogenSubsystem, ismaster::Bool)
    elements = get_elements(db.input)
    ismaster && return copy(elements)
    return copy_elements_iprogtype(elements, get_iprogtype(db.input), get_ifm_replacemap(db.input))
end

function get_subelements(db, subsystem::Union{EVPSubsystem, StochSubsystem}, ismaster::Bool)
    elements = get_elements(db.input)
    ismaster && return copy(elements)
    return copy_elements_iprogtype(elements[subsystem.dataelements], get_iprogtype(db.input), get_ifm_replacemap(db.input))
end

function get_elements_with_horizons_stochastic(db, scenix, subsystem::Union{ExogenSubsystem}, master)
    settings = get_settings(db.input)

    if master
        term = MasterTermName
        startduration = Millisecond(0)
        endduration = parse_duration(settings["horizons"]["master"], "termduration")
    else
        term = SubTermName
        startduration = parse_duration(settings["horizons"]["master"], "termduration")
        endduration = startduration + parse_duration(settings["horizons"]["sub"], "termduration")
    end
    subelements, numperiods_powerhorizon, horizons = get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term, true, master)

    return subelements, numperiods_powerhorizon, horizons
end

function get_elements_with_horizons_stochastic(db, scenix, subsystem::Union{EVPSubsystem, StochSubsystem}, master)
    settings = get_settings(db.input)

    term_ppp = get_horizonterm_stoch(subsystem)
    if master
        startduration = Millisecond(0)
        endduration = parse_duration(settings["horizons"]["clearing"], "termduration")
    else 
        startduration = parse_duration(settings["horizons"]["clearing"], "termduration")
        endduration = get_duration_stoch(subsystem)
    end
    subelements, numperiods_powerhorizon, horizons = get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term_ppp, true, master)

    return subelements, numperiods_powerhorizon, horizons
end

function get_shortenedhorizon(horizons::Dict{Tuple{ScenarioIx, TermName, CommodityName}, TuLiPa.Horizon}, scenix::ScenarioIx, term::TermName, commodity::CommodityName, startduration::Millisecond, endduration::Millisecond, stochastic::Bool, numscen_sim::Int, numscen_stoch::Int)
    if commodity == "Battery"
        subhorizon = horizons[(scenix, term, "Power")]
    else
        subhorizon = horizons[(scenix, term, commodity)]
    end
    if startduration.value == 0
        startperiod = 1
    else
        startperiod = TuLiPa.getendperiodfromduration(subhorizon, startduration) + 1
    end
    endperiod = TuLiPa.getendperiodfromduration(subhorizon, endduration)
    if !(subhorizon isa TuLiPa.ExternalHorizon) && !(term in [MasterTermName, SubTermName])
        subhorizon = TuLiPa.ExternalHorizon(subhorizon)
    end
    shortenedhorizon = TuLiPa.ShortenedHorizon(subhorizon, startperiod, endperiod)
    if stochastic && (numscen_sim != numscen_stoch)
        # TODO: Could replace this with functionality that ignores mustupdate and shrinkatleast if scenario has changed between steps
        shortenedhorizon = TuLiPa.IgnoreMustupdateMayshiftfromHorizon(shortenedhorizon)
    end
    return shortenedhorizon
end

function initialize_cuts!(modelobjects::Vector, cutobjects::Vector, maxcuts::Int, lb::Float64, numscen::Int)
    # Make a cutid
    cutname = TuLiPa.getinstancename(TuLiPa.getid(cutobjects[1]))
    cutid = TuLiPa.Id(TuLiPa.BOUNDARYCONDITION_CONCEPT,"StorageCuts" * cutname)
    
    # Make cut modelobject
    probabilities = [1/numscen for i in 1:numscen]
    cuts = TuLiPa.SimpleSingleCuts(cutid, cutobjects, probabilities, maxcuts, lb)
    push!(modelobjects, cuts)
    return cuts
end

function getcutobjects(modelobjects::Vector)
    cutobjects = Vector{Any}()
    for obj in modelobjects
        if TuLiPa.hasstatevariables(obj)
            if length(TuLiPa.getstatevariables(obj)) > 1
                error("Not supported") # TODO
            else
                push!(cutobjects,obj)
            end
        end
    end
    return cutobjects
end