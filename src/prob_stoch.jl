function create_mp(db::LocalDB, subix::SubsystemIx)
    scenix = 1 # TODO: Which scenario should be represented in the master problem? Not important due to phasein?
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = Millisecond(Hour(get_settings(db)["time"]["steplength_hours"]))
    modelobjects = make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, true)

    maxcuts = settings["problems"]["stochastic"]["maxcuts"] # preallocate fixed number of cuts, no cut selection
    lb = settings["problems"]["stochastic"]["lb"] # lower bound of the future value in the first iteration
    numscen_stoch = get_numscen_stoch(db.input)
    cutobjects = getcutobjects(masterobjects)
    cuts = initialize_cuts!(masterobjects, cutobjects, maxcuts, lb, numscen_stoch)
    states = getstates(cutobjects) # state variables in master and subs for boundary reservoirs

    probmethod = parse_methods(settings["problems"]["stochastic"]["master"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.mp[subix] = MasterProblem(prob, cuts, states)

    return
end

function create_sp(db::LocalDB, scenix::ScenarioIx, subix::SubsystemIx)
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(Hour(get_settings(db)["time"]["steplength_hours"]))
    endduration = get_duration_stoch(subsystem)
    modelobjects = make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, false)

    probmethod = parse_methods(settings["problems"]["stochastic"]["sub"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.sp[(scenix, subix)] = ScenarioProblem(prob)

    return
end

function solve_mp(t, delta, stepnr)
    db = get_local_db()

    update_startstates_mp(stepnr, t)
    update_endstates_sp(stepnr, t)
    perform_scenmod_sp()
    update_prices_mp(stepnr)
    update_prices_sp(stepnr)
    update_statedependent_mp(stepnr)
    update_mp(t)
    update_sp(t)
    solve_benders(stepnr)
    final_solve_mp(t)
end

# Util functions for solve_mp ----------------------------------------------------------------------------------------------

function final_solve_mp(t::ProbTime)
    db = get_local_db()
    settings = get_settings(db)

    if get_headlosscost(settings["problems"]["stochastic"]["master"])
        for (subix, core) in db.dist_mp
            if core == db.core
                mp = db.mp[subix]

                updateheadlosscosts!(ReservoirCurveSlopeMethod(), mp.prob, [mp.prob], t)
                solve!(mp.prob)
                resetheadlosscosts!(mp.prob)
            end
        end
    end
end

function solve_benders(stepnr)
    db = get_local_db()
    settings = get_settings(db)

    for (subix, core) in db.dist_mp
        if core == db.core
            mp = db.mp[subix]

            count = 0
            cutreuse = false
            ub = 0
            lb = mp.cuts.lower_bound
            reltol = settings["problems"]["stochastic"]["reltol"] # relative tolerance

            while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)

                count == 0 && setwarmstart!(mp.prob, false)
        
                if cutreuse # try to reuse cuts from last time step
                    try
                        solve!(mp.prob)
                    catch
                        count == 0 && println("Retrying first iteration without cuts from last time step")
                        count > 0 && println("Restarting iterations without cuts from last time step")
                        clearcuts!(mp.prob, mp.cuts)
                        solve!(mp.prob)
                        cutreuse = false
                    end
                else
                    solve!(mp.prob)
                end
        
                lb = getvarvalue(mp.prob, getfuturecostvarid(mp.cuts), 1)
                ub = 0
        
                count == 0 && setwarmstart!(mp.prob, true)
                (count == 0 && cutreuse) && clearcuts!(mp.prob, mp.cuts) # reuse cuts in first iteration
                
                getoutgoingstates!(mp.prob, mp.states)
                cutix = oldcutix + 1
                if cutix > maxcuts
                    cutix = 1
                end

                @sync for (_scenix, _subix, _core) in db.dist_sp
                    if _subix == subix
                        @spawnat _core solve_sp(_scenix, _subix, mp.states)
                    end
                end

                for (_scenix, _subix, _core) in db.dist_sp
                    if _subix == subix
                        future = @spawnat _core get_data_sp(_scenix, _subix)
                        objectivevalue, scenslopes, scenconstants = fetch(future)

                        ub += objectivevalue*mp.cuts.probabilities[_scenix]
                        mp.cuts.scenslopes[_scenix, cutix, :] .= scenscenslopes
                        mp.cuts.scenconstants[_scenix, cutix] = scenconstant
                    end
                end
        
                updatecutparameters!(mp.prob, mp.cuts)
                if (count == 0 && cutreuse) 
                    updatecuts!(mp.prob, mp.cuts)
                else
                    updatelastcut!(mp.prob, mp.cuts)
                end
                count += 1
                # display(ub)
                # display(abs((lb-ub)/lb))
                # display(abs(ub-lb))
                # display(cuts.slopes)
            end
        end
    end
    return
end

function get_data_sp(scenix, subix)
    db = get_local_db()

    sp = db.sp[(scenix, subix)]

    objectivevalue = getobjectivevalue(sp.prob)
    scenslopes = sp.scenslopes
    scenconstant = sp.scenconstant

    return (objectivevalue, scenslopes, scenconstant)
end

function solve_sp(scenix, subix, states)
    db = get_local_db()

    sp = db_sp[(scenix, subix)]
    setingoingstates!(sp.prob, states)  
    solve!(sp.prob)
    get_scencutparameters!(sp, states)
end

# TODO: Simplify TuLiPa version of getscencutparameters?
function get_scencutparameters!(sp::ScenarioProblem, states::Dict{StateVariableInfo, Float64})
    sp.scenconstant = getobjectivevalue(sp.prob)
    
    for (i, (statevar, value)) in enumerate(states)
        (id, ix) = getvarin(statevar)
        slope = getfixvardual(sp.prob, id, ix)
        sp.scenconstant -= slope * value
        sp.scenslopes[i] = slope
    end

    return
end

function update_mp(t)
    db = get_local_db()

    for (subix, core) in db.dist_mp
        if core == db.core
            mp = db.mp[subix]
            update!(mp.prob, t)
        end
    end
    return
end

function update_sp(t)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_sp
        if core == db.core
            sp = db_sp[(scenix, subix)]
            scentime = get_scentphasein(t, get_scenarios(db.scenmod_stoch)[scenix], db.input)
            update!(sp.prob, scentime)
        end
    end
    return
end

function update_statedependent_mp(stepnr)
    db = get_local_db()
    settings = get_settings(db)

    init = false
    if stepnr == 1
        init = true
    end

    for (subix, core) in db.dist_mp
        if core == db.core
            mp = db.mp[subix]
            get_statedependentprod(settings["problems"]["stochastic"]["master"]) && statedependentprod!(mp.prob, db.startstates, init=init)
            get_statedependentpump(settings["problems"]["stochastic"]["master"]) && statedependentpump!(mp.prob, db.startstates)
        end
    end
    return
end

function update_prices_mp(stepnr)
    db = get_local_db()
    scenix = 1 # Which price to use for master problem?

    for (subix, core) in db.dist_mp
        if core == db.core
            mp = db.mp[subix]
            subsystem = db.subsystems[subix]
            duration_stoch = get_duration_stoch(subsystem)
            for obj in getobjects(mp.prob)
                update_prices_obj(db, scenix, subix, stepnr, obj, duration_stoch)
            end
        end
    end

    return
end

function update_prices_sp(stepnr)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_sp
        if core == db.core
            sp = db_sp[(scenix, subix)]
            subsystem = db.subsystems[subix]
            duration_stoch = get_duration_stoch(subsystem)
            for obj in getobjects(sp.prob)
                update_prices_obj(db, scenix, subix, stepnr, obj, duration_stoch)
            end
        end
    end

    return
end

function update_prices_obj(db, scenix, subix, stepnr, obj, duration)
    if obj isa ExogenBalance
        term_ppp = get_term_ppp(db, subix, scenix, duration)
        periods = get_periods(gethorizon(obj)) # TODO: Implement get_periods
        bid = getid(obj)

        isupdated = isupdated_prices(db, scenix, term_ppp, bid, stepnr)
        !isupdated && update_local_price(db, scenix, term_ppp, bid)

        updated, allvalues = db.prices_ppp[(scenix, term_ppp, bid)]
        obj.price.values .= allvalues[periods]
    end
end

function update_local_price(db, scenix, term_ppp, bid)
    core_ppp = get_core_ppp(db, scenix)
    future = @spawnat core_ppp get_prices_from_core(scenix, term_ppp, bid)
    db.prices_ppp[(scenix, term_ppp, bid)] = fetch(future) # TODO: Should we collect all prices or just relevant periods?
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
        updated, allvalues = db.prices_ppp[(scenix, term_ppp, bid)]
        if updated != stepnr
            return false
        else
            return true
        end
    else
        return false
    end
end

function find_obj_by_id(objects::Any, bid::Id)
    for obj in objects
        if getid(obj) == bid
            return obj
        end
    end
    return nothing
end

function get_prices_from_core(scenix, term_ppp, bid)
    db = get_local_db()

    ppp = get_ppp_term(db.ppp[scenix], term_ppp)

    obj = find_obj_by_id(ppp, bid)
    horizon = gethorizon(obj)

    return [getcondual(ppp, bid, t) for t in 1:getnumperiods(horizon)]
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

function perform_scenmod_sp()
    db = get_local_db()

    scenmod_stoch = get_scenmod_stoch(db)
    for (scenix, subix, core) in db.dist_sp
        if core == db.core
            sp = db_sp[(scenix, subix)]
            perform_scenmod!(scenmod_stoch, scenix, getobjects(sp))
        end
    end

    return
end

function update_endstates_sp(stepnr, t)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_sp
        if core == db.core
            sp = db_sp[(scenix, subix)]
            subsystem = get_subsystems(db)[subix]
            endvaluemethod_sp = get_endvaluemethod_sp(subsystem)

            storages = getstorages(getobjects(sp.prob))
            if endvaluemethod_sp == "monthly_price"
                exogenprice = findfirstprice(getobjects(sp.prob))
                scentime = get_scentphasein(t, get_scenarios(db.scenmod_stoch)[scenix], db.input)
                scenprice = getparamvalue(exogenprice, scentime + getduration(gethorizon(storages[1])), MsTimeDelta(Week(4))) 

                for obj in storages
                    enddual = scenprice * getbalance(obj).metadata[GLOBALENEQKEY]
                    T = getnumperiods(gethorizon(getbalance(obj)))
                    setobjcoeff!(sp.prob, getid(obj), T, -enddual)
                end
            elseif endvaluemethod_sp == "startequalstop"
                setendstates!(sp.prob, storages, startstates)
            elseif endvaluemethod_sp == "evp"
                for obj in storages
                    commodity = getcommodity(getbalance(obj))
                    duration_stoch = get_duration_stoch(subsystem)
                    term_ppp = get_term_ppp(db, subix, scenix, duration_stoch)
                    horizon_evp = db.horizons[(scenix, term_ppp, commodity)]
                    period = getendperiodfromduration(horizon_evp, duration_stoch)

                    core_evp = get_core_evp(db, scenix, subix)
                    bid = getid(getbalance(storage))
                    future = @spawnat core_evp get_balancedual_evp(scenix, subix, bid, period)
                    dual_evp = fetch(future)

                    setobjcoeff!(p, getid(obj), period, dual_evp)
                end
            end # TODO: Endvalue from ppp
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

function get_balancedual_evp(scenix, subix, bid, period)
    db = get_local_db()

    evp = db.evp[(scenix, subix)]
    dual_evp = -getcondual(evp.prob, bid, period)
    return dual_evp
end

function update_startstates_mp(stepnr, t)
    db = get_local_db()

    # TODO: Check if any of the scenarios are on this core first
    if stepnr == 1 # TODO: Might already be done by evp
        get_startstates_stoch_from_input(db, t)
    else # TODO: Copies all startstates
        if stepnr != db.stepnr_startstates
            get_startstates_from_cp(db)
            db.stepnr_startstates = stepnr
        end
    end 
        
    # TODO: set nonstorage startstates
    for (subix, core) in db.dist_mp
        if core == db.core
            mp = db.mp[subix]
            set_startstates!(mp.prob, get_storages(mp.prob), db.startstates)
        end
    end
end

function get_startstates_stoch_from_input(db, t)
    settings = get_settings(db)
    dummystorages = getstorages(first(db.dummyobjects))
    get_startstates!(db.startstates, settings["problems"]["stochastic"], get_dataset(db), first(db.dummyobjects), dummystorages, t)
    startstates_max!(dummystorages, t, db.startstates)
    return
end

# Util function under create_mp, create_sp -------------------------------------------------------------------------------------------------
function make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, master)
    subsystem = get_subsystems(db)[subix]
    subelements, numperiods_powerhorizon = get_elements_with_horizons(db, scenix, subix, startduration, endduration)

    aggzonecopl = get_aggzonecopl(get_settings(db.input))
    change_elements!(subelements, aggzonecopl)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = getmodelobjects(subelements, validate=false)

    if master == true
        # Removes spills from upper and lower storages in PHS, to avoid emptying reservoirs in master problem
        for id in keys(modelobjects)
            instance = getinstancename(id)
            if occursin("Spill", instance) && occursin("_PHS_", instance)
                delete!(modelobjects, id)
            end
        end
    end

    return modelobjects
end

# Aggregate modelobjects and remove modelobjects not relevant for subsystems
function change_elements!(elements::Vector{DataElement}; aggzonecopl::Dict=Dict()) # TODO: Replace with more user settings    
    delix = []
    powerbasebalances = []
    for (i,element) in enumerate(elements)
        
        # Power balance Arrows
        if element.conceptname == "Arrow"
            if element.value["Balance"] in keys(aggzonecopl)
                value = copy(element.value)
                value["Balance"] = aggzonecopl[element.value["Balance"]]
                elements[i] = DataElement(element.conceptname, element.typename, element.instancename, value)
            end
        end

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
    priceareas_added = []

    for element in elements
        if element.conceptname == "Arrow"
            if element.value["Balance"] in keys(aggzonecopl)
                pricearea = element.value["Balance"]
            else
                pricearea = aggzonecopl[element.value["Balance"]]
            end
            if !(pricearea in priceareas_added)
                push!(elements, getelement(BALANCE_CONCEPT, "ExogenBalance", pricearea, 
                (COMMODITY_CONCEPT, "Power"),
                (PRICE_CONCEPT, "Price_" * pricearea)))
                push!(elements, getelement(PRICE_CONCEPT, "VectorPrice", "Price_" * pricearea,
                ("Vector", zeros(Float64, numperiods_powerhorizon))))

                push!(priceareas_added, pricearea)
            end
            
        end
    end

    return 
end

function get_elements_with_horizons(db, scenix, subix, startduration, endduration)
    subsystem = get_subsystems(db)[subix]
    horizons = get_horizons(db)
    subelements = get_subelements(db, subsystem)
    duration_stoch = get_duration_stoch(subsystem)
    term_ppp = get_term_ppp(db, subix, scenix, duration_stoch)
    for commodity in get_commodities(subsystem)
        horizon = get_shortenedhorizon(horizons, scenix, term_ppp, commodity, startduration, endduration)
        set_horizon!(subelements, commodity, horizon)
        if commodity == "Power"
            numperiods_powerhorizon = getnumperiods(horizon)
        end
    end
    return subelements, numperiods_powerhorizon
end

get_subelements(db, subsystem::ExogenSubsystem) = copy(get_elements(db))
function get_subelements(db, subsystem::Union{EVPSubsystem, StochSubsystem})
    elements = get_elements(db)
    return copy(elements[subsystem.dataelements])
end

# Which time resolution (short, med, long) should we use horizons and prices from
# TODO: Should we use different terms for master and subproblems?
function get_term_ppp(db::LocalDB, subix::SubsystemIx, scenix::ScenarioIx, duration::Millisecond)
    subsystem = get_subsystems(db)[subix]
    horizons = get_horizons(db)

    dummycommodity = get_commodities(subsystem)[1] # all of them have the same length

    horizon_short = horizons[(scenix, ShortTermName, dummycommodity)]
    if duration < getduration(horizon_short) # TODO: also account for slack in case of reuse of watervalues
        return ShortTermName
    end
    horizon_med = horizons[(scenix, MedTermName, dummycommodity)]
    if duration < getduration(horizon_med) # TODO: also account for slack in case of reuse of watervalues
        return MedTermName
    end
    horizon_long = horizons[(scenix, LongTermName, dummycommodity)]
    @assert duration < getduration(horizon_long) # TODO: also account for slack in case of reuse of watervalues
    return LongTermName   
end

function get_shortenedhorizon(horizons::Dict{Tuple{ScenarioIx, TermName, CommodityName}, Horizon}, scenix::ScenarioIx, term::TermName, commodity::CommodityName, startduration::Millisecond, endduration::Millisecond)
    subhorizon = horizons[(scenix, term, commodity)]
    if startduration.value == 0
        startperiod = 1
    else
        startperiod = getendperiodfromduration(subhorizon, startduration) + 1
    end
    endperiod = getendperiodfromduration(subhorizon, endduration)

    return ShortenedHorizon(subhorizon, startperiod, endperiod)
end

function initialize_cuts!(modelobjects::Vector, cutobjects::Vector, maxcuts::Int, lb::Float64, numscen::Int)
    # Make a cutid
    cutname = getinstancename(getid(modelobjects[1]))
    cutid = Id(BOUNDARYCONDITION_CONCEPT,"StorageCuts" * cutname)
    
    # Make cut modelobject
    probabilities = [1/numscen for i in numscen]
    cuts = SimpleSingleCuts(cutid, cutobjects, probabilities, maxcuts, lb)
    push!(modelobjects, cuts)
    return cuts
end

function getcutobjects(modelobjects::Vector)
    cutobjects = Vector{Any}()
    for obj in modelobjects
        if hasstatevariables(obj)
            if length(getstatevariables(obj)) > 1
                error("Not supported") # TODO
            else
                push!(cutobjects,obj)
            end
        end
    end
    return cutobjects
end


