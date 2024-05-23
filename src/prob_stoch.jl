function create_mp(db::LocalDB, subix::SubsystemIx)
    scenix = 1 # TODO: Which scenario should be represented in the master problem? Not important due to phasein?
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = parse_duration(settings["horizons"]["clearing"], "termduration")
    modelobjects = make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, true)

    maxcuts = settings["problems"]["stochastic"]["maxcuts"] # preallocate fixed number of cuts, no cut selection
    lb = settings["problems"]["stochastic"]["lb"] # lower bound of the future value in the first iteration
    numscen_stoch = get_numscen_stoch(db.input)
    cutobjects = getcutobjects(modelobjects)
    cuts = initialize_cuts!(modelobjects, cutobjects, maxcuts, lb, numscen_stoch)
    states = get_states(cutobjects) # state variables in master and subs for boundary reservoirs

    probmethod = parse_methods(settings["problems"]["stochastic"]["master"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.mp[subix] = MasterProblem(prob, cuts, states, Dict())

    return
end

function create_sp(db::LocalDB, scenix::ScenarioIx, subix::SubsystemIx)
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = parse_duration(settings["horizons"]["clearing"], "termduration")
    endduration = get_duration_stoch(subsystem)
    modelobjects = make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, false)
    cutobjects = getcutobjects(modelobjects)
    states = get_states(cutobjects)

    probmethod = parse_methods(settings["problems"]["stochastic"]["subs"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.sp[(scenix, subix)] = ScenarioProblem(prob, zeros(length(states)), -1.0, Dict())

    return
end

function solve_stoch(t, delta, stepnr, skipmed)
    db = get_local_db()

    for (subix, core) in db.dist_mp
        if core == db.core
            subsystem = db.subsystems[subix]
            if skipmed_check(subsystem, skipmed)
                mp = db.mp[subix]

                update_probabilities(mp.cuts, db.scenmod_stoch) # TODO: Add possibility for scenario modelling per subsystem
                set_startstates!(mp.prob, getstorages(getobjects(mp.prob)), db.startstates)
                update_prices_mp(stepnr, subix)
                update_statedependent_mp(stepnr, subsystem, mp.prob, db.startstates, get_settings(db))
                update!(mp.prob, t)

                @sync for _core in get_cores(db)
                    @spawnat _core update_sps(t, stepnr, subix)
                end
                    
                println("solve_benders")
                solve_benders(stepnr, subix)
                println("final_solve_mp")
                final_solve_mp(t, mp.prob)
            end
        end
    end
end

function update_sps(t, stepnr, subix)
    db = get_local_db()

    for (_scenix, _subix, _core) in db.dist_sp
        if (_core == db.core) && (subix == _subix)
            update_endstates_sp(_scenix, subix)
            perform_scenmod_sp(_scenix, subix)
            update_prices_sp(stepnr, _scenix, subix)
            update_sp(t, _scenix, subix)
        end
    end
end

# Util functions for solve_mp ----------------------------------------------------------------------------------------------

function final_solve_mp(t::ProbTime, prob)
    db = get_local_db()
    settings = get_settings(db)

    if get_headlosscost(settings["problems"]["stochastic"]["master"])
        updateheadlosscosts!(ReservoirCurveSlopeMethod(), prob, [prob], t)
        solve!(prob)
        resetheadlosscosts!(prob)
    end
end    

function solve_benders(stepnr, subix)
    db = get_local_db()
    settings = get_settings(db)

    mp = db.mp[subix]

    count = 0
    cutreuse = true
    if stepnr == 1
        cutreuse = false
    end
    ub = 0.0
    lb = mp.cuts.lower_bound
    reltol = settings["problems"]["stochastic"]["reltol"] # relative tolerance
    println(getid(getobjects(mp.prob)[1]))

    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)
        println(count)
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
        ub = 0.0

        count == 0 && setwarmstart!(mp.prob, true)
        (count == 0 && cutreuse) && clearcuts!(mp.prob, mp.cuts) # reuse cuts in first iteration
        
        getoutgoingstates!(mp.prob, mp.states)
        cutix = getcutix(mp.cuts) + 1
        if cutix > getmaxcuts(mp.cuts)
            cutix = 1
        end
        println(collect(values(mp.states)))

        @sync for (_scenix, _subix, _core) in db.dist_sp
            if _subix == subix
                @spawnat _core solve_sp(_scenix, _subix, mp.states)
            end
        end

        for (_scenix, _subix, _core) in db.dist_sp
            if _subix == subix
                future = @spawnat _core get_data_sp(_scenix, _subix)
                objectivevalue, scenslopes, scenconstant = fetch(future)

                ub += objectivevalue*mp.cuts.probabilities[_scenix]
                mp.cuts.scenslopes[_scenix, cutix, :] .= scenslopes
                mp.cuts.scenconstants[_scenix, cutix] = scenconstant

                # println(objectivevalue)
                # println(scenslopes)
                # println(scenconstant)
            end
        end
        println(ub)
        println(lb)

        updatecutparameters!(mp.prob, mp.cuts)
        if (count == 0 && cutreuse) 
            updatecuts!(mp.prob, mp.cuts)
        else
            updatelastcut!(mp.prob, mp.cuts)
        end
        count += 1

        if count == 20
            for obj in getobjects(mp.prob)
                println(getid(obj))
            end
            error("Not converging")
        end
        # display(ub)
        # display(abs((lb-ub)/lb))
        # display(abs(ub-lb))
        # display(cuts.slopes)
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

    sp = db.sp[(scenix, subix)]
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

function update_sp(t, scenix, subix)
    db = get_local_db()

    sp = db.sp[(scenix, subix)]
    scentime = get_scentphasein(t, get_scenarios(db.scenmod_stoch)[scenix], db.input)
    update!(sp.prob, scentime)
    return
end

function update_statedependent_mp(stepnr, subsystem, prob, startstates, settings)
    init = false
    if stepnr == 1
        init = true
    end

    get_statedependentprod(settings["problems"]["stochastic"]["master"]) && statedependentprod!(prob, startstates, init=init)
    get_statedependentpump(settings["problems"]["stochastic"]["master"]) && statedependentpump!(prob, startstates)
    return
end

function update_prices_mp(stepnr, subix)
    db = get_local_db()
    subsystem = db.subsystems[subix]
    mp = db.mp[subix]
    scenix = 1 # Which price to use for master problem?

    term_ppp = get_horizonterm_stoch(subsystem)
    for obj in getobjects(mp.prob)
        update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    end
    return
end

function update_prices_sp(stepnr, scenix, subix)
    db = get_local_db()
    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]
              
    term_ppp = get_horizonterm_stoch(subsystem)
    for obj in getobjects(sp.prob)
        update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    end
    return
end

function update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    if obj isa ExogenBalance
        periods = getperiods(gethorizon(obj))
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
    db.prices_ppp[(scenix, term_ppp, bid)] = (true, fetch(future)) # TODO: Should we collect all prices or just relevant periods?
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

    obj = find_obj_by_id(ppp.objects, bid)
    horizon = gethorizon(obj)

    return [-getcondual(ppp, bid, t) for t in 1:getnumperiods(horizon)]
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
    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]

    scenmod_stoch = get_scenmod_stoch(db)
    perform_scenmod!(scenmod_stoch, scenix, getobjects(sp.prob))
    return
end

function update_endstates_sp(scenix, subix)
    db = get_local_db()

    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]

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
            commodityname = getinstancename(getid(getcommodity(getbalance(obj))))
            horizon_sp = gethorizon(getbalance(obj))
            duration_stoch = getdurationtoend(horizon_sp)
            term_evp = get_horizonterm_evp(subsystem)
            horizon_evp = db.horizons[(scenix, term_evp, commodityname)]
            period_evp = getendperiodfromduration(horizon_evp, duration_stoch)
            core_evp = get_core_evp(db, scenix, subix)
            bid = getid(getbalance(obj))
            future = @spawnat core_evp get_balancedual_evp(scenix, subix, bid, period_evp)
            dual_evp = fetch(future)
            # println("dual evp $(dual_evp)")

            period_sp = getnumperiods(horizon_sp)
            setobjcoeff!(sp.prob, getid(obj), period_sp, dual_evp)
        end
    end # TODO: Endvalue from ppp

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
    # for p in 1:period
    #     println("$(getinstancename(bid)) $(p) $(getcondual(evp.prob, bid, period))")
    # end
    dual_evp = getcondual(evp.prob, bid, period)
    return dual_evp
end

function update_probabilities(cuts, scenmod)
    cuts.probabilities -= [get_probability(scenario) for scenario in scenmod.scenarios]
end

# Util function under create_mp, create_sp -------------------------------------------------------------------------------------------------
function make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, master)
    subsystem = get_subsystems(db)[subix]
    term_ppp = get_horizonterm_stoch(subsystem)
    subelements, numperiods_powerhorizon = get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term_ppp)

    aggzonecopl = get_aggzonecopl(get_aggzone(get_settings(db.input)))
    change_elements!(subelements, aggzonecopl=aggzonecopl)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = getmodelobjects(subelements, validate=false)

    if master
        # Removes spills from upper and lower storages in PHS, to avoid emptying reservoirs in master problem
        for id in keys(modelobjects)
            instance = getinstancename(id)
            if occursin("Spill", instance) && occursin("_PHS_", instance)
                delete!(modelobjects, id)
            end
        end
    end

    return collect(values(modelobjects))
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

        push!(elements, getelement(BALANCE_CONCEPT, "ExogenBalance", pricearea, 
        (COMMODITY_CONCEPT, "Power"),
        (PRICE_CONCEPT, "Price_" * pricearea)))
        push!(elements, getelement(PRICE_CONCEPT, "VectorPrice", "Price_" * pricearea,
        ("Vector", zeros(Float64, numperiods_powerhorizon))))
    end
    return 
end

function get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term_ppp)
    horizons = get_horizons(db)
    subelements = get_subelements(db, subsystem)
    local numperiods_powerhorizon::Int

    for commodity in get_commodities(subsystem)
        horizon = get_shortenedhorizon(horizons, scenix, term_ppp, commodity, startduration, endduration)
        set_horizon!(subelements, commodity, horizon)
        if commodity == "Power"
            numperiods_powerhorizon = getnumperiods(horizon)
        end
    end

    # Needed for inflow_models.includeModeledInflow! to work
    add_scenix_to_ModeledInflow_elements(subelements, scenix)

    return subelements, numperiods_powerhorizon
end

get_subelements(db, subsystem::ExogenSubsystem) = copy(get_elements(db.input))
function get_subelements(db, subsystem::Union{EVPSubsystem, StochSubsystem})
    elements = get_elements(db.input)
    return copy(elements[subsystem.dataelements])
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
    probabilities = [1/numscen for i in 1:numscen]
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


