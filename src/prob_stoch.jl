function create_mp(db::LocalDB, subix::SubsystemIx)
    scenix = 1 # TODO: Which scenario should be represented in the master problem? Not important due to phasein?
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = parse_duration(settings["horizons"]["clearing"], "termduration")

    # TODO: Not use ifm here as in clearnig?
    modelobjects = make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, true)

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
    div[MainTiming] = zeros(4)

    db.mp[subix] = MasterProblem(prob, cuts, states, div)

    return
end

function create_sp(db::LocalDB, scenix::ScenarioIx, subix::SubsystemIx)
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = parse_duration(settings["horizons"]["clearing"], "termduration")
    endduration = get_duration_stoch(subsystem)
    modelobjects = make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, false)
    states = get_states(modelobjects) # different order than mp.cuts.statevars, so only use length

    probmethod = parse_methods(settings["problems"]["stochastic"]["subs"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    div = Dict()
    div[MainTiming] = zeros(3)

    db.sp[(scenix, subix)] = ScenarioProblem(prob, zeros(length(states)), -1.0, div)

    return
end

function solve_stoch(t, stepnr, skipmed)
    db = get_local_db()

    for (subix, core) in db.dist_mp
        if core == db.core
            subsystem = db.subsystems[subix]
            if skipmed_check(subsystem, skipmed)
                mp = db.mp[subix]
                maintiming = mp.div[MainTiming]

                maintiming[4] = @elapsed begin
                    update_probabilities(mp.cuts, db.scenmod_stoch) # TODO: Add possibility for scenario modelling per subsystem
                    set_startstates!(mp.prob, TuLiPa.getstorages(TuLiPa.getobjects(mp.prob)), db.startstates)
                    update_prices_mp(stepnr, subix)
                    update_statedependent_mp(stepnr, subsystem, mp.prob, db.startstates, get_settings(db))
                end
                maintiming[1] = @elapsed TuLiPa.update!(mp.prob, t)

                @sync for _core in get_cores(db)
                    @spawnat _core update_sps(t, stepnr, subix)
                end
                    
                solve_benders(stepnr, subix)
                maintiming[3] = @elapsed final_solve_mp(t, mp.prob)
            end
        end
    end
end

function update_sps(t, stepnr, subix)
    db = get_local_db()

    for (_scenix, _subix, _core) in db.dist_sp
        if (_core == db.core) && (subix == _subix)
            maintiming = db.sp[(_scenix, subix)].div[MainTiming]
            maintiming[3] = @elapsed begin
                update_endconditions_sp(_scenix, subix)
                perform_scenmod_sp(_scenix, subix)
                update_prices_sp(stepnr, _scenix, subix)
            end
            maintiming[1] = @elapsed update_sp(t, _scenix, subix)
        end
    end
end

# Util functions for solve_mp ----------------------------------------------------------------------------------------------

function final_solve_mp(t::TuLiPa.ProbTime, prob)
    db = get_local_db()
    settings = get_settings(db)

    if get_headlosscost(settings["problems"]["stochastic"]["master"])
        TuLiPa.updateheadlosscosts!(TuLiPa.ReservoirCurveSlopeMethod(), prob, [prob], t)
        TuLiPa.solve!(prob)
        TuLiPa.resetheadlosscosts!(prob)
    end
end    

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

    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)
        count == 0 && TuLiPa.setwarmstart!(mp.prob, false)

        maintiming[2] += @elapsed begin
            if cutreuse # try to reuse cuts from last time step
                try
                    TuLiPa.solve!(mp.prob)
                catch
                    count == 0 && println("Retrying first iteration without cuts from last time step")
                    count > 0 && println("Restarting iterations without cuts from last time step")
                    TuLiPa.clearcuts!(mp.prob, mp.cuts)
                    TuLiPa.solve!(mp.prob)
                    cutreuse = false
                end
            else
                TuLiPa.solve!(mp.prob)
            end
        end

        maintiming[4] += @elapsed begin
            lb = TuLiPa.getvarvalue(mp.prob, TuLiPa.getfuturecostvarid(mp.cuts), 1)
            ub = 0.0

            count == 0 && TuLiPa.setwarmstart!(mp.prob, true)
            (count == 0 && cutreuse) && TuLiPa.clearcuts!(mp.prob, mp.cuts) # reuse cuts in first iteration
            
            TuLiPa.getoutgoingstates!(mp.prob, mp.states)
            cutix = TuLiPa.getcutix(mp.cuts) + 1
            if cutix > TuLiPa.getmaxcuts(mp.cuts)
                cutix = 1
            end
        end

        futures = []
        @sync for (_scenix, _subix, _core) in db.dist_sp
            if _subix == subix
                f = @spawnat _core solve_sp(_scenix, _subix, mp.states)
                push!(futures, f)
            end
        end

        maintiming[4] += @elapsed begin
            for future in futures
                scenix, objectivevalue, scenslopes, scenconstant = fetch(future)

                ub += objectivevalue*mp.cuts.probabilities[scenix]
                mp.cuts.scenslopes[scenix, cutix, :] .= scenslopes
                mp.cuts.scenconstants[scenix, cutix] = scenconstant
            end
        
            TuLiPa.updatecutparameters!(mp.prob, mp.cuts)
            TuLiPa.updatelastcut!(mp.prob, mp.cuts)
            count += 1
        end
    end
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
    init = false
    if stepnr == 1
        init = true
    end

    get_statedependentprod(settings["problems"]["stochastic"]["master"]) && TuLiPa.statedependentprod!(prob, startstates, init=init)
    get_statedependentpump(settings["problems"]["stochastic"]["master"]) && TuLiPa.statedependentpump!(prob, startstates)
    return
end

function update_prices_mp(stepnr, subix)
    db = get_local_db()
    subsystem = db.subsystems[subix]
    mp = db.mp[subix]
    scenix = 1 # Which price to use for master problem?

    term_ppp = get_horizonterm_stoch(subsystem)
    for obj in TuLiPa.getobjects(mp.prob)
        update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    end
    return
end

function update_prices_sp(stepnr, scenix, subix)
    db = get_local_db()
    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]
              
    term_ppp = get_horizonterm_stoch(subsystem)
    for obj in TuLiPa.getobjects(sp.prob)
        update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    end
    return
end

function update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    if obj isa TuLiPa.ExogenBalance
        periods = TuLiPa.getperiods(TuLiPa.gethorizon(obj))
        bid = TuLiPa.getid(obj)

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
    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]

    scenmod_stoch = get_scenmod_stoch(db)
    perform_scenmod!(scenmod_stoch, scenix, TuLiPa.getobjects(sp.prob))
    return
end

function update_endconditions_sp(scenix, subix)
    db = get_local_db()

    subsystem = db.subsystems[subix]
    sp = db.sp[(scenix, subix)]

    endvaluemethod_sp = get_endvaluemethod_sp(subsystem)

    storages = TuLiPa.getstorages(TuLiPa.getobjects(sp.prob))
    if endvaluemethod_sp == "monthly_price"
        exogenprice = findfirstprice(TuLiPa.getobjects(sp.prob))
        scentime = get_scentphasein(t, get_scenarios(db.scenmod_stoch)[scenix], db.input)
        scenprice = TuLiPa.getparamvalue(exogenprice, scentime + TuLiPa.getduration(TuLiPa.gethorizon(storages[1])), TuLiPa.MsTimeDelta(Week(4))) 

        for obj in storages
            enddual = scenprice * TuLiPa.getbalance(obj).metadata[TuLiPa.GLOBALENEQKEY]
            T = TuLiPa.getnumperiods(TuLiPa.gethorizon(TuLiPa.getbalance(obj)))
            TuLiPa.setobjcoeff!(sp.prob, TuLiPa.getid(obj), T, -enddual)
        end
    elseif endvaluemethod_sp == "startequalstop"
        TuLiPa.setendstates!(sp.prob, storages, db.startstates)
    elseif endvaluemethod_sp == "evp" # TODO: Store bid and period in sp (or subsystem?)
        core_evp = get_core_evp(db, scenix, subix)
        for obj in storages
            commodityname = TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(TuLiPa.getbalance(obj))))
            horizon_sp = TuLiPa.gethorizon(TuLiPa.getbalance(obj))
            duration_stoch = TuLiPa.getdurationtoend(horizon_sp)
            term_evp = get_horizonterm_evp(subsystem)
            horizon_evp = db.horizons[(scenix, term_evp, commodityname)]
            period_evp = TuLiPa.getendperiodfromduration(horizon_evp, duration_stoch)
            bid = TuLiPa.getid(TuLiPa.getbalance(obj))
            future = @spawnat core_evp get_balancedual_evp(scenix, subix, bid, period_evp)
            dual_evp = fetch(future)

            period_sp = TuLiPa.getnumperiods(horizon_sp)
            TuLiPa.setobjcoeff!(sp.prob, TuLiPa.getid(obj), period_sp, dual_evp)
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
    dual_evp = TuLiPa.getcondual(evp.prob, bid, period)
    return dual_evp
end

function update_probabilities(cuts, scenmod)
    cuts.probabilities = [get_probability(scenario) for scenario in scenmod.scenarios]
end

# Util function under create_mp, create_sp -------------------------------------------------------------------------------------------------
function make_modelobjects_stochastic(db, scenix, subix, startduration, endduration, master)
    subsystem = get_subsystems(db)[subix]
    term_ppp = get_horizonterm_stoch(subsystem)
    subelements, numperiods_powerhorizon = get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term_ppp, true)

    aggzonecopl = get_aggzonecopl(get_aggzone(get_settings(db.input)))
    change_elements!(subelements, aggzonecopl=aggzonecopl)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = TuLiPa.getmodelobjects(subelements, validate=false)

    if master && (get_horizonterm_stoch(subsystem) == ShortTermName) 
        # Removes spills from upper and lower storages in PHS, to avoid emptying reservoirs in master problem. TODO: Find better solution
        for id in keys(modelobjects)
            instance = TuLiPa.getinstancename(id)
            if occursin("Spill", instance) && occursin("_PHS_", instance)
                pop!(modelobjects, id)
            end
        end
    end

    return collect(values(modelobjects))
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

function get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term_ppp, stochastic)
    horizons = get_horizons(db)
    subelements = get_subelements(db, subsystem)
    numscen_sim = get_numscen_sim(db.input)
    numscen_stoch = get_numscen_stoch(db.input)
    local numperiods_powerhorizon::Int

    for commodity in get_commodities(subsystem)
        horizon = get_shortenedhorizon(horizons, scenix, term_ppp, commodity, startduration, endduration, stochastic, numscen_sim, numscen_stoch)
        set_horizon!(subelements, commodity, horizon)
        if commodity == "Power"
            numperiods_powerhorizon = TuLiPa.getnumperiods(horizon)
        end
    end

    add_scenix_to_InflowParam(subelements, scenix)

    return subelements, numperiods_powerhorizon
end

function get_subelements(db, subsystem::ExogenSubsystem)
    return copy_elements_iprogtype(get_elements(db.input), get_iprogtype(db.input), get_ifm_replacemap(db.input))
end

function get_subelements(db, subsystem::Union{EVPSubsystem, StochSubsystem})
    elements = get_elements(db.input)
    return copy_elements_iprogtype(elements[subsystem.dataelements], get_iprogtype(db.input), get_ifm_replacemap(db.input))
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


