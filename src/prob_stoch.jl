struct MasterProblem
    prob::Prob
    cuts::SimpleSingleCuts
    states::Dict{StateVariableInfo, Float64}
end

struct ScenarioProblem
    prob::Prob
end

function create_mp(db::LocalDB, subix::SubsystemIx)
    scenix = 1 # TODO: Which scenario should be represented in the master problem? Not important due to phasein?
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = Millisecond(Hour(get_settings(db)["time"]["steplength_hours"]))
    modelobjects = make_stochastic_modelobjects(db, scenix, subix, startduration, endduration, true)

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
    modelobjects = make_stochastic_modelobjects(db, scenix, subix, startduration, endduration, false)

    probmethod = parse_methods(settings["problems"]["stochastic"]["sub"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.sp[(scenix, subix)] = ScenarioProblem(prob)

    return
end

function solve_mp(T, t, delta, stepnr, thiscore)
    db = get_local_db()

    for (scenix, mp) in db.mp # Should db be input to all these functions?
        update_startstates_mp(stepnr, t)
        update_endstates_sp(stepnr, t)
        perform_scenmod_sp(stepnr)
        update_prices_mp(stepnr)
        update_prices_sp(stepnr)
        update_statedependent_mp(stepnr)
        update_mp(t)
        update_sp(t)
        solve_benders(stepnr)
        final_solve_mp()
    return
end

# Util functions for solve_mp ----------------------------------------------------------------------------------------------

function update_statedependent_mp(stepnr)
    db = get_local_db()
    settings = get_settings(db)

    init = false
    if stepnr == 1
        init = true
    end

    for (subix, mp) in db.mp
        getstatedependentprod(settings["problems"]["stochastic"]["master"]) && statedependentprod!(mp, db.startstates, init=init)
        getstatedependentpump(settings["problems"]["stochastic"]["master"]) && statedependentpump!(mp, db.startstates)
    end
    return
end

function update_prices_mp(stepnr)
    db = get_local_db()
    scenix = 1 # Which price to use for master problem?

    for (subix, mp) in db.mp
        for obj in getobjects(mp)
            update_prices_obj(db, scenix, subix, stepnr, obj)
        end
    end

    return
end

function update_prices_sp(stepnr)
    db = get_local_db()

    for ((scenix, subix), sp) in db.sp
        for obj in getobjects(sp)
            update_prices_obj(db, scenix, subix, stepnr, obj)
        end
    end

    return
end

function update_prices_obj(db, scenix, subix, stepnr, obj)
    if obj isa ExogenBalance
        term_ppp = get_term_ppp(db, subix, scenix)
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

get_ppp_term(ppp, ::LongTermName) = ppp.longprob
get_ppp_term(ppp, ::MedTermName) = ppp.medprob      
get_ppp_term(ppp, ::ShortTermName) = ppp.shortprob

function perform_scenmod_sp()
    db = get_local_db()

    scenmod_stoch = get_scenmod_stoch(db)
    for ((scenix, subix), sp) in db.sp
        perform_scenmod!(scenmod_stoch, scenix, getobjects(sp))
    end

    return
end

function update_endstates_sp(stepnr, t)
    db = get_local_db(scenmod_stoch)

    for ((scenix, subix), sp) in db.sp
        subsystem = get_subsystems(db)[subix]
        endvaluemethod_sp = get_endvaluemethod_sp(subsystem)

        storages = getstorages(sp.objects)
        if endvaluemethod_sp == "monthly_price"
            exogenprice = findfirstprice(sp.objects)
            scentime = get_scenariotime(t, get_scenarios(db.scenmod_stoch)[scenix], db.input, "normaltime")
            scenprice = getparamvalue(exogenprice, scentime + getduration(gethorizon(storages[1])), MsTimeDelta(Week(4))) 

            for obj in storages
                enddual = scenprice * getbalance(obj).metadata[GLOBALENEQKEY]
                T = getnumperiods(gethorizon(getbalance(obj)))
                setobjcoeff!(sp, getid(obj), T, -enddual)
            end
        elseif endvaluemethod_sp == "startequalstop"
            setendstates!(sp, storages, startstates)
        elseif endvaluemethod_sp == "evp"
            for obj in storages
                commodity = getcommodity(getbalance(obj))
                term_ppp = get_term_ppp(db, subix, scenix)
                horizon_evp = db.horizons[(scenix, term_ppp, commodity)]
                period = getendperiodfromduration(horizon_evp, get_duration_stoch(subsystem))

                core_evp = get_core_evp(db, scenix, subix)
                bid = getid(getbalance(storage))
                future = @spawnat core_evp get_balancedual_evp(scenix, subix, bid, period)
                dual_evp = fetch(future)

                setobjcoeff!(p, getid(obj), period, dual_evp)
            end
        end # TODO: Endvalue from ppp
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

    if stepnr == 1 # TODO: Might already be done by evp
        get_startstates_stoch_from_input(db, t)
    else # TODO: Copies all startstates
        if stepnr != db.stepnr_startstates
            get_startstates_from_cp(db)
            db.stepnr_startstates = stepnr
        end
    end
    for (subix, mp) in db.mp # TODO: set nonstorage startstates
        set_startstates!(mp.prob, get_storages(mp.prob), db.startstates)
    end
end

function get_startstates_stoch_from_input(db, t)
    settings = get_settings(db)
    dummystorages = getstorages(db.dummyobjects)
    get_startstates!(db.startstates, settings["problems"]["stochastic"], get_dataset(db), db.dummyobjects, dummystorages, t)
    startstates_max!(dummystorages, t, db.startstates)
    return
end

# Util function under create_mp, create_sp -------------------------------------------------------------------------------------------------
function make_stochastic_modelobjects(db, scenix, subix, startduration, endduration, master)
    subsystem = get_subsystems(db)[subix]
    subelements, numperiods_powerhorizon = get_elements_with_horizons(db, scenix, subix, startduration, endduration)

    add_prices!(subelements, subsystem, numperiods_powerhorizon)

    aggzone = get_aggzone(get_settings(db.input))
    change_elements!(subelements, aggzone)

    modelobjects = getmodelobjects(elements, validate=false)

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
function change_elements!(elements::Vector{DataElement}; aggzone::Dict=Dict()) # TODO: Replace with more user settings
    aggzonecopl = Dict()
    for (k,v) in aggzone
        for vv in v
            aggzonecopl["PowerBalance_" * vv] = "PowerBalance_" * k
        end
    end
    
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

add_prices!(elements, subsystem::ExogenSubsystem, numperiods_powerhorizon) = nothing
function add_prices!(elements, subsystem, numperiods_powerhorizon)
    priceareas = get_priceareas(subsystem)
    for area in priceareas
        push!(elements, getelement(BALANCE_CONCEPT, "ExogenBalance", "PowerBalance_" * area, 
        (COMMODITY_CONCEPT, "Power"),
        (PRICE_CONCEPT, "Price_" * area)))
        push!(elements, getelement(PRICE_CONCEPT, "VectorPrice", "Price_" * area,
        ("Vector", zeros(Float64, numperiods_powerhorizon))))
    end
    return 
end

function get_elements_with_horizons(db, scenix, subix, startduration, endduration)
    subsystem = get_subsystems(db)[subix]
    horizons = get_horizons(db)
    subelements = get_subelements(db, subsystem)
    term_ppp = get_term_ppp(db, subix, scenix)
    for commodity in get_commodities(subsystem)
        horizon = get_shortendhorizon_mp(horizons, scenix, term_ppp, commodity, startduration, endduration)
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
function get_term_ppp(db::LocalDB, subix::SubsystemIx, scenix::ScenarioIx)
    subsystem = get_subsystems(db)[subix]
    horizons = get_horizons(db)

    dummycommodity = get_commodities(subsystem)[1] # all of them have the same length
    duration_stoch = get_duration_stoch(subsystem)

    horizon_short = horizons[(scenix, ShortTermName, dummycommodity)]
    if duration_stoch < getduration(horizon_short) # TODO: also account for slack in case of reuse of watervalues
        return ShortTermName
    end
    horizon_med = horizons[(scenix, MedTermName, dummycommodity)]
    if duration_stoch < getduration(horizon_med) # TODO: also account for slack in case of reuse of watervalues
        return MedTermName
    end
    horizon_long = horizons[(scenix, LongTermName, dummycommodity)]
    @assert duration_stoch < getduration(horizon_long) # TODO: also account for slack in case of reuse of watervalues
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


