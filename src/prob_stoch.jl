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

# Util function under create_mp, create_sp -------------------------------------------------------------------------------------------------
function make_stochasticmodelobjects(db, scenix, subix, startduration, endduration, master)
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
        horizon = get_shortendhorizon_mp(horizons, scenix, term_ppp, commodity, duration_mp)
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

    horizon_short = horizons[(scenix, "Short", dummycommodity)]
    if duration_stoch < getduration(horizon_short) # TODO: also account for slack in case of reuse of watervalues
        return "Short"
    end
    horizon_med = horizons[(scenix, "Med", dummycommodity)]
    if duration_stoch < getduration(horizon_med) # TODO: also account for slack in case of reuse of watervalues
        return "Med"
    end
    horizon_long = horizons[(scenix, "Long", dummycommodity)]
    @assert duration_stoch < getduration(horizon_long) # TODO: also account for slack in case of reuse of watervalues
    return "Long"    
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

# Initialize cuts
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

# Get cutobjects
function getcutobjects(modelobjects::Vector)
    cutobjects = Vector{Any}()
    for obj in modelobjects
        if hasstatevariables(obj)
            if length(getstatevariables(obj)) > 1
                error("Not supported")
            else
                push!(cutobjects,obj)
            end
        end
    end
    return cutobjects
end


