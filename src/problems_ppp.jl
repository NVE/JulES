struct PricePrognosisProblem
    longprob::Prob
    medprob::Prob
    shortprob::Prob
    long_prices::
    med_prices::
    short_prices::
    short_nonstoragestates::Dict{StateVariableInfo, Float64}
end

function create_ppp(db::LocalDB, scenarioix::Int)
    lhh = db.horizons[(scenarioix, "long", "Hydro")]
    lph = db.horizons[(scenarioix, "long", "Power")]
    mhh = db.horizons[(scenarioix, "med", "Hydro")]
    mph = db.horizons[(scenarioix, "med", "Power")]
    shh = db.horizons[(scenarioix, "short", "Hydro")]
    sph = db.horizons[(scenarioix, "short", "Power")]

    aggzone = getaggzone(db.settings)
    aggsupplyn = db.settings["problems"]["prognosis"]["aggsupplyn"]
    removestoragehours = db.settings["problems"]["shorttermstoragecutoff_hours"]
    residualarealist = db.settings["problems"]["prognosis"]["residualarealist"]

    # Create long problems
    longobjects = make_obj(db.input.progelements, lhh, lph)
    simplify!(longobjects; aggzone=aggzone, aggsupplyn=aggsuplyn, removestoragehours=removestoragehours, residualarealist=residualarealist)
    addPowerUpperSlack!(longobjects)

    longprobmethod = parse_methods(db.settings["problems"]["prognosis"]["long"]["prob"])
    longprob = buildprob(longprobmethod, longobjects)

    # Create med problems
    medobjects = make_obj(db.input.progelements, mhh, mph)
    simplify!(medobjects; aggzone=aggzone, aggsupplyn=aggsuplyn, removestoragehours=removestoragehours, residualarealist=residualarealist)
    addPowerUpperSlack!(medobjects)

    medprobmethod = parse_methods(db.settings["problems"]["prognosis"]["med"]["prob"])
    medprob = buildprob(medprobmethod, medobjects)

    # TODO: Possible to improve? Not needed to make endvaluesobj and store it in probobjects
    medstorages = getstorages(getobjects(medprob))
    medendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
    medendvaluesobj = EndValues(medendvaluesid, medstorages) # initialize endvalues object
    push!(medprob.objects, medendvaluesobj) # push end values object to med problem objects

    # Create short problems
    shortobjects = make_obj(db.input.progelements, shh, sph)
    simplify!(shortobjects; aggzone=aggzone, removestartup=false)
    addPowerUpperSlack!(shortobjects)

    shortprobmethod = parse_methods(db.settings["problems"]["prognosis"]["short"]["prob"])
    shortprob = buildprob(shortprobmethod, shortobjects)

    shortstorages = getstorages(getobjects(shortprob))
    shortendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
    shortendvaluesobj = EndValues(shortendvaluesid, longtermstorages)
    push!(shortprob.objects, shortendvaluesobj)

    nonstoragestates = getnonstoragestatevariables(shortprob.objects) 

    db.ppp[scenarioix] = PricePrognosisProblem(longprob, medprob, shortprob, nonstoragestates)
    return
end

function make_obj(elements::Vector{DataElement}, hydro_horizon::Horizon, power_horizon::Horizon; validate::Bool=false)
    elements1 = copy(elements)

    set_horizon!(elements1, "Power", power_horizon)
    set_horizon!(elements1, "Battery", power_horizon)
    set_horizon!(elements1, "Hydro", hydro_horizon)

    modelobjects = getmodelobjects(elements1; validate=validate)

    return modelobjects
end

function simplify!(modelobjects::Dict; aggzone::Dict=Dict(), removestartup::Bool=true, removetransmissionramping::Bool=true, aggsupplyn::Int=0, removestoragehours::Int=0, residualarealist::Vector=[])
    # Aggregate price areas and add power balance slack variables
    # For the new area FRACHE, the transmission line FRA-CHE is transformed into a demand based on the loss and utilization of the line
    aggzonedict = Dict()
    for (k,v) in aggzone
        aggzonedict[Id(BALANCE_CONCEPT,"PowerBalance_" * k)] = [modelobjects[Id(BALANCE_CONCEPT,"PowerBalance_" * vv)] for vv in v]
    end
    aggzone!(modelobjects, aggzonedict)

    # Start-up-costs are not compatible with aggregatesupplycurve! or AdaptiveHorizon
    removestartup && remove_startupcosts!(modelobjects)

    # Transmissionramping is not compatible with aggregatesupplycurve! or AdaptiveHorizon and slows down short problem
    removetransmissionramping && remove_transmissionramping!(modelobjects)

    # Aggregate all simple plants (only connected to power market, mostly thermal) for each area into 4 equivalent plants
    aggsupplyn > 0 && aggregatesupplycurve!(modelobjects, aggsupplyn)

    # Short-term storage systems are only needed when the horizon is fine 
    removestoragehours > 0 && removestoragesystems!(modelobjects, Hour(removestoragehours))

    # Only calculate AdaptiveHorizon based on residual loads in these areas
    length(residualarealist) > 0 && residualloadareas!(modelobjects, residualarealist)
end

function solve_ppp(T, t, delta, stepnr, skipmed, thiscore)
    update_startstates_ppp!(stepnr, t) # TODO: A bit uncessesary to to update startstates for long and med if skipmed != 0
    skipmed.value == 0 && update_endstates_longppp()
    solve_local_ppp(t, skipmed)
    skipmed.value == 0 && syncronize_horizons(thiscore)
    return
end

function solve_local_ppp(t, skipmed)
    db = get_local_db()
    for (scenix, p) in db.ppp
        if skipmed.value == 0
            update!(p.longprob, t)
            solve!(p.longprob)

            update!(p.medprob, t)
            lhh = db.horizons[(scenix, "long", "Hydro")]
            mhh = db.horizons[(scenix, "med", "Hydro")]
            transferduals!(p.longprob, lhh, p.medprob, mhh)
            solve!(p.medprob)
        end

        update!(p.shortprob, t)
        if skipmed.value == 0 # cannot update if medprob not updated. Assume reuse of watervalues not important for short. TODO: Solve medprob at every step? Split second week in 2 day intervals? Don't reuse watervalues?
            shh = db.horizons[(scenix, "short", "Hydro")]
            transferduals!(p.medprob, mhh, p.shortprob, shh)
        end
        solve!(p.shortprob)

        sph = db.horizons[(scenix, "short", "Power")]
        update_nonstoragestates!(p, db, sph)
    end
end

function update_startstates_ppp(stepnr, t)
    db = get_local_db()

    if stepnr == 1
        get_startstates_ppp(db, t)
    else
        if stepnr == db.stepnr_startstates
            get_startstates_ppp!(db)
        else
            get_startstates_ppp_from_cp!(db)
            db.stepnr_startstates = stepnr
        end
    end
    for p in db.ppp
        setstartstates!(p.longprob, db.startstates)
        setstartstates!(p.medprob, db.startstates)
        setstartstates!(p.shortprob, db.startstates)
    end
end

function get_startstates_ppp!(db)
    dummyprogstorages = getstorages(db.dummyprogobjects)
    getstartstates!(db.startstates, db.settings["problems"], "prognosis", db.input.dataset, db.dummyprogobjects, dummyprogstorages, t)
    startstates_max!(dummyprogstorages, t, db.startstates)
    return
end

function get_startstates_ppp_from_cp!(db)
    f = @spawnat db.cp_core get_startstates_from_cp()
    startstates_cp = fetch(f)

    for (k, v) in startstates_cp
        db.startstates[k] = v
    end
    return 
end

function update_endstates_longppp()
    db = get_local_db()

    for p in db.ppp
        setstartstates!(p.longprob, db.startstates)
    end
end

function setstartstates!(p::Prob, startstates::Dict{String, Float64})
    storages = getstorages(getobjects(p))
    setstartstates!(p, storages, startstates)
    return
end

function getstartstates!(startstates::Dict, problemsconfig::Dict, problem::String, dataset::Dict, objects::Dict, storages::Vector, tnormal::ProbTime)
    startstorages = problemsconfig[problem]["startstorages"]
    if startstorages["function"] == "percentages"
        shorttermstorages = getshorttermstorages(collect(values(objects)), Hour(problemsconfig["shorttermstoragecutoff_hours"]))
        longtermstorages = setdiff(storages, shorttermstorages)
        merge!(startstates, getstartstoragepercentage(shorttermstorages, tnormal, startstorages["shortpercentage"]))
        merge!(startstates, getstartstoragepercentage(longtermstorages, tnormal, startstorages["longpercentage"]))
    elseif startstorages["function"] == "percentage"
        merge!(startstates, getstartstoragepercentage(storages, tnormal, startstorages["percentage"]))
    elseif haskey(dataset, startstorages["function"])
        merge!(startstates, [startstorages["function"]])
    end
end

# Dual values from long problem used as end values for med problem
function transferduals!(giverprob, giverhorizon, takerprob, takerhorizon)
    storages = getstorages(getobjects(takerprob))
    period = getendperiodfromduration(giverhorizon, getduration(takerhorizon)) # which period in long problem correspond to end period in medium problem
    endvalues = getinsideduals(giverprob, storages, period) # get dual values from long problem at period which correspond to end period in medium problem
    
    endvaluesobj = takerprob.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), takerprob.objects)]
    updateendvalues!(takerprob, endvaluesobj, endvalues) # update end values in problem object and in problem formulation
    return
end

# Market clearing problem uses end state values from short problem for non-storage state variables,
# TODO: Only getoutgoingstates!, init should handle changeendtoinsidestates!
function update_nonstoragesstates!(shortprob, db, sph, stepnr)
    if stepnr == 1
        clearingperiod = getendperiodfromduration(sph, db.input.steplength) # which period in short problem correspond to end period in market clearing problem
        changeendtoinsidestates!(shortprob, shortprob.short_nonstoragestates, clearingperiod) # change outgoing state variable to outgoing state in market clearing problem, and collect value
    else
        getoutgoingstates!(shortprob, shortprob.short_nonstoragestates)
    end
end

function syncronize_horizons(thiscore)
    db = get_local_db()

    owner_scenarios = [s for (s, c) in db.ppp_dist if c == thiscore]

    for ((this_scen, term, commodity), horizon) in db.horizons
        if !(this_scen in owner_scenarios)
            continue
        end

        changes = getchanges(horizon)

        if length(changes) > 0
            @sync for (other_scen, other_core) in db.ppp_dist
                if !(other_scen in owner_scenarios)
                    @spawnat other_core transfer_horizon_changes(other_scen, term, commodity, changes)
                end
            end        
        end
    end
    return
end

function transfer_horizon_changes(s::ScenarioIx, t::TermName, c::CommodityName, changes)
    db = get_local_db()
    h = db.horizons[(s, t, c)]
    setchanges!(h, changes)
    return
end

