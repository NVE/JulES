# Functions for making modelobjects from main elements and horizon settings
function make_obj(elements::Vector{DataElement}, horizonduration::Millisecond, hydroperiodduration::Millisecond, rhsdata::AdaptiveHorizonData, method::AdaptiveHorizonMethod, clusters::Int, unit_duration::Millisecond)
    hydroperiods = ceil(Int64, horizonduration/hydroperiodduration)
    hydro_horizon = SequentialHorizon(hydroperiods, hydroperiodduration)
    power_horizon = AdaptiveHorizon(clusters, unit_duration, rhsdata, method, hydroperiods, hydroperiodduration)

    return _make_obj(elements, hydro_horizon, power_horizon)
end

function make_shrinkable_obj(elements::Vector{DataElement}, horizonduration::Millisecond, hydroperiodduration::Millisecond, rhsdata::AdaptiveHorizonData, method::AdaptiveHorizonMethod, clusters::Int, unit_duration::Millisecond, startafter::Millisecond, shrinkatleast::Millisecond, minperiod::Millisecond)
    hydroperiods = ceil(Int64, horizonduration/hydroperiodduration)
    hydro_horizon = ShrinkableHorizon(SequentialHorizon(1, startafter, hydroperiods, hydroperiodduration), startafter, shrinkatleast, minperiod)
    power_horizon = ShrinkableHorizon(AdaptiveHorizon(clusters, unit_duration, rhsdata, method, 1, startafter, hydroperiods, hydroperiodduration), startafter, shrinkatleast, minperiod)

    return _make_obj(elements, hydro_horizon, power_horizon)
end

function make_obj(elements::Vector{DataElement}, horizonduration::Millisecond, hydroperiodduration::Millisecond, powerparts::Int)
    hydroperiods = ceil(Int64, horizonduration/hydroperiodduration)
    hydro_horizon = SequentialHorizon(hydroperiods, hydroperiodduration)
    power_horizon = SequentialHorizon(hydro_horizon, powerparts)

    return _make_obj(elements, hydro_horizon, power_horizon)
end

function _make_obj(elements::Vector{DataElement}, hydro_horizon::Horizon, power_horizon::Horizon)
    elements1 = copy(elements)
    set_horizon!(elements1, "Power", power_horizon)
    set_horizon!(elements1, "Battery", power_horizon)
    set_horizon!(elements1, "Hydro", hydro_horizon)

    modelobjects = getmodelobjects(elements1)

    return modelobjects, hydro_horizon, power_horizon
end

# Simplify modelobjects
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

# Initialize price prognosis problems for specific scenario
function prognosis_init!(probmethods::Vector, objects::Tuple, horizons::Tuple, startstates::Dict{String, Float64}, tnormal::ProbTime, tphasein::ProbTime, phaseinoffset::Millisecond, medprice::Dict, shortprice::Dict)
    (longobjects, medobjects, shortobjects) = objects
    (lhh, lph, mhh, mph, shh, sph) = horizons
    
    longprob = buildprob(probmethods[1], copy(longobjects))

    longstorages = getstorages(getobjects(longprob))
    setstartstates!(longprob, longstorages, startstates)
    setendstates!(longprob, longstorages, startstates)
    
    update!(longprob, tnormal)
    solve!(longprob)

    medprob = buildprob(probmethods[2], copy(medobjects))
    
    medstorages = getstorages(getobjects(medprob))
    setstartstates!(medprob, medstorages, startstates)
    
    # Dual values from long problem used as end values for med problem, initialize
    longperiod = getendperiodfromduration(lhh, getduration(mhh)) # which period in long problem correspond to end period in medium problem
    medendvalues = getinsideduals(longprob, medstorages, longperiod) # get dual values from long problem at period which correspond to end period in medium problem
    medendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"MedEndValue")
    medendvaluesobj = EndValues(medendvaluesid, medstorages) # initialize endvalues object
    push!(medprob.objects, medendvaluesobj) # push end values object to med problem objects
    updateendvalues!(medprob, medendvaluesobj, medendvalues) # update end values in problem object and in problem formulation

    update!(medprob, tphasein)
    solve!(medprob)

    # Collect prices that will be used in stochastic subsystem problems
    getareaprices!(medprice, medprob, mph, tnormal)
    
    shortprob= buildprob(probmethods[3], copy(shortobjects))
    
    shortstorages = getstorages(getobjects(shortprob))
    shorttermstorages = getshorttermstorages(getobjects(shortprob), Hour(10))
    longtermstorages = setdiff(shortstorages, shorttermstorages)
    setstartstates!(shortprob, shortstorages, startstates) # set startstates for all storages
    setendstates!(shortprob, shorttermstorages, startstates) # set endstates for shortterm storages
    
    # Dual values from med problem used as end values for short problem, initialize
    medperiod = getendperiodfromduration(mhh, getduration(shh))
    shortendvalues = getinsideduals(medprob, longtermstorages, medperiod)
    shortendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"ShortEndValue")
    shortendvaluesobj = EndValues(shortendvaluesid, longtermstorages)
    push!(shortprob.objects, shortendvaluesobj)
    updateendvalues!(shortprob, shortendvaluesobj, shortendvalues)

    update!(shortprob, tphasein)
    solve!(shortprob)

    # Market clearing problem uses end state values from short problem for non-storage state variables, 
    clearingperiod = getendperiodfromduration(sph, phaseinoffset) # which period in short problem correspond to end period in market clearing problem
    nonstoragestates = getnonstoragestatevariables(shortprob.objects) # get non-storage state variables (e.g. variables used for thermal power ramping or startup costs)
    changeendtoinsidestates!(shortprob, nonstoragestates, clearingperiod) # change outgoing state variable to outgoing state in market clearing problem, and collect value

    # Collect prices that will be used in stochastic subsystem problems
    getareaprices!(shortprice, shortprob, sph, tnormal)

    return (longprob, medprob, shortprob, medendvaluesobj, nonstoragestates)
end

# Initialize price prognosis models for all scenarios in parallel
function pl_prognosis_init!(probmethods::Vector, probs::Tuple{DArray, DArray, DArray}, objects::Tuple, horizons::Tuple, input::Tuple, output::Tuple)
    (longprobs, medprobs, shortprobs) = probs
    (numcores, scenarios, phaseinoffset, startstates) = input
    (medprices, shortprices, medendvaluesobjs, nonstoragestates) = output
    
    # Execute each scenario in parallel on different cores
    @sync @distributed for core in 1:max(numcores-1,1)

        # Local version of distributed arrays only consist of elements that are assignet to this specific core
        scenario = localpart(scenarios)
        longprob = localpart(longprobs)
        medprob = localpart(medprobs)
        shortprob = localpart(shortprobs)
        medprice = localpart(medprices)
        shortprice = localpart(shortprices)
        medendvaluesobj = localpart(medendvaluesobjs)
        nonstoragestate = localpart(nonstoragestates)

        localix = 0
        for range in localindices(longprobs) # for each scenariorange on this specific core
            for ix in range # for each scenario on this specific core
                localix += 1
                (tnormal, tphasein, scen) = scenario[localix]
                probs = prognosis_init!(probmethods, objects, horizons, startstates, tnormal, tphasein, phaseinoffset, medprice[localix], shortprice[localix])
                longprob[localix], medprob[localix], shortprob[localix], medendvaluesobj[localix], nonstoragestate[localix] = probs
            end
        end
    end
end

# Run price prognosis model for a specific scenario
function prognosis!(longprob::Prob, medprob::Prob, shortprob::Prob, medprice::Dict, shortprice::Dict, nonstoragestates::Dict{StateVariableInfo, Float64}, startstates::Dict, tnormal::ProbTime, tphasein::ProbTime, skipmed::Millisecond)

    # Collect hydro and power horizons for each problem
    mh = medprob.horizons
    mhp = [getnumperiods(h) for h in mh]
    sh = shortprob.horizons
    shp = [getnumperiods(h) for h in sh]
    if skipmed.value == 0
        lh = longprob.horizons
        lhp = [getnumperiods(h) for h in lh]
        lhh = lh[argmin(lhp)]
        mph = mh[argmax(mhp)]
    end
    mhh = mh[argmin(mhp)]
    sph = sh[argmax(shp)]
    shh = sh[argmin(shp)]

    # Skipmed inidcates if we reuse water values for this time step, and therefore does not have to run the long and medium problems
    if skipmed.value == 0
        longstorages = getstorages(getobjects(longprob))
        
        setstartstates!(longprob, longstorages, startstates)
        setendstates!(longprob, longstorages, startstates)
        
        update!(longprob, tnormal)
        solve!(longprob)
        
        medstorages = getstorages(getobjects(medprob))
        setstartstates!(medprob, medstorages, startstates)
        
        update!(medprob, tphasein) # horizons needs to be updated before we can calculate longperiod

        longperiod = getendperiodfromduration(lhh, getduration(mhh))
        medendvalues = getinsideduals(longprob, medstorages, longperiod)
        medendvaluesobj = medprob.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"MedEndValue"), medprob.objects)]
        updateendvalues!(medprob, medendvaluesobj, medendvalues)

        solve!(medprob)

        updateareaprices!(medprice, medprob, mph, tnormal)
    end
    
    shorttermstorages = getshorttermstorages(getobjects(shortprob), Hour(10))
    allstorages = getstorages(getobjects(shortprob))
    longtermstorages = setdiff(allstorages, shorttermstorages)
    nonstorageobjects = getnonstorageobjects(getobjects(shortprob))
    
    setstartstates!(shortprob, allstorages, startstates)
    setendstates!(shortprob, shorttermstorages, startstates)
    setstartstates!(shortprob, nonstorageobjects, startstates) # NB! Assumes same resolution in shortprob as market clearing
    
    medperiod = getendperiodfromduration(mhh, getduration(shh))
    shortendvalues = getinsideduals(medprob, longtermstorages, medperiod)
    shortendvaluesobj = shortprob.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"ShortEndValue"), shortprob.objects)]
    updateendvalues!(shortprob, shortendvaluesobj, shortendvalues)

    update!(shortprob, tphasein)
    solve!(shortprob)

    getoutgoingstates!(shortprob, nonstoragestates)

    updateareaprices!(shortprice, shortprob, sph, tnormal)
end

# Run price prognosis models in parallel
function pl_prognosis!(numcores::Int, longprobs::DArray, medprobs::DArray, shortprobs::DArray, medprices::DArray, shortprices::DArray, nonstoragestates::DArray, startstates::Dict, scenarios::DArray, skipmed::Millisecond)
    
    @sync @distributed for core in 1:max(numcores-1,1)
        scenario = localpart(scenarios)
        longprob = localpart(longprobs)
        medprob = localpart(medprobs)
        shortprob = localpart(shortprobs)
        medprice = localpart(medprices)
        shortprice = localpart(shortprices)
        nonstoragestate = localpart(nonstoragestates)

        localix = 0
        for range in localindices(longprobs)
            for ix in range
                localix += 1
                (tnormal, tphasein, scen) = scenario[localix]
                prognosis!(longprob[localix], medprob[localix], shortprob[localix], medprice[localix], shortprice[localix], nonstoragestate[localix], startstates, tnormal, tphasein, skipmed)
            end
        end
    end
end

# Collect prices in each period for each area
# - Only supports AdaptiveHorizons
function getareaprices!(price::Dict, prob::Prob, horizon::Horizon, t::ProbTime)
    price["steprange"] = getscenariosteprange(horizon, t)
    price["names"] = []

    for obj in getobjects(prob)
        if obj isa BaseBalance
            if getinstancename(getid(getcommodity(obj))) == "Power"
                push!(price["names"], split(getinstancename(getid(obj)),"PowerBalance_")[2])
            end
        end
    end

    price["matrix"] = zeros(Float64, length(price["steprange"]), length(price["names"]))

    for (i, name) in enumerate(price["names"])
        id = Id("Balance", "PowerBalance_" * name)
        price["matrix"][:,i] = -getconduals(prob, horizon, id)
    end
end

function updateareaprices!(price::Dict, prob::Prob, horizon::Horizon, t::ProbTime) # could specify this to each case (Shrinkable, Sequential, Adaptive)
    (nr,nc) = size(price["matrix"])
    numconduals = -1

    for (i, name) in enumerate(price["names"])
        id = Id("Balance", "PowerBalance_" * name)
        conduals = -getconduals(prob, horizon, id)
        price["matrix"][(nr - length(conduals) + 1):end,i] .= conduals

        if i == 1
            numconduals = length(conduals)
        end
    end

    # For SequentialPeriod we could check if horizonduration equals steprange duration
    if numconduals == nr
        price["steprange"] = getscenariosteprange(horizon, t)
    end
end

# Get steprange that represent a horizons periods
function getscenariosteprange(horizon::SequentialHorizon, start::ProbTime)
    @assert length(horizon.periods.data) == 1 # error if different period durations in horizon
    scenariostart = getscenariotime(getstarttime(horizon, 1, start))
    scenariostop = getscenariotime(getstarttime(horizon, getnumperiods(horizon), start))
    periodduration = getduration(gettimedelta(horizon, 1))
    return StepRange(scenariostart, periodduration, scenariostop)
end

function getscenariosteprange(horizon::AdaptiveHorizon, start::ProbTime)
    scenariostart = getscenariotime(start)
    scenariostop = getscenariotime(start + getduration(horizon) - horizon.unit_duration)
    periodduration = horizon.unit_duration
    return StepRange(scenariostart, periodduration, scenariostop)
end

getscenariosteprange(horizon::ShrinkableHorizon, start::ProbTime) = getshrinkablescenariosteprange(horizon.subhorizon, start)

getshrinkablescenariosteprange(horizon::Horizon, start::ProbTime) = error("Not supported")
# Could make steprange for SequentialPeriod based on minperiod
getshrinkablescenariosteprange(horizon::AdaptiveHorizon, start::ProbTime) = getscenariosteprange(horizon, start)


# Collect dual values for a modelobject id for all periods - TODO: Move to TuLiPa
getconduals(prob::Prob, horizon::SequentialHorizon, id::Id) = [getcondual(prob, id, s) for s in 1:getnumperiods(horizon)]

function getconduals(prob::Prob, horizon::AdaptiveHorizon, id::Id)
    values = [getcondual(prob, id, s) for s in 1:getnumperiods(horizon)]

    points = []
    for t in 1:getnumperiods(horizon)
        value = values[t]
        timedelta = gettimedelta(horizon, t)
        startduration = getstartduration(horizon, t)
        for unit_range in timedelta.units
            for i in unit_range
                point_dur = startduration + ((i-1) * timedelta.unit_duration)
                point = (point_dur, value)
                push!(points, point)
            end
        end
    end
    sort!(points)
    return [p[2] for p in points]
end

getconduals(prob::Prob, horizon::ShrinkableHorizon, id::Id) = getshrinkableconduals(prob, horizon.subhorizon, id)

getshrinkableconduals(prob::Prob, horizon::Horizon, id::Id) = error("Horizon not supported") 
# Could split conduals for SequentialHorizon into minperiod
getshrinkableconduals(prob::Prob, horizon::AdaptiveHorizon, id::Id) = getconduals(prob, horizon, id)
