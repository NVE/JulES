# Functions for making prob from main elements and horizon settings
function make_prob(probmethod::ProbMethod, elements::Vector{DataElement}, horizonduration::Millisecond, hydroperiodduration::Millisecond, rhsdata::AdaptiveHorizonData, method::AdaptiveHorizonMethod, clusters::Int, unit_duration::Millisecond, simplify::String)

    hydroperiods = ceil(Int64, horizonduration/hydroperiodduration)
    hydro_horizon = SequentialHorizon(hydroperiods, hydroperiodduration)
    power_horizon = AdaptiveHorizon(clusters, unit_duration, rhsdata, method, hydroperiods, hydroperiodduration)

    return _make_prob(probmethod, elements, hydro_horizon, power_horizon, simplify)
end

function make_prob(probmethod::ProbMethod, elements::Vector{DataElement}, horizonduration::Millisecond, hydroperiodduration::Millisecond, powerparts::Int, simplify::String)

    hydroperiods = ceil(Int64, horizonduration/hydroperiodduration)
    hydro_horizon = SequentialHorizon(hydroperiods, hydroperiodduration)
    power_horizon = SequentialHorizon(hydro_horizon, powerparts)

    return _make_prob(probmethod, elements, hydro_horizon, power_horizon, simplify)
end

function _make_prob(probmethod::ProbMethod, elements::Vector{DataElement}, hydro_horizon::Horizon, power_horizon::Horizon, simplify::String)
    elements1 = copy(elements)
    set_horizon!(elements1, "Power", power_horizon)
    set_horizon!(elements1, "Battery", power_horizon)
    set_horizon!(elements1, "Hydro", hydro_horizon)

    modelobjects = getmodelobjects(elements1)
    
    simplify == "long" && simplify!(modelobjects; aggsupplyn=4, removestoragehours=10, residualarealist=["DEU","NLDBEL","GBR","NOS","NON","SEN","DMK"]) # TODO: Replace with user settings
    simplify == "short" && simplify!(modelobjects; removestartup=false, residualarealist=["DEU","NLDBEL","GBR","NOS","NON","SEN","DMK"])
    addPowerUpperSlack!(modelobjects)

    prob = buildprob(probmethod, modelobjects)

    return prob, hydro_horizon, power_horizon
end

# Simplify modelobjects
function simplify!(modelobjects::Dict; aggzone::Bool=true, removestartup::Bool=true, removetransmissionramping::Bool=true, aggsupplyn::Int=0, removestoragehours::Int=0, residualarealist::Vector=[])
    # Aggregate price areas and add power balance slack variables
    # For the new area FRACHE, the transmission line FRA-CHE is transformed into a demand based on the loss and utilization of the line
    if aggzone
         aggzoneareadict = Dict("NLDBEL" => ["NLD","HUB_NLD","BEL","HUB_BEL"],  # TODO: Replace with user settings
        "FRACHE" => ["FRA","CHE"],
        "AUTCZE" => ["AUT","CZE"],
        "BAL" => ["LTU","LVA","EST","HUB_OST"],
        "DMK" => ["DK1","HUB_DK1","DK2","HUB_DK2"],
        "NOS" => ["NO1","NO2","NO5"],
        "NON" => ["NO3","NO4"],
        "SEN" => ["SE1","SE2"],
        "SES" => ["SE3","SE4"])
        aggzonedict = Dict()
        for (k,v) in aggzoneareadict
            aggzonedict[Id(BALANCE_CONCEPT,"PowerBalance_" * k)] = [modelobjects[Id(BALANCE_CONCEPT,"PowerBalance_" * vv)] for vv in v]
        end

        aggzone!(modelobjects, aggzonedict)
    end

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
function prognosis_init!(probmethods::Vector, elements::Vector{DataElement}, longinput::Tuple, medinput::Tuple, shortinput::Tuple, tnormal::ProbTime, tphasein::ProbTime, phaseinoffset::Millisecond, medprice::Dict, shortprice::Dict)
    longprob, lhh, lph  = make_prob(probmethods[1], elements, longinput...)
    
    longstorages = getstorages(getobjects(longprob))
    
    setstartstoragepercentage!(longprob, longstorages, tnormal, 65.0) # replace with user settings
    setendstoragepercentage!(longprob, longstorages, tnormal, 65.0)
    
    update!(longprob, tnormal)
    solve!(longprob)

    medprob, mhh, mph = make_prob(probmethods[1], elements, medinput...)
    
    medstorages = getstorages(getobjects(medprob))
    
    setstartstoragepercentage!(medprob, medstorages, tnormal, 65.0) # TODO: Replace with user settings
    
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
    
    shortprob, shh, sph = make_prob(probmethods[1], elements, shortinput...)
    
    shorttermstorages = getshorttermstorages(getobjects(shortprob), Hour(10))
    allstorages = getstorages(getobjects(shortprob))
    longtermstorages = setdiff(allstorages, shorttermstorages)
    
    setstartstoragepercentage!(shortprob, shorttermstorages, tnormal, 50.0)
    setendstoragepercentage!(shortprob, shorttermstorages, tnormal, 50.0)
    setstartstoragepercentage!(shortprob, longtermstorages, tnormal, 65.0)
    
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
function pl_prognosis_init!(probmethods::Vector, probs::Tuple{DArray, DArray, DArray}, allinput::Tuple, longinput::Tuple, medinput::Tuple, shortinput::Tuple, output::Tuple)
    (numcores, elements, scenarios, phaseinoffset) = allinput
    (longprobs, medprobs, shortprobs) = probs
    (medprices, shortprices, medendvaluesobjs, nonstoragestates) = output
    
    # Execute each scenario in parallel on different cores
    @sync @distributed for core in 1:(numcores-1)

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
                probs = prognosis_init!(probmethods, elements, longinput, medinput, shortinput, tnormal, tphasein, phaseinoffset, medprice[localix], shortprice[localix])
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
        
        longperiod = getendperiodfromduration(lhh, getduration(mhh))
        medendvalues = getinsideduals(longprob, medstorages, longperiod)
        medendvaluesobj = medprob.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"MedEndValue"), medprob.objects)]
        updateendvalues!(medprob, medendvaluesobj, medendvalues)

        update!(medprob, tphasein)
        solve!(medprob)

        updateareaprices!(medprice, medprob, mph, tnormal)
    end
    
    shorttermstorages = getshorttermstorages(getobjects(shortprob), Hour(10))
    allstorages = getstorages(getobjects(shortprob))
    longtermstorages = setdiff(allstorages, shorttermstorages)
    nonstorageobjects = getnonstorageobjects(getobjects(shortprob))
    
    setstartstates!(shortprob, shorttermstorages, startstates)
    setendstates!(shortprob, shorttermstorages, startstates)
    setstartstates!(shortprob, longtermstorages, startstates)
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
    
    @sync @distributed for core in 1:(numcores-1)
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

function updateareaprices!(price::Dict, prob::Prob, horizon::Horizon, t::ProbTime)
    price["steprange"] = getscenariosteprange(horizon, t)

    for (i, name) in enumerate(price["names"])
        id = Id("Balance", "PowerBalance_" * name)
        price["matrix"][:,i] = -getconduals(prob, horizon, id)
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
