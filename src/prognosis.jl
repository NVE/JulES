# Function for making prob from main elements and horizon setting
function make_prob(scenario::Int, t::ProbTime, elements::Vector{DataElement}, horizonstart::Int, horizonend::DateTime, hydroperiodduration::Millisecond, rhsdata::AdaptiveHorizonData, method::AdaptiveHorizonMethod, clusters::Int, unit_duration::Millisecond, simplify::String)

    offset = TimeDeltaOffset(MsTimeDelta(getisoyearstart(horizonstart + scenario - 1) - getisoyearstart(horizonstart)))

    hydro_horizon = SequentialHorizon(getisoyearstart(horizonstart), horizonend, hydroperiodduration; offset)
    power_horizon = AdaptiveHorizon(clusters, unit_duration, rhsdata, method, getisoyearstart(horizonstart), horizonend, hydroperiodduration; offset)

    elements1 = copy(elements)
    elements1 = set_horizon!(elements1, "Power", power_horizon)
    set_horizon!(elements1, "Battery", power_horizon)
    set_horizon!(elements1, "Hydro", hydro_horizon)

    modelobjects = getmodelobjects(elements1)
    simplify == "long" && simplify!(modelobjects; aggsupplyn=4, removestoragehours=10, residualarealist=["DEU","NLDBEL","GBR","NOS","NON","SEN","DMK"])
    simplify == "short" && simplify!(modelobjects; removestartup=false, residualarealist=["DEU","NLDBEL","GBR","NOS","NON","SEN","DMK"])
    addPowerUpperSlack!(modelobjects)
    
    # model = Model(HiGHS.Optimizer)
    # set_silent(model)
    # prob = JuMP_Prob(modelobjects, model)
    prob = HiGHS_Prob(modelobjects)
    
    return prob, hydro_horizon, power_horizon
end

# Function for making prob from main elements and horizon setting
function make_prob(scenario::Int, t::ProbTime, elements::Vector{DataElement}, horizonstart::Int, horizonend::DateTime, hydroperiodduration::Millisecond, powerparts::Int, simplify::String)

    offset = TimeDeltaOffset(MsTimeDelta(getisoyearstart(horizonstart + scenario - 1) - getisoyearstart(horizonstart)))

    hydro_horizon = SequentialHorizon(getisoyearstart(horizonstart), horizonend, hydroperiodduration; offset)
    power_horizon = SequentialHorizon(hydro_horizon, powerparts)

    elements1 = copy(elements)
    elements1 = set_horizon!(elements1, "Power", power_horizon)
    set_horizon!(elements1, "Battery", power_horizon)
    set_horizon!(elements1, "Hydro", hydro_horizon)

    modelobjects = getmodelobjects(elements1)
    
    simplify == "long" && simplify!(modelobjects; aggsupplyn=4, removestoragehours=10, residualarealist=["DEU","NLDBEL","GBR","NOS","NON","SEN","DMK"])
    simplify == "short" && simplify!(modelobjects; removestartup=false, residualarealist=["DEU","NLDBEL","GBR","NOS","NON","SEN","DMK"])
    addPowerUpperSlack!(modelobjects)

    # model = Model(HiGHS.Optimizer)
    # set_silent(model)
    # prob = JuMP_Prob(modelobjects, model)
    prob = HiGHS_Prob(modelobjects)

    return prob, hydro_horizon, power_horizon
end

# Function updating and solving prob
function update_solve!(prob::Prob, t::ProbTime)
    
    update!(prob, t)

    solve!(prob)
    #     display(string(i, ": ", getobjectivevalue(prob)))
    return prob
end

function prognosis_init!(elements, horizonstart, longinput, medinput, shortinput, scenario, tnormal, tphasein, phaseinoffsetdays, medprice, shortprice)
    longprob, lhh, lph  = make_prob(scenario, tnormal, elements, horizonstart, longinput...)
    
    longstorages = getstorages(getobjects(longprob))
    
    setstartstoragepercentage!(longprob, longstorages, tnormal, 65)
    setendstoragepercentage!(longprob, longstorages, tnormal, 65)
    
    update_solve!(longprob, tnormal)

    medprob, mhh, mph = make_prob(scenario, tphasein, elements, horizonstart, medinput...)
    
    medstorages = getstorages(getobjects(medprob))
    
    setstartstoragepercentage!(medprob, medstorages, tnormal, 65)
    
    longperiod = getendperiodfromduration(lhh, getduration(mhh))
    medendvalues = getinsideduals(longprob, medstorages, longperiod)
    medendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"MedEndValue")
    medendvaluesobj = EndValues(medendvaluesid, medstorages)
    push!(medprob.objects, medendvaluesobj)
    updateendvalues!(medprob, medendvaluesobj, medendvalues)

    # # Alternative with transfering end storage instead of end storage value (might give non-optimal solution)
    # medstates = getstatevariables(medstorages)
    # longperiod = getendperiodfromduration(lhh, getduration(mhh))
    # getinsidestates!(longprob, medstates, longperiod)
    # setoutgoingstates!(medprob, medstates)

    update_solve!(medprob, tphasein)

    getareaprices!(medprice, medprob, mph, tnormal)
    
    shortprob, shh, sph = make_prob(scenario, tphasein, elements, horizonstart, shortinput...)
    
    shorttermstorages = getshorttermstorages(getobjects(shortprob), Hour(10))
    allstorages = getstorages(getobjects(shortprob))
    longtermstorages = setdiff(allstorages, shorttermstorages)
    
    setstartstoragepercentage!(shortprob, shorttermstorages, tnormal, 50)
    setendstoragepercentage!(shortprob, shorttermstorages, tnormal, 50)
    setstartstoragepercentage!(shortprob, longtermstorages, tnormal, 65)
    
    medperiod = getendperiodfromduration(mhh, getduration(shh))
    shortendvalues = getinsideduals(medprob, longtermstorages, medperiod)
    shortendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"ShortEndValue")
    shortendvaluesobj = EndValues(shortendvaluesid, longtermstorages)
    push!(shortprob.objects, shortendvaluesobj)
    updateendvalues!(shortprob, shortendvaluesobj, shortendvalues)

    # # Alternative with transfering end storage instead of end storage value (might give non-optimal solution)
    # shortstates = getstatevariables(longtermstorages)
    # medperiod = getendperiodfromduration(mhh, getduration(shh))
    # getinsidestates!(medprob, shortstates, medperiod)
    # setoutgoingstates!(shortprob, shortstates)

    update_solve!(shortprob, tphasein)

    clearingperiod = getendperiodfromduration(sph, Millisecond(Day(phaseinoffsetdays)))
    nonstoragestates = getnonstoragestatevariables(shortprob.objects)
    changeendtoinsidestates!(shortprob, nonstoragestates, clearingperiod)

    getareaprices!(shortprice, shortprob, sph, tnormal)

    return (longprob, medprob, shortprob, medendvaluesobj, nonstoragestates)
end

function getareaprices!(price::Dict, prob::Prob, horizon::Horizon, t::ProbTime)
    price["steprange"] = getscenariosteprange(horizon, t)
    price["names"] = []
    ids = Id[]

    for obj in getobjects(prob)
        if obj isa BaseBalance
            if getinstancename(getid(getcommodity(obj))) == "Power"
                id = getid(obj)
                push!(ids, id)
                push!(price["names"], split(getinstancename(id),"PowerBalance_")[2])
            end
        end
    end

    price["matrix"] = zeros(Float64, length(price["steprange"]), length(price["names"]))

    for (i, id) in enumerate(ids)
        price["matrix"][:,i] = getconduals(prob, horizon, id)
    end
end

function updateareaprices!(price::Dict, prob::Prob, horizon::Horizon, t::ProbTime)
    price["steprange"] = getscenariosteprange(horizon, t)

    for (i, name) in enumerate(price["names"])
        id = Id("Balance", "PowerBalance_" * name)
        price["matrix"][:,i] = getconduals(prob, horizon, id)
    end
end

function getscenariosteprange(horizon::SequentialHorizon, start::ProbTime)
    @assert length(horizon.periods.data) == 1
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

function pl_prognosis_init!(probs, allinput, longinput, medinput, shortinput, output)
    (numcores, elements, horizonstart, tnormal, tphasein, phaseinoffsetdays) = allinput
    (longprobs, medprobs, shortprobs) = probs
    (medprices, shortprices, medendvaluesobjs, nonstoragestates) = output
    
    @sync @distributed for core in 1:(numcores-1)
        longprob = localpart(longprobs)
        medprob = localpart(medprobs)
        shortprob = localpart(shortprobs)
        medprice = localpart(medprices)
        shortprice = localpart(shortprices)
        medendvaluesobj = localpart(medendvaluesobjs)
        nonstoragestate = localpart(nonstoragestates)

        localix = 0
        for range in localindices(longprobs)
            for ix in range
                localix += 1
                probs = prognosis_init!(elements, horizonstart, longinput, medinput, shortinput, ix, tnormal, tphasein, phaseinoffsetdays, medprice[localix], shortprice[localix])
                longprob[localix], medprob[localix], shortprob[localix], medendvaluesobj[localix], nonstoragestate[localix] = probs
            end
        end
    end
end

function simplify!(modelobjects; aggzone=true, removestartup=true, aggsupplyn=0, removestoragehours=0, residualarealist=[])
    # Aggregate price areas and add power balance slack variables
    # For the new area FRACHE, the transmission line FRA-CHE is transformed into a demand based on the loss and utilization of the line
    if aggzone
         aggzoneareadict = Dict("NLDBEL" => ["NLD","HUB_NLD","BEL","HUB_BEL"],
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

    # Aggregate all simple plants (only connected to power market, mostly thermal) for each area into 4 equivalent plants
    aggsupplyn > 0 && aggregatesupplycurve!(modelobjects, aggsupplyn)

    # Short-term storage systems are only needed when the horizon is fine 
    removestoragehours > 0 && removestoragesystems!(modelobjects, Hour(removestoragehours))

    # Only calculate AdaptiveHorizon based on residual loads in these areas
    length(residualarealist) > 0 && residualloadareas!(modelobjects, residualarealist)
end

function prognosis!(longprob, medprob, shortprob, medprice, shortprice, nonstoragestates, startstates, tnormal, tphasein, skipmed)

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

    if skipmed.value == 0
        longstorages = getstorages(getobjects(longprob))
        
        setstartstates!(longprob, longstorages, startstates)
        setendstates!(longprob, longstorages, startstates)
        
        update_solve!(longprob, tnormal)
        
        medstorages = getstorages(getobjects(medprob))
        
        setstartstates!(medprob, medstorages, startstates)
        
        longperiod = getendperiodfromduration(lhh, getduration(mhh))
        medendvalues = getinsideduals(longprob, medstorages, longperiod)
        medendvaluesobj = medprob.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"MedEndValue"), medprob.objects)]
        updateendvalues!(medprob, medendvaluesobj, medendvalues)

        update_solve!(medprob, tphasein)

        updateareaprices!(medprice, medprob, mph, tnormal)
    end
    
    shorttermstorages = getshorttermstorages(getobjects(shortprob), Hour(10))
    allstorages = getstorages(getobjects(shortprob))
    longtermstorages = setdiff(allstorages, shorttermstorages)
    
    setstartstates!(shortprob, shorttermstorages, startstates)
    setendstates!(shortprob, shorttermstorages, startstates)
    setstartstates!(shortprob, longtermstorages, startstates)
    
    medperiod = getendperiodfromduration(mhh, getduration(shh))
    shortendvalues = getinsideduals(medprob, longtermstorages, medperiod)
    shortendvaluesobj = shortprob.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"ShortEndValue"), shortprob.objects)]
    updateendvalues!(shortprob, shortendvaluesobj, shortendvalues)

    update_solve!(shortprob, tphasein)

    getoutgoingstates!(shortprob, nonstoragestates)

    updateareaprices!(shortprice, shortprob, sph, tnormal)
end

function pl_prognosis!(numcores, longprobs, medprobs, shortprobs, medprices, shortprices, nonstoragestates, startstates, tnormal, tphasein, skipmed)
    
    @sync @distributed for core in 1:(numcores-1)
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
                prognosis!(longprob[localix], medprob[localix], shortprob[localix], medprice[localix], shortprice[localix], nonstoragestate[localix], startstates, tnormal, tphasein, skipmed)
            end
        end
    end
end

function getmedendvaluesdict(medendvaluesobjs)
    scenario = 1
    instance = [getinstancename(getid(obj)) for obj in medendvaluesobjs[scenario].objects]
    endvalues = medendvaluesobjs[scenario].values
    medendvaluesdict = Dict(instance .=> endvalues)

    for (k,v) in detailedrescopl
        medendvaluesdict["Reservoir_" * k] = medendvaluesdict["Reservoir_" * v * "_hydro_reservoir"]
    end
    
    return medendvaluesdict
end
