using DistributedArrays

function makestochasticobjects(elements, days::Int, offset::Union{Offset,Nothing}, scenario::Int, prices::DArray, short::Bool, master::Bool)
    if short
        # Add horizons to the dataset
        battery_horizon = SequentialHorizon(days*24, Hour(1); offset)
        hydro_horizon = SequentialHorizon(days*24, Hour(1); offset)
        power_horizon = SequentialHorizon(days*24, Hour(1); offset)
    else
        if days < 7
            battery_horizon = SequentialHorizon(days, Hour(24*days); offset)
            hydro_horizon = SequentialHorizon(days, Hour(24*days); offset)
            power_horizon = SequentialHorizon(days, Hour(24*days); offset)
        else
            battery_horizon = SequentialHorizon(ceil(Int64, days/7), Hour(168); offset) # TODO: Quickfix
            hydro_horizon = SequentialHorizon(ceil(Int64, days/7), Hour(168); offset)
            power_horizon = SequentialHorizon(ceil(Int64, days/7), Hour(168); offset)
        end
    end
    push!(elements, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Power", 
            (HORIZON_CONCEPT, power_horizon)))
    push!(elements, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Hydro", 
            (HORIZON_CONCEPT, hydro_horizon)))
    push!(elements, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Battery", 
            (HORIZON_CONCEPT, battery_horizon)))
    
    price = prices[scenario]

    push!(elements, DataElement(TIMEINDEX_CONCEPT, "RangeTimeIndex", "ShortTermTimeIndex", price["steprange"]))
    
    for (i, area) in enumerate(price["names"])
        push!(elements, getelement(TIMEVALUES_CONCEPT, "VectorTimeValues", "ValuesPrice_" * area,
            ("Vector", price["matrix"][:, i])))
        push!(elements, getelement(TIMEVECTOR_CONCEPT, "InfiniteTimeVector", "ProfilePrice_" * area,
            (TIMEINDEX_CONCEPT, "ShortTermTimeIndex"), (TIMEVALUES_CONCEPT, "ValuesPrice_" * area)))
        push!(elements, getelement(PARAM_CONCEPT, "MeanSeriesIgnorePhaseinParam", "Price_" * area,
            ("Level", 1),
            ("Profile", "ProfilePrice_" * area)))
        push!(elements, getelement(BALANCE_CONCEPT, "ExogenBalance", "PowerBalance_" * area, 
            (COMMODITY_CONCEPT, "Power"),
            (PRICE_CONCEPT, "Price_" * area)))
    end

    modelobjects = getmodelobjects(elements)
    if short
        if master # removes spills from upper and lower storages in PHS, to avoid emptying reservoirs in master problem
            for id in keys(modelobjects)
                instance = getinstancename(id)
                if occursin("Spill", instance) && occursin("_PHS_", instance)
                    delete!(modelobjects, id)
                end
            end
        end
        storagesystems = getstoragesystems_full!(getshorttermstoragesystems(getstoragesystems(modelobjects), Hour(10)))
    else
        storagesystems = getstoragesystems_full!(getlongtermstoragesystems(getstoragesystems(modelobjects), Hour(10)))
    end

    for storagesystem in storagesystems
        sort!(storagesystem, by = x -> getinstancename(getid(x)))
    end

    sort!(storagesystems, by = x -> getinstancename(getid(first(x))))

    return storagesystems
end

function makemastersubobjects!(inputs, mastersubobjects, shorts)
    (elements, totaldays, numscen, horizonstart, phaseinoffsetdays, prices, short) = inputs

    elements1 = copy(elements)
    removeelements!(elements1, short)

    masterobjects = makestochasticobjects(copy(elements1), phaseinoffsetdays, nothing, 1, prices, short, true) # what scenario to use here? use of phasein of scenarios gives similar prices in the start of all scenarios?

    subscenarioobjects = []
    for scenario in 1:numscen
        offset = TimeDeltaOffset(MsTimeDelta(getisoyearstart(horizonstart + scenario - 1) - getisoyearstart(horizonstart)))
        push!(subscenarioobjects, makestochasticobjects(copy(elements1), totaldays*7 - phaseinoffsetdays, offset, scenario, prices, short, false))
    end

    for (i, masterobject) in enumerate(masterobjects)
        subobjects = [v[i] for v in subscenarioobjects]
        push!(mastersubobjects, (masterobject, subobjects))
        push!(shorts, short)
    end
end

function removeelements!(elements, short)
    aggzoneareadict = Dict("NLDBEL" => ["NLD","HUB_NLD","BEL","HUB_BEL"],
    "FRACHE" => ["FRA","CHE"],
    "AUTCZE" => ["AUT","CZE"],
    "BAL" => ["LTU","LVA","EST","HUB_OST"],
    "DMK" => ["DK1","HUB_DK1","DK2","HUB_DK2"],
    "NOS" => ["NO1","NO2","NO5"],
    "NON" => ["NO3","NO4"],
    "SEN" => ["SE1","SE2"],
    "SES" => ["SE3","SE4"])
    aggzonecopl = Dict()
    for (k,v) in aggzoneareadict
        for vv in v
            aggzonecopl["PowerBalance_" * vv] = "PowerBalance_" * k
        end
    end
    
    delix = []
    powerbasebalances = []
    for (i,element) in enumerate(elements)
        # BaseBalances
        if element.typename == "BaseBalance"
            if element.value["Commodity"] == "Power"
                push!(delix,i)
                push!(powerbasebalances, element.instancename)
            end
        end
        
        # Residualhints
        if element.typename == "Residualhint"
            push!(delix,i)
        end
        
        # Power balance RHSTerms
        if element.conceptname == "RHSTerm"
            if element.value["Balance"] in powerbasebalances
                push!(delix,i)
            end
        end
        
        # Power balance Arrows
        if element.conceptname == "Arrow"
            if element.value["Balance"] in keys(aggzonecopl)
                value = copy(element.value)
                value["Balance"] = aggzonecopl[element.value["Balance"]]
                elements[i] = DataElement(element.conceptname, element.typename, element.instancename, value)
            end
        end
    end
    deleteat!(elements, delix)
    
    return elements
end

# Get cutobjects
function getcutobjects(modelobjects)
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

# Initialize cuts
function initialize_cuts!(modelobjects, cutobjects, maxcuts, lb, numscen)
    # Make a cutid
    cutid = Id(BOUNDARYCONDITION_CONCEPT,"StorageCuts")
    
    # Probability of each subproblem / second stage scenario
    probabilities = [1/numscen for i in 1:numscen]
    
    # Make cut modelobject
    cuts = SimpleSingleCuts(cutid, cutobjects, probabilities, maxcuts, lb)
    push!(modelobjects, cuts)
    return cuts
end

function getstatevariables(cutobjects::Vector{Any})
    states = Dict{StateVariableInfo, Float64}()
    
    for obj in cutobjects
        for statevariable in getstatevariables(obj)
            states[statevariable] = 0.0
        end
    end
    return states
end

# Set start reservoir as a percentage of capacity
function setstartstorage!(prob, start, percentage)
    for obj in prob.objects
        if obj isa Storage

            dummydelta = MsTimeDelta(Millisecond(0))
            startreservoir = getparamvalue(getub(obj), start, dummydelta)*percentage/100
            
            ingoingstates = Dict{StateVariableInfo, Float64}()
            for statevariable in getstatevariables(obj)
                ingoingstates[statevariable] = startreservoir
            end
            # display(ingoingstates)
            
            setingoingstates!(prob, ingoingstates)
        end
    end
end

function setendstorage!(prob, endtime, percentage)
    for obj in prob.objects
        if obj isa Storage

            dummydelta = MsTimeDelta(Millisecond(0))
            endreservoir = getparamvalue(getub(obj), endtime, dummydelta)*percentage/100
            
            outgoingstates = Dict{StateVariableInfo, Float64}()
            for statevariable in getstatevariables(obj)
                outgoingstates[statevariable] = endreservoir
            end
            # display(outgoingstates)
            
            setoutgoingstates!(prob, outgoingstates)
        end
    end
end

# Transfer master problem end reservoir to subproblems (start reservoir)
function transferboundarystorage!(master, sub, states)
    states = getoutgoingstates!(master, states)
    # display(states)
    setingoingstates!(sub, states)
end

function stochastic_init(masterobjects, subobjects, short, storageinfo, numscen, lb, maxcuts, reltol, t)
    @time begin
        shortstartstorage, medstartstorage, medendvaluesdict = storageinfo
        if short
            startstorage = shortstartstorage
        else
            startstorage = medstartstorage
        end

        cutobjects = getcutobjects(masterobjects)
        cuts = initialize_cuts!(masterobjects, cutobjects, maxcuts, lb, numscen);
        states = getstatevariables(cutobjects) # state variables in master and subs for boundary reservoirs


        # master = JuMP_Prob(masterobjects, Model(HiGHS.Optimizer))
        # subs = [JuMP_Prob(subobject, Model(HiGHS.Optimizer)) for subobject in subobjects] # initialize subproblems
        master = HiGHS_Prob(masterobjects)
        subs = [HiGHS_Prob(subobject) for subobject in subobjects] # initialize subproblems

        # Init cutparameters
        cutparameters = Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}(undef, length(subs)) # preallocate for cutparameters from subproblems

        # Update master
        update!(master, t)
        setstartstorage!(master, t, startstorage)

        # Update subs
        for (i,sub) in enumerate(subs)
            update!(sub, t) # update parameters given problem start time of scenario

            if short
                setendstorage!(sub, t, startstorage) # set end reservoir
            else
                storages = getstorages(getobjects(sub))
                subendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
                subendvaluesobj = EndValues(subendvaluesid, storages)
                push!(sub.objects, subendvaluesobj)
                subendvalues = [medendvaluesdict[getinstancename(getid(obj))] for obj in storages]
                updateendvalues!(sub, subendvaluesobj, subendvalues)
            end
        end

        ub = 0
        while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)
            solve!(master)
            lb = getvarvalue(master, getfuturecostvarid(cuts),1)
#             display(lb)
            ub = 0
            for (i,sub) in enumerate(subs)

                transferboundarystorage!(master, sub, states)

                solve!(sub)

                ub += getobjectivevalue(sub)
                cutparameters[i] = getcutparameters(sub, states)
            end
            ub /= numscen
            updatecuts!(master, cuts, cutparameters)
#             display(ub)
#             display(abs((lb-ub)/lb))
#             display(cuts.slopes)
        end

        return (master, subs, states, cuts)
    end
end

function pl_stochastic_init!(numcores, storagesystemobjects, shorts, masters_, subs_, states_, cuts_, storageinfo, numscen, lb, maxcuts, reltol, t)
    @sync @distributed for core in 1:(numcores-1)
        storagesystemobject = localpart(storagesystemobjects)
        short = localpart(shorts)
        masters = localpart(masters_)
        subs = localpart(subs_)
        states = localpart(states_)
        cuts = localpart(cuts_)

        localix = 0
        for range in localindices(storagesystemobjects)
            for ix in range
                localix += 1
                masterobjects, subobjects = storagesystemobject[localix]
                masters[localix], subs[localix], states[localix], cuts[localix] = stochastic_init(masterobjects, subobjects, short[localix], storageinfo, numscen, lb, maxcuts, reltol, t)
            end
        end
    end
end

function updatestochasticprices!(prob::Prob, prices::Vector{Dict}, scenario::Int)
    price = prices[scenario]
    for obj in prob.objects
        if obj isa ExogenBalance
            priceix = findfirst(x -> x == split(getinstancename(getid(obj)), "PowerBalance_")[2], price["names"])
            
            obj.price.param.profile.index = price["steprange"]
            obj.price.param.profile.values .= price["matrix"][:, priceix]
        end
    end
end

@everywhere function stochastic!(master, subs, states, cuts, startstates, medprices, shortprices, medendvaluesdict, short, numscen, reltol, t)
        
    # Init cutparameters
    cutparameters = Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}(undef, length(subs)) # preallocate for cutparameters from subproblems

    # Update master
    masterstorages = getstorages(master.objects)
    setstartstates!(master, masterstorages, startstates)
    short && updatestochasticprices!(master, shortprices, 1)
    !short && updatestochasticprices!(master, medprices, 1)

    update!(master, t)

    # Update subs
    for (i,sub) in enumerate(subs)
        
        substorages = getstorages(getobjects(sub))
        
        if short
            setendstates!(sub, substorages, startstates) # set end reservoir
            updatestochasticprices!(sub, shortprices, i) 
        else
            subendvaluesobj = sub.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), sub.objects)]
            subendvalues = [medendvaluesdict[getinstancename(getid(obj))] for obj in substorages]
            updateendvalues!(sub, subendvaluesobj, subendvalues)
            updatestochasticprices!(sub, medprices, i) 
        end

        update!(sub, t) # update parameters given problem start time of scenario
    end

    lb = cuts.lower_bound
    ub = 0

    count = 0
    cutreuse = true
    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)

        # TODO: Other way to fix this?
        if count == 0 # try to reuse cuts from last iteration
            try
                solve!(master)
            catch
                println("Retrying first iteration without cuts")
                clearcuts!(master, cuts)
                solve!(master)
                cutreuse = false
            end
        else
            try
                solve!(master)
            catch
                println("Restarting iterations without cuts")
                cutreuse && clearcuts!(master, cuts)
                solve!(master)
                cutreuse = false
            end
        end

        lb = getvarvalue(master, getfuturecostvarid(cuts),1)
        ub = 0
        
        for (i,sub) in enumerate(subs)

            transferboundarystorage!(master, sub, states)
            
            solve!(sub)

            ub += getobjectivevalue(sub)
            cutparameters[i] = getcutparameters(sub, states)
        end
        ub /= numscen

        count += 1
        count == 1 && clearcuts!(master, cuts) # reuse cuts in first iteration
        updatecuts!(master, cuts, cutparameters)
#             display(ub)
#             display(abs((lb-ub)/lb))
#             display(cuts.slopes)
    end
end

function pl_stochastic!(numcores, masters_, subs_, states_, cuts_, startstates, medprices, shortprices, medendvaluesdict, shorts, numscen, reltol, t, skipmed)
    @sync @distributed for core in 1:(numcores-1)
        masters = localpart(masters_)
        subs = localpart(subs_)
        states = localpart(states_)
        cuts = localpart(cuts_)
        short = localpart(shorts)

        localix = 0
        for range in localindices(shorts)
            for ix in range
                localix += 1
                if !(!short[localix] && (skipmed.value > 0))
                    stochastic!(masters[localix], subs[localix], states[localix], cuts[localix], startstates, medprices, shortprices, medendvaluesdict, short[localix], numscen, reltol, t)
                end
            end
        end
        if skipmed.value == 0
            GC.gc() # regulary garbage collection
        end
    end
end




