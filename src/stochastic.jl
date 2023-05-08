# Make modelobjects for stochastic subsystems, group into subsystems
function makestochasticobjects(elements, days::Int, offset::Union{Offset,Nothing}, scenario::Int, prices::Vector{Dict}, short::Bool, master::Bool)
    
    # Add horizons to elements
    if short
        battery_horizon = SequentialHorizon(days*12, Hour(2); offset) # replace with user settings
        hydro_horizon = SequentialHorizon(days*12, Hour(2); offset)
        power_horizon = SequentialHorizon(days*12, Hour(2); offset)
    else
        if days < 7
            battery_horizon = SequentialHorizon(days, Hour(24); offset)
            hydro_horizon = SequentialHorizon(days, Hour(24); offset)
            power_horizon = SequentialHorizon(days, Hour(24); offset)
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
    
    # Add prices to elements
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

    # Make modelobjets from elements and group into subsystems
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

    # Sort storagesystems
    for storagesystem in storagesystems
        sort!(storagesystem, by = x -> getinstancename(getid(x)))
    end
    sort!(storagesystems, by = x -> getinstancename(getid(first(x))))

    return storagesystems
end

# Make master and subproblem objects for each subsystem
function makemastersubobjects!(inputs, mastersubobjects, shorts)
    (elements, totaldays, numscen, horizonstart, phaseinoffsetdays, prices, short) = inputs

    elements1 = copy(elements)
    removeelements!(elements1, short)

    masterobjects = makestochasticobjects(copy(elements1), phaseinoffsetdays, nothing, 1, prices, short, true) # TODO: what scenario price to use here? random? now 1, use of phasein of scenarios gives similar prices in the start of all scenarios?

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

# Aggregate modelobjects and remove modelobjects not relevant for subsystems
function removeelements!(elements, short) # TODO: Replace with user settings
    # Aggregate areas
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

        if element.typename == "HydroRampingWithout"
            push!(delix,i)
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
    cutname = getinstancename(getid(modelobjects[1]))
    cutid = Id(BOUNDARYCONDITION_CONCEPT,"StorageCuts" * cutname)
    
    # Probability of each subproblem / second stage scenario
    probabilities = [1/numscen for i in 1:numscen]
    
    # Make cut modelobject
    cuts = SimpleSingleCuts(cutid, cutobjects, probabilities, maxcuts, lb)
    push!(modelobjects, cuts)
    return cuts
end

# Transfer master problem end states to subproblem start states (currently only storage states)
function transferboundarystates!(master, sub, states)
    states = getoutgoingstates!(master, states)
    # display(states)
    setingoingstates!(sub, states)
end

# Initialize stochastic subsystem problems and solve for first time step
function stochastic_init(masterobjects, subobjects, short, storageinfo, numscen, lb, maxcuts, reltol, t)
    @time begin
        shortstartstorage, medstartstorage, medendvaluesdicts = storageinfo
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
        setstartstoragepercentage!(master, getstorages(getobjects(master)), t, startstorage)

        # Update subs
        for (i,sub) in enumerate(subs)
            update!(sub, t) # update parameters given problem start time of scenario

            storages = getstorages(getobjects(sub))
            if short
                setendstoragepercentage!(sub, storages, t, startstorage) # set end reservoir
            else
                subendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
                subendvaluesobj = EndValues(subendvaluesid, storages)
                push!(sub.objects, subendvaluesobj)
                subendvalues = [medendvaluesdicts[i][getinstancename(getid(obj))] for obj in storages]
                updateendvalues!(sub, subendvaluesobj, subendvalues)
            end
        end

        ub = 0
        cutreuse = false
        iterate_convergence!(master, subs, cuts, cutparameters, states, numscen, cutreuse, lb, ub, reltol)

        return (master, subs, states, cuts)
    end
end

# Iterate until convergence between master and subproblems
function iterate_convergence!(master, subs, cuts, cutparameters, states, numscen, cutreuse, lb, ub, reltol)
    count = 0

    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)

        if cutreuse # try to reuse cuts from last iteration
            try
                solve!(master)
            catch
                count == 0 && println("Retrying first iteration without cuts from last time step")
                count > 0 && println("Restarting iterations without cuts from last time step")
                clearcuts!(master, cuts)
                solve!(master)
                cutreuse = false
            end
        else
            solve!(master)
        end

        lb = getvarvalue(master, getfuturecostvarid(cuts),1)
        ub = 0
        
        for (i,sub) in enumerate(subs)

            transferboundarystates!(master, sub, states)
            
            solve!(sub)

            ub += getobjectivevalue(sub)
            cutparameters[i] = getcutparameters(sub, states)
        end
        ub /= numscen

        count += 1
        (count == 1 && cutreuse) && clearcuts!(master, cuts) # reuse cuts in first iteration
        updatecuts!(master, cuts, cutparameters)
#             display(ub)
#             display(abs((lb-ub)/lb))
#             display(cuts.slopes)
    end
end

# Initialize stochastic subsystem problems in parallel
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

# Update prices in exogen prices
function updatestochasticprices!(prob::Prob, prices::Vector{Dict}, scenario::Int)
    price = prices[scenario]
    for obj in prob.objects
        if obj isa ExogenBalance
            priceix = findfirst(x -> x == split(getinstancename(getid(obj)), "PowerBalance_")[2], price["names"])
            
            obj.price.param.profile.index = price["steprange"] # TODO: API function in TuLiPa
            obj.price.param.profile.values .= price["matrix"][:, priceix]
        end
    end
end

# Run stochastic subsystem problem
@everywhere function stochastic!(master, subs, states, cuts, startstates, medprices, shortprices, medendvaluesdicts, short, numscen, reltol, t)
        
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
            subendvalues = [medendvaluesdicts[i][getinstancename(getid(obj))] for obj in substorages]
            updateendvalues!(sub, subendvaluesobj, subendvalues)
            updatestochasticprices!(sub, medprices, i) 
        end

        update!(sub, t) # update parameters given problem start time of scenario
    end

    lb = cuts.lower_bound
    ub = 0
    cutreuse = true
    iterate_convergence!(master, subs, cuts, cutparameters, states, numscen, cutreuse, lb, ub, reltol)
end

# Run stochastic subsystem problems in parallel
function pl_stochastic!(numcores, masters_, subs_, states_, cuts_, startstates, medprices, shortprices, medendvaluesdicts, shorts, numscen, reltol, t, skipmed)
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
                    stochastic!(masters[localix], subs[localix], states[localix], cuts[localix], startstates, medprices, shortprices, medendvaluesdicts, short[localix], numscen, reltol, t)
                end
            end
        end
        if skipmed.value == 0
            GC.gc() # regulary garbage collection, else there is always one processor core for each time step that the others have to wait for
        end
    end
end

# Distribute subsystems after number of objects in subsystem, equal amount of subsystems on each core
# Consist of several elements:
# - Testsort to see how an array of this size would be distributed by default by DistributedArrays
# - Make rules for how to distribute subsystems 
#   - Example: 6 elements on 3 cores, assume it is already sorted by most modelobjects in subsystems
#   - Should distribute: 
#       -  Core 1: Index 1 and 6
#       -  Core 2: Index 2 and 5
#       -  Core 3: Index 3 and 4
# TODO: Replace with low level Distributed functions that are more flexible than DistributedArrays
function distribute_subsystems(ustoragesystemobjects::Vector, ushorts::Vector)
    objectcount = [length(first(objs)) for objs in ustoragesystemobjects]
    sortedindex = sortperm(objectcount, rev=true)
    testsort = distribute(collect(1:length(ustoragesystemobjects)))
    storagesystemobjects = Vector(undef,length(objectcount))
    shorts = Vector(undef,length(ushorts))
    direction = 1
    rangeix = 1
    localix = 1
    for i in 1:length(objectcount)
        (range,) = testsort.indices[rangeix]
        if localix > length(range) # possible for last localix
            rangeix = 1
            direction = 1
            (range,) = testsort.indices[rangeix]
        end

        storagesystemobjects[collect(range)[localix]] = ustoragesystemobjects[sortedindex][i]
        shorts[collect(range)[localix]] = ushorts[sortedindex][i]
        rangeix += direction

        if rangeix > length(testsort.indices)
            direction = -1
            rangeix += direction
            localix += 1
        elseif rangeix == 0
            direction = 1
            rangeix += direction
            localix += 1
        end
    end
    storagesystemobjects = distribute(storagesystemobjects)
    shorts = distribute(shorts, storagesystemobjects)
    return storagesystemobjects, shorts
end




