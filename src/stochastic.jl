# Make modelobjects for stochastic subsystems, group into subsystems
function makestochasticobjects(elements::Vector{DataElement}, problemduration::Millisecond, pdp::Millisecond, pdh::Millisecond, offset::Union{Offset,Nothing}, scenario::Int, prices::Vector{Dict}, short::Bool, master::Bool)
    # Add horizons to elements
    battery_horizon = SequentialHorizon(ceil(Int64, problemduration/pdp), pdp; offset) # TODO: Quickfix
    hydro_horizon = SequentialHorizon(ceil(Int64, problemduration/pdh), pdh; offset)
    power_horizon = SequentialHorizon(ceil(Int64, problemduration/pdp), pdp; offset)

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
        push!(elements, getelement(TIMEVECTOR_CONCEPT, "MutableInfiniteTimeVector", "ProfilePrice_" * area,
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
function makemastersubobjects!(inputs::Tuple{Vector{DataElement}, Millisecond, Millisecond, Millisecond, Millisecond, Millisecond, Vector{Tuple{Any, Any, Int64}}, Millisecond, Vector{Dict}, Bool}, mastersubobjects::Vector{Tuple{Vector, Vector{Vector}}}, shorts::Vector{Bool})
    (elements, totalduration, mpdp, mpdh, spdp, spdh, scenarios, phaseinoffset, prices, short) = inputs

    elements1 = copy(elements)
    removeelements!(elements1, short)

    masterobjects = makestochasticobjects(copy(elements1), phaseinoffset, mpdp, mpdh, nothing, 1, prices, short, true) # TODO: what price scenario price to use here? random? now 1, use of phasein of scenarios gives similar prices in the start of all scenarios?

    subscenarioobjects = []
    for (tnormal, tphasein, scenario) in scenarios
        offset = TimeDeltaOffset(MsTimeDelta(phaseinoffset))
        push!(subscenarioobjects, makestochasticobjects(copy(elements1), totalduration - phaseinoffset, spdp, spdh, offset, scenario, prices, short, false))
    end

    for (i, masterobject) in enumerate(masterobjects)
        subobjects = [v[i] for v in subscenarioobjects]
        push!(mastersubobjects, (masterobject, subobjects))
        push!(shorts, short)
    end
end

# Aggregate modelobjects and remove modelobjects not relevant for subsystems
function removeelements!(elements::Vector{DataElement}, short) # TODO: Replace with user settings
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

# Initialize cuts
function initialize_cuts!(modelobjects::Vector, cutobjects::Vector, maxcuts::Int, lb::Float64, numscen::Int)
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
function transferboundarystates!(master::Prob, sub::Prob, states::Dict{StateVariableInfo, Float64})
    states = getoutgoingstates!(master, states)
    # display(states)
    setingoingstates!(sub, states)
end

# Initialize stochastic subsystem problems and solve for first time step
function stochastic_init(probmethods::Vector, masterobjects::Vector, subobjects::Vector{Vector}, short::Bool, storageinfo::Tuple{Float64,Float64,Vector{Dict}}, lb::Float64, maxcuts::Int, reltol::Float64, scenarios::Vector{Tuple{Any, Any, Int64}})
    shortstartstorage, medstartstorage, medendvaluesdicts = storageinfo
    if short
        startstorage = shortstartstorage
    else
        startstorage = medstartstorage
    end

    cutobjects = getcutobjects(masterobjects)
    cuts = initialize_cuts!(masterobjects, cutobjects, maxcuts, lb, length(scenarios));
    states = getstates(cutobjects) # state variables in master and subs for boundary reservoirs

    # master = JuMP_Prob(masterobjects, Model(HiGHS.Optimizer))
    # subs = [JuMP_Prob(subobject, Model(HiGHS.Optimizer)) for subobject in subobjects] # initialize subproblems
    master = buildprob(probmethods[1], masterobjects)
    subs = [buildprob(probmethods[2], subobject) for subobject in subobjects] # initialize subproblems

    # Init cutparameters
    cutparameters = Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}(undef, length(subs)) # preallocate for cutparameters from subproblems

    # Update master
    (tnormalmaster, tphaseinmaster, scenario) = scenarios[1] # tphasein is the same for all scenarios before phasein
    update!(master, tphaseinmaster) 
    setstartstoragepercentage!(master, getstorages(getobjects(master)), tphaseinmaster, startstorage)

    # Update subs
    for (i,sub) in enumerate(subs)
        (tnormalsub, tphaseinsub, scenario) = scenarios[i]
        update!(sub, tphaseinsub) # update parameters given problem start time of scenario

        storages = getstorages(getobjects(sub))
        if short
            setendstoragepercentage!(sub, storages, tphaseinsub, startstorage) # set end reservoir
        else
            subendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
            subendvaluesobj = EndValues(subendvaluesid, storages)
            push!(sub.objects, subendvaluesobj)
            subendvalues = [medendvaluesdicts[scenario][getinstancename(getid(obj))] for obj in storages]
            updateendvalues!(sub, subendvaluesobj, subendvalues)
        end
    end

    ub = 0
    cutreuse = false
    iterate_convergence!(master, subs, cuts, cutparameters, states, cutreuse, lb, ub, reltol)

    # Move solution from HiGHS instance to HiGHS_Prob. Has to be done on parallel processor because 
    # HiGHS API cannot be used when master is moved to local process.
    # We use the dual values of the master problem to calculate the headloss costs
    # Possible TODO
    if master isa Union{HiGHS_Prob, CPLEX_Prob}
        setconduals!(master)
        master.iscondualsupdated = true
    end

    return (master, subs, states, cuts)
end

# Iterate until convergence between master and subproblems
function iterate_convergence!(master::Prob, subs::Vector, cuts::SimpleSingleCuts, cutparameters::Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}, states::Dict{StateVariableInfo, Float64}, cutreuse::Bool, lb::Float64, ub::Int, reltol::Float64)
    count = 0

    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)

        if (count == 0) && (master isa HiGHS_Prob) # implement this for other solvers
            master.warmstart = false
        elseif (count == 0) && (master isa CPLEX_Prob)
            setparam!(master, "CPXPARAM_Advance", 0)
        end
        if cutreuse # try to reuse cuts from last time step
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
        if (count == 0) && (master isa HiGHS_Prob) 
            master.warmstart = true
        elseif (count == 0) && (master isa CPLEX_Prob)
            setparam!(master, "CPXPARAM_Advance", 1) # TODO: What if setting is 2?
        end

        lb = getvarvalue(master, getfuturecostvarid(cuts),1)
        ub = 0
        
        for (i,sub) in enumerate(subs)

            transferboundarystates!(master, sub, states)
            
            solve!(sub)

            ub += getobjectivevalue(sub)*cuts.probabilities[i]
            cutparameters[i] = getcutparameters(sub, states)
        end

        count += 1
        (count == 1 && cutreuse) && clearcuts!(master, cuts) # reuse cuts in first iteration
        updatecuts!(master, cuts, cutparameters)
#             display(ub)
#             display(abs((lb-ub)/lb))
#             display(cuts.slopes)
    end
end

# Initialize stochastic subsystem problems in parallel
function pl_stochastic_init!(probmethods::Vector, numcores::Int, storagesystemobjects::DArray, shorts::DArray, masters_::DArray, subs_::DArray, states_::DArray, cuts_::DArray, storageinfo::Tuple{Float64, Float64, Vector{Dict}}, lb::Float64, maxcuts::Int, reltol::Float64, scenarios::Vector{Tuple{Any, Any, Int}})
    @sync @distributed for core in 1:max(numcores-1,1)
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
                masters[localix], subs[localix], states[localix], cuts[localix] = stochastic_init(probmethods, masterobjects, subobjects, short[localix], storageinfo, lb, maxcuts, reltol, scenarios)
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
function stochastic!(master::Prob, subs::Vector, states::Dict{StateVariableInfo, Float64}, cuts::SimpleSingleCuts, startstates::Dict, medprices::Vector{Dict}, shortprices::Vector{Dict}, medendvaluesdicts::Vector{Dict}, short::Bool, reltol::Float64, scenarios::Vector{Tuple{Any, Any, Int64}})
        
    # Init cutparameters
    cutparameters = Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}(undef, length(subs)) # preallocate for cutparameters from subproblems

    # Update master
    masterstorages = getstorages(master.objects)
    setstartstates!(master, masterstorages, startstates)
    short && updatestochasticprices!(master, shortprices, 1) # TODO: what price scenario price to use here? random? now 1, use of phasein of scenarios gives similar prices in the start of all scenarios?
    !short && updatestochasticprices!(master, medprices, 1)

    (tnormalmaster, tphaseinmaster) = scenarios[1] # tphasein is the same for all scenarios before phasein
    update!(master, tphaseinmaster)

    # Update subs
    for (i,sub) in enumerate(subs)
        (tnormalsub, tphaseinsub, scenario) = scenarios[i]
        substorages = getstorages(getobjects(sub))
        
        if short
            setendstates!(sub, substorages, startstates) # set end reservoir
            updatestochasticprices!(sub, shortprices, scenario)
        else
            subendvaluesobj = sub.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), sub.objects)]
            subendvalues = [medendvaluesdicts[scenario][getinstancename(getid(obj))] for obj in substorages]
            updateendvalues!(sub, subendvaluesobj, subendvalues)
            updatestochasticprices!(sub, medprices, scenario) 
        end
        
        update!(sub, tnormalsub) # update parameters given problem start time of scenario
    end

    lb = cuts.lower_bound
    ub = 0
    cutreuse = true
    iterate_convergence!(master, subs, cuts, cutparameters, states, cutreuse, lb, ub, reltol)

    # Move solution from HiGHS instance to HiGHS_Prob. Has to be done on parallel processor because 
    # HiGHS API cannot be used when master is moved to local process.
    # We use the dual values of the master problem to calculate the headloss costs
    # Possible TODO
    if master isa Union{HiGHS_Prob, CPLEX_Prob}
        setconduals!(master)
        master.iscondualsupdated = true
    end
end

# Run stochastic subsystem problems in parallel
function pl_stochastic!(numcores::Int, masters_::DArray, subs_::DArray, states_::DArray, cuts_::DArray, startstates::Dict{String, Float64}, medprices::Vector{Dict}, shortprices::Vector{Dict}, medendvaluesdicts::Vector{Dict}, shorts::DArray, reltol::Float64, scenarios::Vector{Tuple{Any, Any, Int64}}, skipmed::Millisecond)
    @sync @distributed for core in 1:max(numcores-1,1)
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
                    stochastic!(masters[localix], subs[localix], states[localix], cuts[localix], startstates, medprices, shortprices, medendvaluesdicts, short[localix], reltol, scenarios)
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
function distribute_subsystems(ustoragesystemobjects::Vector{Tuple{Vector, Vector{Vector}}}, ushorts::Vector)
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




