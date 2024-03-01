# Make modelobjects for stochastic subsystems, group into subsystems
function makestochasticobjects(elements::Vector{DataElement}, problemduration::Millisecond, pdp::Millisecond, pdh::Millisecond, offset::Union{Offset,Nothing}, scenario::Int, prices::Union{Vector{Dict},Nothing}, short::Bool, master::Bool, validate::Bool, settings::Dict)
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
    if prices != nothing
        price = prices[scenario]

        push!(elements, DataElement(TIMEINDEX_CONCEPT, "RangeTimeIndex", "ShortTermTimeIndex", price["steprange"]))
        
        for (i, area) in enumerate(price["names"])
            push!(elements, getelement(TIMEVALUES_CONCEPT, "VectorTimeValues", "ValuesPrice_" * area,
                ("Vector", price["matrix"][:, i])))
            push!(elements, getelement(TIMEVECTOR_CONCEPT, "MutableInfiniteTimeVector", "ProfilePrice_" * area,
                (TIMEINDEX_CONCEPT, "ShortTermTimeIndex"), (TIMEVALUES_CONCEPT, "ValuesPrice_" * area)))
            push!(elements, getelement(PARAM_CONCEPT, "MeanSeriesIgnorePhaseinParam", "SeriesPrice_" * area,
                ("Level", 1.0),
                ("Profile", "ProfilePrice_" * area)))
            push!(elements, getelement(PARAM_CONCEPT, "StatefulParam", "Price_" * area,
                (PARAM_CONCEPT, "SeriesPrice_" * area)))
            push!(elements, getelement(BALANCE_CONCEPT, "ExogenBalance", "PowerBalance_" * area, 
                (COMMODITY_CONCEPT, "Power"),
                (PRICE_CONCEPT, "Price_" * area)))
        end

        # Make modelobjets from elements and group into subsystems
        modelobjects = getmodelobjects(elements, validate=validate)
        if short
            if master # removes spills from upper and lower storages in PHS, to avoid emptying reservoirs in master problem
                for id in keys(modelobjects)
                    instance = getinstancename(id)
                    if occursin("Spill", instance) && occursin("_PHS_", instance)
                        delete!(modelobjects, id)
                    end
                end
            end
            storagesystems = getstoragesystems_full!(getshorttermstoragesystems(getstoragesystems(modelobjects), Hour(settings["problems"]["shorttermstoragecutoff_hours"])))
        else
            storagesystems = getstoragesystems_full!(getlongtermstoragesystems(getstoragesystems(modelobjects), Hour(settings["problems"]["shorttermstoragecutoff_hours"])))
        end

        # Sort storagesystems
        for storagesystem in storagesystems
            sort!(storagesystem, by = x -> getinstancename(getid(x)))
        end
        sort!(storagesystems, by = x -> getinstancename(getid(first(x))))

        return storagesystems, modelobjects
    else
        modelobjects = getmodelobjects(elements, validate=validate)
        storagesystems = [collect(values(modelobjects))]
        return storagesystems, modelobjects
    end
end

# Make master and subproblem objects for each subsystem
function makemastersubobjects!(inputs::Tuple{Vector{DataElement}, Millisecond, Millisecond, Millisecond, Millisecond, Millisecond, Vector{Tuple{Any, Any, Int64}}, Millisecond, Union{Vector{Dict}, Nothing}, Bool}, mastersubobjects::Vector{Tuple{Vector, Vector{Vector}}}, shorts::Vector{Bool}, settings::Dict)
    (elements, totalduration, mpdp, mpdh, spdp, spdh, scenarios, phaseinoffset, prices, short) = inputs

    # Only validate dataelements once
    masterobjects, modelobjects = makestochasticobjects(copy(elements), phaseinoffset, mpdp, mpdh, nothing, 1, prices, short, true, false, settings) # TODO: what price scenario price to use here? random? now 1, use of phasein of scenarios gives similar prices in the start of all scenarios?

    subscenarioobjects = []
    for (tnormal, tphasein, scenario) in scenarios
        offset = TimeDeltaOffset(MsTimeDelta(phaseinoffset))
        subobject, dummyobjects = makestochasticobjects(copy(elements), totalduration - phaseinoffset, spdp, spdh, offset, scenario, prices, short, false, false, settings)
        push!(subscenarioobjects, subobject)
    end

    for (i, masterobject) in enumerate(masterobjects)
        subobjects = [v[i] for v in subscenarioobjects]
        push!(mastersubobjects, (masterobject, subobjects))
        push!(shorts, short)
    end
    return modelobjects
end

# Aggregate modelobjects and remove modelobjects not relevant for subsystems
function removeelements!(elements::Vector{DataElement}; aggzone::Dict=Dict(), rm_basebalances::Bool=true) # TODO: Replace with more user settings
    aggzonecopl = Dict()
    for (k,v) in aggzone
        for vv in v
            aggzonecopl["PowerBalance_" * vv] = "PowerBalance_" * k
        end
    end
    
    delix = []
    powerbasebalances = []
    for (i,element) in enumerate(elements)
        if rm_basebalances
            # BaseBalances
            if element.typename == "BaseBalance"
                if element.value["Commodity"] == "Power"
                    push!(delix,i)
                    push!(powerbasebalances, element.instancename)
                end
            end

            # Power balance RHSTerms
            if element.conceptname == "RHSTerm"
                if element.value["Balance"] in powerbasebalances
                    push!(delix,i)
                end
            end
            
            # Residualhints
            if element.typename == "Residualhint"
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
    for deli in sort(delix; rev=true)
        popat!(elements, deli)
    end
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
function initialize_cuts!(modelobjects::Vector, cutobjects::Vector, maxcuts::Int, lb::Float64, numscen::Int, probabilities::Vector)
    # Make a cutid
    cutname = getinstancename(getid(modelobjects[1]))
    cutid = Id(BOUNDARYCONDITION_CONCEPT,"StorageCuts" * cutname)
    
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
function stochastic_init!(probmethods::Vector, masterobjects::Vector, subobjects::Vector{Vector}, storagevalues::Array, short::Bool, storageinfo::Tuple{Dict{String, Float64},Vector{Dict}}, lb::Float64, maxcuts::Int, reltol::Float64, tnormal::ProbTime, scenarios::ScenarioModellingMethod, settings::Dict)
    startstates, medendvaluesdicts = storageinfo

    cutobjects = getcutobjects(masterobjects)
    cuts = initialize_cuts!(masterobjects, cutobjects, maxcuts, lb, length(scenarios.scentimes), scenarios.weights);
    states = getstates(cutobjects) # state variables in master and subs for boundary reservoirs

    master = buildprob(probmethods[1], masterobjects)
    subs = [buildprob(probmethods[2], subobject) for subobject in subobjects] # initialize subproblems

    # Init cutparameters
    cutparameters = Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}(undef, length(subs)) # preallocate for cutparameters from subproblems

    getstatedependentprod(settings["problems"]["stochastic"]["master"]) && statedependentprod!(master, startstates, init=true)
    getstatedependentpump(settings["problems"]["stochastic"]["master"]) && statedependentpump!(master, startstates)

    # Update master
    update!(master, tnormal)
    setstartstates!(master, getstorages(getobjects(master)), startstates)

    # Find price
    medendvaluesdicts == Dict[] && global exogenprice = findfirstprice(master.objects)

    # Update subs
    for (i,sub) in enumerate(subs)
        (tnormalsub, tphaseinsub, scenario) = scenarios.scentimes[i]
        scaleinflow!(scenarios, i, getobjects(sub)) # scale inflow to average of represented scenarios
        update!(sub, tphaseinsub) # update parameters given problem start time of scenario

        storages = getstorages(getobjects(sub))
        if short
            setendstates!(sub, storages, startstates) # set end reservoir
        else
            subendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
            subendvaluesobj = EndValues(subendvaluesid, storages)
            push!(getobjects(sub), subendvaluesobj)
            if medendvaluesdicts != Dict[]
                subendvalues = [medendvaluesdicts[scenario][getinstancename(getid(obj))] for obj in storages]
                updateendvalues!(sub, subendvaluesobj, subendvalues)
            else
                # Alternative end reservoir conditions:
                # - Stop equals start reservoir, and long enough horizon (can crash if many start reservoirs full or empty)
                # - End value = 0, and long enough horizon
                # - List of endvalues from other calculation, maybe SDP
                # - End value equals average price in this price area for different scenarios (could also adjust based on reservoir flexibility)

                # setendstoragepercentage!(sub, storages, scentphasein, startstorage) # set end reservoir

                # Set monthly price at end of horizon as end value
                scenprice = getparamvalue(exogenprice, tnormalsub + getduration(gethorizon(storages[1])), MsTimeDelta(Week(4))) 
                subendvalues = scenprice .* [getbalance(obj).metadata[GLOBALENEQKEY] for obj in storages] # end value also depends on the energy equivalent of the reservoir, the energy equivalent is stored as metadata
                updateendvalues!(sub, subendvaluesobj, subendvalues)
            end
        end
    end

    ub = 0
    cutreuse = false
    count, mastertime, subtime = iterate_convergence!(master, subs, cuts, cutparameters, states, cutreuse, lb, ub, reltol)

    # Init storagevalues results # TODO move to function
    if settings["results"]["storagevalues"]
        for i in 1:length(subs)
            (constant, slopes) = cutparameters[i]
            for (j, (state, value)) in enumerate(slopes)
                minslope = 0
                maxslope = -1e9
                for k in 1:getnumcuts(cuts)
                    val = cuts[i][k][state]
                    minslope = min(minslope, val)
                    maxslope = max(maxslope, val)
                end
                balance = getbalance(cuts.objects[j])
                if haskey(balance.metadata, GLOBALENEQKEY)
                    minslope = minslope / balance.metadata[GLOBALENEQKEY]
                    maxslope = maxslope / balance.metadata[GLOBALENEQKEY]
                end
                storagevalues[1, (i-1)*2+1, j] = minslope
                storagevalues[1, (i-1)*2+2, j] = maxslope
            end
        end

        for (j, obj) in enumerate(cuts.objects) # operative water values
            balance = getbalance(obj)
            storagevalues[1, length(subs)*2 + 1, j] = getcondual(master, getid(balance), getnumperiods(gethorizon(balance)))
            if haskey(balance.metadata, GLOBALENEQKEY)
                storagevalues[1, length(subs)*2 + 1, j] = storagevalues[1, length(subs)*2 + 1, j] / balance.metadata[GLOBALENEQKEY]
            end
        end
    end

    # Final master run with headloss adjusted storagevalues 
    if getheadlosscost(settings["problems"]["stochastic"]["master"])
        updateheadlosscosts!(ReservoirCurveSlopeMethod(), master, [master], tnormal)
        solve!(master)
        resetheadlosscosts!(master)

        if storagevalues != [0]
            for (j, obj) in enumerate(cuts.objects) # operative water values after headlosscost
                balance = getbalance(obj)
                storagevalues[1, length(subs)*2 + 2, j] = getcondual(master, getid(balance), getnumperiods(gethorizon(balance)))
                if haskey(balance.metadata, GLOBALENEQKEY)
                    storagevalues[1, length(subs)*2 + 2, j] = storagevalues[1, length(subs)*2 + 2, j] / balance.metadata[GLOBALENEQKEY]
                end
            end
        end
    end

    # Move solution from HiGHS instance to HiGHS_Prob. Has to be done on parallel processor because 
    # HiGHS API cannot be used when master is moved to local process.
    # We use the dual values of the master problem to calculate the headloss costs
    # Possible TODO
    if master isa HiGHS_Prob || is_CPLEX_Prob(master)
        setconduals!(master)
        master.iscondualsupdated = true
    end

    return (master, subs, states, cuts)
end

# Iterate until convergence between master and subproblems
function iterate_convergence!(master::Prob, subs::Vector, cuts::SimpleSingleCuts, cutparameters::Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}, states::Dict{StateVariableInfo, Float64}, cutreuse::Bool, lb::Float64, ub::Int, reltol::Float64)
    count = 0
    mastertime = 0
    subtime = 0

    while !((abs((ub-lb)/ub) < reltol) || abs(ub-lb) < 1)

        count == 0 && setwarmstart!(master, false)

        mastertime += @elapsed begin
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
        end
        
        count == 0 && setwarmstart!(master, true)

        lb = getvarvalue(master, getfuturecostvarid(cuts),1)
        ub = 0
        
        for (i,sub) in enumerate(subs)

            transferboundarystates!(master, sub, states)
            
            subtime += @elapsed solve!(sub)

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
    
    return count, mastertime, subtime
end

# Initialize stochastic subsystem problems in parallel
function pl_stochastic_init!(probmethods::Vector, numcores::Int, storagesystemobjects::DArray, shorts::DArray, masters_::DArray, subs_::DArray, states_::DArray, cuts_::DArray, storagevalues_::DArray, storageinfo::Tuple{Dict{String, Float64}, Vector{Dict}}, lb::Float64, maxcuts::Int, reltol::Float64, tnormal::ProbTime, scenarios::ScenarioModellingMethod, settings::Dict)
    @sync @distributed for core in 1:max(numcores-1,1)
        storagesystemobject = localpart(storagesystemobjects)
        short = localpart(shorts)
        masters = localpart(masters_)
        subs = localpart(subs_)
        states = localpart(states_)
        cuts = localpart(cuts_)
        storagevalues = localpart(storagevalues_)

        localix = 0
        for range in localindices(storagesystemobjects)
            for ix in range
                localix += 1
                masterobjects, subobjects = storagesystemobject[localix]
                masters[localix], subs[localix], states[localix], cuts[localix] = stochastic_init!(probmethods, masterobjects, subobjects, storagevalues[localix], short[localix], storageinfo, lb, maxcuts, reltol, tnormal, scenarios, settings)
            end
        end
    end
end

# Update prices in exogen prices
function updatestochasticprices!(prob::Prob, prices::Vector{Dict}, scenario::Int)
    price = prices[scenario]
    for obj in getobjects(prob)
        if obj isa ExogenBalance
            priceix = findfirst(x -> x == split(getinstancename(getid(obj)), "PowerBalance_")[2], price["names"])
            
            obj.price.param.param.profile.index = price["steprange"] # TODO: API function in TuLiPa
            obj.price.param.param.profile.values .= price["matrix"][:, priceix]
        end
    end
end

# Run stochastic subsystem problem
function stochastic!(master::Prob, subs::Vector, states::Dict{StateVariableInfo, Float64}, cuts::SimpleSingleCuts, storagevalues::Array, startstates::Dict, medprices::Union{Vector{Dict}, Nothing}, shortprices::Union{Vector{Dict}, Nothing}, medendvaluesdicts::Vector{Dict}, short::Bool, reltol::Float64, tnormal::ProbTime, scenarios::ScenarioModellingMethod, stochastictimes::Matrix{Float64}, stepnr::Int, settings::Dict)
    stochastictimes[stepnr-1, 9] = @elapsed begin   
        # Init cutparameters
        cutparameters = Vector{Tuple{Float64, Dict{StateVariableInfo, Float64}}}(undef, length(subs)) # preallocate for cutparameters from subproblems

        # Update master
        masterstorages = getstorages(getobjects(master))
        setstartstates!(master, masterstorages, startstates)

        if medendvaluesdicts != Dict[]
            short && updatestochasticprices!(master, shortprices, 1) # TODO: what price scenario price to use here? random? now 1, use of phasein of scenarios gives similar prices in the start of all scenarios?
            !short && updatestochasticprices!(master, medprices, 1)
        else
            global exogenprice = findfirstprice(master.objects)
        end

        getstatedependentprod(settings["problems"]["stochastic"]["master"]) && statedependentprod!(master, startstates)
        getstatedependentpump(settings["problems"]["stochastic"]["master"]) && statedependentpump!(master, startstates)

        stochastictimes[stepnr-1, 1] = @elapsed update!(master, tnormal)

        # Update subs
        for (i,sub) in enumerate(subs)
            (tnormalsub, tphaseinsub, scenario) = scenarios.scentimes[i]
            substorages = getstorages(getobjects(sub))
            
            if short
                setendstates!(sub, substorages, startstates) # set end reservoir
                medendvaluesdicts != Dict[] && updatestochasticprices!(sub, shortprices, scenario)
            else
                if medendvaluesdicts != Dict[]
                    subendvaluesobj = getobjects(sub)[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), getobjects(sub))]
                    subendvalues = [medendvaluesdicts[scenario][getinstancename(getid(obj))] for obj in substorages]
                    updateendvalues!(sub, subendvaluesobj, subendvalues)
                    medendvaluesdicts != Dict[] && updatestochasticprices!(sub, medprices, scenario) 
                else
                    # Set monthly price at end of horizon as end value
                    scenprice = getparamvalue(exogenprice, tnormalsub + getduration(gethorizon(substorages[1])), MsTimeDelta(Week(4))) 
                    subendvaluesobj = getobjects(sub)[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), sub.objects)]
                    subendvalues = scenprice .* [getbalance(obj).metadata[GLOBALENEQKEY] for obj in substorages] # end value also depends on the energy equivalent of the reservoir, the energy equivalent is stored as metadata
                    updateendvalues!(sub, subendvaluesobj, subendvalues)
                end
            end
            
            scaleinflow!(scenarios, i, getobjects(sub)) # scale inflow to average of represented scenarios
            stochastictimes[stepnr-1, 2] += @elapsed update!(sub, tnormalsub) # update parameters given problem start time of scenario
        end

        cuts.probabilities = scenarios.weights
        lb = cuts.lower_bound
        ub = 0
        cutreuse = true
        stochastictimes[stepnr-1, 3] = @elapsed begin
            (count, mastertime, subtime) = iterate_convergence!(master, subs, cuts, cutparameters, states, cutreuse, lb, ub, reltol)
        end
        stochastictimes[stepnr-1, 4] = count
        stochastictimes[stepnr-1, 5] = mastertime
        stochastictimes[stepnr-1, 6] = subtime

        firstwwtime = @elapsed begin
            # Init storagevalues results
            wwcount = 0
            if storagevalues != [0]
                for i in 1:length(subs)
                    (constant, slopes) = cutparameters[i]
                    for (j, (state, value)) in enumerate(slopes[1])
                        minslope = 0
                        maxslope = -1e9
                        for k in 1:getnumcuts(cuts)
                            val = cuts[i][k][state]
                            minslope = min(minslope, val)
                            maxslope = max(maxslope, val)
                        end
                        balance = getbalance(cuts.objects[j])
                        if haskey(balance.metadata, GLOBALENEQKEY)
                            minslope = minslope / balance.metadata[GLOBALENEQKEY]
                            maxslope = maxslope / balance.metadata[GLOBALENEQKEY]
                        end
                        storagevalues[stepnr, (i-1)*2+1, j] = minslope
                        storagevalues[stepnr, (i-1)*2+2, j] = maxslope
                    end
                end
        
                for (j, obj) in enumerate(cuts.objects) # operative water values
                    balance = getbalance(obj)
                    storagevalues[stepnr, length(subs)*2 + 1, j] = getcondual(master, getid(balance), getnumperiods(gethorizon(balance)))
                    if haskey(balance.metadata, GLOBALENEQKEY)
                        storagevalues[stepnr, length(subs)*2 + 1, j] = storagevalues[stepnr, length(subs)*2 + 1, j] / balance.metadata[GLOBALENEQKEY]
                    end
                end
            end
        end

        # Final master run with headloss adjusted storagevalues 
        stochastictimes[stepnr-1, 7] = @elapsed begin
            if getheadlosscost(settings["problems"]["stochastic"]["master"])
                updateheadlosscosts!(ReservoirCurveSlopeMethod(), master, [master], tnormal)
                solve!(master)
                resetheadlosscosts!(master)
            end
        end

        secondwwtime = @elapsed begin
            if getheadlosscost(settings["problems"]["stochastic"]["master"])
                if storagevalues != [0]
                    for (j, obj) in enumerate(cuts.objects) # operative water values after headlosscost
                        balance = getbalance(obj)
                        storagevalues[stepnr, length(subs)*2 + 2, j] = getcondual(master, getid(balance), getnumperiods(gethorizon(balance)))
                        if haskey(balance.metadata, GLOBALENEQKEY)
                            storagevalues[stepnr, length(subs)*2 + 2, j] = storagevalues[stepnr, length(subs)*2 + 2, j] / balance.metadata[GLOBALENEQKEY]
                        end
                    end
                end
            end
        end
        stochastictimes[stepnr-1, 8] = firstwwtime + secondwwtime

        # Move solution from HiGHS instance to HiGHS_Prob. Has to be done on parallel processor because 
        # HiGHS API cannot be used when master is moved to local process.
        # We use the dual values of the master problem to calculate the headloss costs
        # Possible TODO
        if master isa HiGHS_Prob || is_CPLEX_Prob(master)
            setconduals!(master)
            master.iscondualsupdated = true
        end
    end
end

function pl_stochastic!(numcores::Int, masters_::DArray, subs_::DArray, states_::DArray, cuts_::DArray, storagevalues_::DArray, startstates::Dict{String, Float64}, medprices::Union{Vector{Dict}, Nothing}, shortprices::Union{Vector{Dict}, Nothing}, medendvaluesdicts::Vector{Dict}, shorts::DArray, reltol::Float64, tnormal::ProbTime, scenarios::ScenarioModellingMethod, skipmed::Millisecond, stochastictimes::DArray, stepnr::Int, settings::Dict)
    @sync @distributed for core in 1:max(numcores-1,1)
        stochastictime = localpart(stochastictimes)
        masters = localpart(masters_)
        subs = localpart(subs_)
        states = localpart(states_)
        cuts = localpart(cuts_)
        storagevalues = localpart(storagevalues_)
        short = localpart(shorts)

        localix = 0
        for range in localindices(shorts)
            for ix in range
                localix += 1
                if !(!short[localix] && (skipmed.value > 0))
                    stochastic!(masters[localix], subs[localix], states[localix], cuts[localix], storagevalues[localix], startstates, medprices, shortprices, medendvaluesdicts, short[localix], reltol, tnormal, scenarios, stochastictime[localix], stepnr, settings)
                end
            end
        end
        # if skipmed.value == 0
        #     GC.gc() # regulary garbage collection, else there is always one processor core for each time step that the others have to wait for
        # end
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
function distribute_subsystems(ustoragesystemobjects::Vector{Tuple{Vector, Vector{Vector}}}, ushorts::Vector, ustoragevalues::Vector)
    objectcount = [length(first(objs)) for objs in ustoragesystemobjects]
    sortedindex = sortperm(objectcount, rev=true)
    testsort = distribute(collect(1:length(ustoragesystemobjects)))
    storagesystemobjects = Vector(undef,length(objectcount))
    shorts = Vector(undef,length(ushorts))
    storagevalues = Vector(undef,length(ustoragevalues))
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
        storagevalues[collect(range)[localix]] = ustoragevalues[sortedindex][i]
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
    storagesystemobjects1 = distribute(storagesystemobjects)
    shorts1 = distribute(shorts, storagesystemobjects1)
    storagevalues1 = distribute(storagevalues, storagesystemobjects1)
    return storagesystemobjects1, shorts1, storagevalues1
end
function distribute_subsystems_flat(ustoragesystemobjects::Vector{Tuple{Vector, Vector{Vector}}}, ushorts::Vector, ustoragevalues::Vector)
    storagesystemobjects = distribute(ustoragesystemobjects)
    shorts = distribute(ushorts, storagesystemobjects)
    storagevalues = distribute(ustoragevalues, storagesystemobjects)
    return storagesystemobjects, shorts, storagevalues
end



