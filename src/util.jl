# Read json files into dataelements
function json_to_elements(path, filename)
    parsed = JSON.parsefile(joinpath(path, filename))
    return getelements(parsed, path)
end

# Add element for scenariotimeperiod
function addscenariotimeperiod!(elements, start, stop)
    push!(elements, getelement(TIMEPERIOD_CONCEPT, "ScenarioTimePeriod", "ScenarioTimePeriod", 
            ("Start", getisoyearstart(start)), ("Stop", getisoyearstart(stop))))
end

# Insert horizons into commodities. E.g. all batteries will have the power horizon, since they interact with the power market
function set_horizon!(elements, commodity, horizon)
    # If element already exist, replace horizon with new
    for element in elements
        if element.typename == "BaseCommodity"
            if element.instancename == commodity
                element.value[HORIZON_CONCEPT] = horizon
                return
            end
        end
    end
    
    # Else, add commodity to element list
    push!(elements, getelement(COMMODITY_CONCEPT, "BaseCommodity", commodity, 
        (HORIZON_CONCEPT, horizon)))
    elements
end

# The hydropower storages in the dataset needs boundary conditions for the state variables
function addStartEqualStopAllStorages!(modelobjects)
    for obj in values(modelobjects)
        if obj isa BaseStorage
            trait = StartEqualStop(obj)
            modelobjects[getid(trait)] = trait
        end
    end
end

# Power balances needs slack variable for when the inelastic supply (wind, solar, RoR) is higher than the inelastic demand
function addPowerUpperSlack!(modelobjects) # add after object manipulation
    for obj in values(modelobjects)
        if obj isa BaseBalance
            if getid(getcommodity(obj)) == Id("Commodity", "Power")
                balancename = getinstancename(getid(obj))
                
                varname = "SlackVar_" * balancename
                varkey = Id(FLOW_CONCEPT, varname)
                var = BaseFlow(varkey)
                
                sethorizon!(var, gethorizon(obj))
                setlb!(var, LowerZeroCapacity())
                
                arrowname = "SlackArrow_" * balancename
                arrowkey = Id(ARROW_CONCEPT, arrowname) 
                arrow = BaseArrow(arrowkey, obj, BaseConversion(PlusOneParam()), 0)
                addarrow!(var, arrow)
                
                modelobjects[varkey] = var
            end
        end 
    end
end

# Remove start-up costs. Does not make sense to have them when the horizon does not have a fine time-resolution
function remove_startupcosts!(modelobjects)
    for (id,obj) in modelobjects
        if obj isa StartUpCost
            delete!(modelobjects, id)
        end
    end
end

# Set start reservoir as a percentage of capacity
function setstartstoragepercentage!(prob, storages, start, percentage)
    for obj in storages

        dummydelta = MsTimeDelta(Millisecond(0))
        startreservoir = getparamvalue(getub(obj), start, dummydelta)*percentage/100

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startreservoir
        end
        
        setingoingstates!(prob, states)
    end
end

function setendstoragepercentage!(prob, storages, endtime, percentage)
    for obj in storages

        dummydelta = MsTimeDelta(Millisecond(0))
        endreservoir = getparamvalue(getub(obj), endtime, dummydelta)*percentage/100
        
        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = endreservoir
        end
        
        setoutgoingstates!(prob, states)
    end
end

function setstartstates!(prob, storages, startstates)
    for obj in storages

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startstates[getinstancename(first(getvarout(statevariable)))]
        end
        
        setingoingstates!(prob, states)
    end
end

function setendstates!(prob, storages, startstates)
    for obj in storages

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startstates[getinstancename(first(getvarout(statevariable)))]
        end
        
        setoutgoingstates!(prob, states)
    end
end

function getstatevariables(modelobjects::Vector)
    states = Dict{StateVariableInfo, Float64}()
    
    for obj in modelobjects
        if length(getstatevariables(obj)) > 1
            error("Not supported")
        else
            for statevariable in getstatevariables(obj)
                states[statevariable] = 0.0
            end
        end
    end
    return states
end

function getnonstoragestatevariables(modelobjects::Vector)
    states = Dict{StateVariableInfo, Float64}()
    
    for obj in modelobjects
        if getconceptname(getid(obj)) != "Storage"
            if length(getstatevariables(obj)) > 1
                error("Not supported")
            else
                for statevariable in getstatevariables(obj)
                    states[statevariable] = 0.0
                end
            end
        end
    end
    return states
end

function getinsideduals(p::Prob, storages::Vector, t::Int)
    endvalues = zeros(Float64,length(storages))
    for (i, storage) in enumerate(storages)
        bid = getid(getbalance(storage))
        endvalues[i] = getcondual(p, bid, t)
    end
    return endvalues
end

# Distribute subsystems after number of objects in subsystem
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

