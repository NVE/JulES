# Read json files into dataelements
function json_to_elements(path::String, filename::String)
    parsed = JSON.parsefile(joinpath(path, filename))
    return getelements(parsed, path)
end

# Add element for scenariotimeperiod
function addscenariotimeperiod!(elements::Vector{DataElement}, start::Int, stop::Int)
    push!(elements, getelement(TIMEPERIOD_CONCEPT, "ScenarioTimePeriod", "ScenarioTimePeriod", 
            ("Start", getisoyearstart(start)), ("Stop", getisoyearstart(stop))))
end

# Insert horizons into commodities. E.g. all batteries will have the power horizon, since they interact with the power market
function set_horizon!(elements::Vector{DataElement}, commodity::String, horizon::Horizon)
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
function addStartEqualStopAllStorages!(modelobjects::Dict)
    for obj in values(modelobjects)
        if obj isa BaseStorage
            trait = StartEqualStop(obj)
            modelobjects[getid(trait)] = trait
        end
    end
end

# Power balances needs slack variable for when the inelastic supply (wind, solar, RoR) is higher than the inelastic demand
function addPowerUpperSlack!(modelobjects::Dict) # add after object manipulation
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
function remove_startupcosts!(modelobjects::Dict)
    for (id,obj) in modelobjects
        if obj isa StartUpCost
            delete!(modelobjects, id)
        end
    end
end

# Remove start-up costs. Does not make sense to have them when the horizon does not have a fine time-resolution
function remove_transmissionramping!(modelobjects::Dict)
    for (id,obj) in modelobjects
        if obj isa TransmissionRamping
            delete!(modelobjects, id)
        end
    end
end

# Remove hydroramping. Does not make sense to have them when the horizon does not have a fine time-resolution
function remove_hydrorampingwithout!(modelobjects::Dict)
    for (id,obj) in modelobjects
        if obj isa HydroRampingWithout
            delete!(modelobjects, id)
        end
    end
end
function remove_hydroramping!(modelobjects::Dict)
    for (id,obj) in modelobjects
        if obj isa HydroRamping
            delete!(modelobjects, id)
        end
    end
end

# Set start and end reservoir as a percentage of capacity
function setstartstoragepercentage!(prob::Prob, storages::Vector, start::ProbTime, percentage::Int)
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

function setendstoragepercentage!(prob::Prob, storages::Vector, endtime::ProbTime, percentage::Int)
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

# Initialize dict of statevariables from list of modelobjects
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

# Set start and end states for objects with statevariables
function setstartstates!(prob::Prob, objects::Vector, startstates::Dict)
    for obj in objects

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startstates[getinstancename(first(getvarout(statevariable)))]
        end
        
        setingoingstates!(prob, states)
    end
end

function setendstates!(prob::Prob, objects::Vector, startstates::Dict)
    for obj in objects

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startstates[getinstancename(first(getvarout(statevariable)))]
        end
        
        setoutgoingstates!(prob, states)
    end
end

function getnonstorageobjects(modelobjects::Vector)
    nonstorageobjects = []
    for obj in modelobjects
        if getconceptname(getid(obj)) != "Storage"
            if length(getstatevariables(obj)) > 1
                error("Not supported")
            elseif length(getstatevariables(obj)) == 1
                push!(nonstorageobjects, obj)
            end
        end
    end
    return nonstorageobjects
end

# Initialize dict of statevariables that are not storages from list of modelobjects
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

# Get dual value of a storage at a specific time period
function getinsideduals(p::Prob, storages::Vector, t::Int)
    endvalues = zeros(Float64,length(storages))
    for (i, storage) in enumerate(storages)
        bid = getid(getbalance(storage))
        endvalues[i] = -getcondual(p, bid, t)
    end
    return endvalues
end