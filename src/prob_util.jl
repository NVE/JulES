# Add element for scenariotimeperiod
function add_scenariotimeperiod_vector!(elements::Vector{DataElement}, start::Int, stop::Int)
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
end

# The hydropower storages in the dataset needs boundary conditions for the state variables
function add_StartEqualStopAllStorages!(modelobjects::Dict)
    for obj in values(modelobjects)
        if obj isa BaseStorage
            trait = StartEqualStop(obj)
            modelobjects[getid(trait)] = trait
        end
    end
end

# Power balances needs slack variable for when the inelastic supply (wind, solar, RoR) is higher than the inelastic demand
function add_PowerUpperSlack!(modelobjects::Dict) # add after object manipulation
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
function set_startstoragepercentage!(prob::Prob, storages::Vector, start::ProbTime, percentage::Float64)
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

function set_endstoragepercentage!(prob::Prob, storages::Vector, endtime::ProbTime, percentage::Float64)
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
function get_states(modelobjects::Vector)
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

# Startstates-------------------------------------------------------------------------------------
function update_startstates(stepnr, t)
    db = get_local_db()
    settings = get_settings(db)

    if stepnr == 1
        if haskey(settings["problems"], "prognosis")
            dummyobjects_ppp = first(db.dummyobjects_ppp)
            dummystorages_ppp = getstorages(dummyobjects_ppp)
            get_startstates!(db.startstates, settings["problems"]["prognosis"], get_dataset(db), dummyobjects_ppp, dummystorages_ppp, t)
            startstates_max!(dummystorages_ppp, t, db.startstates)
        end
        if haskey(settings["problems"], "endvalue")
            dummyobjects = first(db.dummyobjects)
            dummystorages = getstorages(dummyobjects)
            get_startstates!(db.startstates, settings["problems"]["endvalue"], get_dataset(db), dummyobjects, dummystorages, t)
            startstates_max!(dummystorages, t, db.startstates)
        end
        bothequal = false
        if haskey(settings["problems"], "endvalue") && haskey(settings["problems"], "stochastic")
            if settings["problems"]["endvalue"] == settings["problems"]["stochastic"]
                bothequal = true
            end
        end
        if haskey(settings["problems"], "stochastic") && !bothequal
            dummyobjects = first(db.dummyobjects)
            dummystorages = getstorages(dummyobjects)
            get_startstates!(db.startstates, settings["problems"]["stochastic"], get_dataset(db), dummyobjects, dummystorages, t)
            startstates_max!(dummystorages, t, db.startstates)
        end
    else
        get_startstates_from_cp(db.startstates, db.core_cp)
    end
end

function get_startstates_from_cp(startstates, core)
    f = @spawnat core get_startstates_from_cp()
    startstates_cp = fetch(f)

    for (k, v) in startstates_cp
        startstates[k] = v
    end
    return 
end

function setstartstates!(p::Prob, startstates::Dict{String, Float64})
    storages = getstorages(getobjects(p))
    set_startstates!(p, storages, startstates)
    return
end

function get_startstates!(startstates::Dict, problemconfig::Dict, dataset::Dict, objects::Dict, storages::Vector, tnormal::ProbTime)
    startstorages = problemconfig["startstorages"]
    if startstorages["function"] == "percentages"
        shorttermstorages = getshorttermstorages(collect(values(objects)), Hour(problemconfig["shorttermstoragecutoff_hours"]))
        longtermstorages = setdiff(storages, shorttermstorages)
        merge!(startstates, getstartstoragepercentage(shorttermstorages, tnormal, startstorages["shortpercentage"]))
        merge!(startstates, getstartstoragepercentage(longtermstorages, tnormal, startstorages["longpercentage"]))
    elseif startstorages["function"] == "percentage"
        merge!(startstates, getstartstoragepercentage(storages, tnormal, startstorages["percentage"]))
    elseif haskey(dataset, startstorages["function"])
        merge!(startstates, dataset[startstorages["function"]])
    end
end

# Initialize max startstates and cap at maximum
function startstates_max!(objects::Vector, t::ProbTime, startstates::Dict)
    for obj in objects
        resname = getinstancename(getid(obj))
        if haskey(startstates, resname)
            startstates[resname * "_max"] = getparamvalue(obj.ub, t, MsTimeDelta(Millisecond(0)))
            if startstates[resname] > startstates[resname * "_max"] # TODO: Add warning or logging
                startstates[resname] = startstates[resname * "_max"]
            end
        end
    end
    
    return startstates
end

# Set start and end states for objects with statevariables
function set_startstates!(prob::Prob, objects::Vector, startstates::Dict)
    for obj in objects

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startstates[getinstancename(first(getvarout(statevariable)))]
        end
        
        setingoingstates!(prob, states)
    end
end

function set_endstates!(prob::Prob, objects::Vector, startstates::Dict)
    for obj in objects

        states = Dict{StateVariableInfo, Float64}()
        for statevariable in getstatevariables(obj)
            states[statevariable] = startstates[getinstancename(first(getvarout(statevariable)))]
        end
        
        setoutgoingstates!(prob, states)
    end
end

function get_startstoragepercentage(storages::Vector, start::ProbTime, percentage::Float64)
    startstates = Dict{String, Float64}()

    for obj in storages
        dummydelta = MsTimeDelta(Millisecond(0))
        startreservoir = getparamvalue(getub(obj), start, dummydelta)*percentage/100
        for statevariable in getstatevariables(obj)
            startstates[getinstancename(first(getvarout(statevariable)))] = startreservoir
        end
    end

    return startstates
end

function get_nonstorageobjects(modelobjects::Vector)
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
function get_nonstoragestatevariables(modelobjects::Vector)
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
function get_insideduals(p::Prob, storages::Vector, t::Int)
    endvalues = zeros(Float64,length(storages))
    for (i, storage) in enumerate(storages)
        bid = getid(getbalance(storage))
        endvalues[i] = -getcondual(p, bid, t)
    end
    return endvalues
end

# Find first exogen price in a vector of model objects
function find_firstprice(objects)
    for obj in objects
        if obj isa ExogenBalance
            return getprice(obj)
        end
    end
end;