# Read json files into dataelements
function json_to_elements(path::String, filename::String)
    parsed = JSON.parsefile(joinpath(path, filename))
    return getelements(parsed, path)
end

# Add element for scenariotimeperiod
function addscenariotimeperiod_vector!(elements::Vector{DataElement}, start::Int, stop::Int)
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
function setstartstoragepercentage!(prob::Prob, storages::Vector, start::ProbTime, percentage::Float64)
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

function setendstoragepercentage!(prob::Prob, storages::Vector, endtime::ProbTime, percentage::Float64)
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
function getstates(modelobjects::Vector)
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

function getstartstoragepercentage(storages::Vector, start::ProbTime, percentage::Float64)
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

# Time util functions
function gettnormal(type::String, datatime::DateTime, scenariotime::DateTime)
    if type == "PrognosisTime"
        return PrognosisTime(datatime, datatime, scenariotime)
    elseif type == "FixedDataTwoTime"
        return FixedDataTwoTime(datatime, scenariotime)
    else
        error("$type not implementet in getnormal-function")
    end
end

function gettphasein(type::String, datatime::DateTime, scenariotime::DateTime, uncertaintyscenariotime::DateTime, phaseinoffset::Millisecond, phaseindelta::Millisecond, phaseinsteps::Int64)
    if type == "PhaseinPrognosisTime"
        return PhaseinPrognosisTime(datatime, datatime, scenariotime, uncertaintyscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
    elseif type == "PhaseinFixedDataTwoTime"
        return PhaseinFixedDataTwoTime(datatime, scenariotime, uncertaintyscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
    elseif type == "PrognosisTime"
        return PrognosisTime(datatime, datatime, scenariotime)
    elseif type == "FixedDataTwoTime"
        return FixedDataTwoTime(datatime, scenariotime)
    else
        error("$type not implementet in gettphasein-function")
    end
end

function getprobtimes(datayear::Int64, weekstart::Int64, scenarioyear::Int64, datanumscen::Int64, tnormaltype::String, tphaseintype::String, phaseinoffset::Millisecond, phaseindelta::Millisecond, phaseinsteps::Int64)
    # Standard time for market clearing - perfect information so simple time type
    datatime = getisoyearstart(datayear) + Week(weekstart-1)
    scenariotime = getisoyearstart(scenarioyear) + Week(weekstart-1)
    tnormal = gettnormal(tnormaltype, datatime, scenariotime)

    # Make scenario times for all uncertainty scenarios. List of tuples with tnormal, tphasein and scenarionumber
    datascentimes = []
    for scen in 1:datanumscen
        uncertaintyscenariotime = getisoyearstart(scenarioyear + scen - 1) + Week(weekstart-1)
        scentnormal = gettnormal(tnormaltype, datatime, uncertaintyscenariotime)
        scentphasein = gettphasein(tphaseintype, datatime, scenariotime, uncertaintyscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
       
        push!(datascentimes, (scentnormal, scentphasein, scen))
    end
    datascenmodmethod = NoScenarioModellingMethod(datanumscen, datascentimes)
    return (tnormal, datascenmodmethod)
end

# Parse methods (alternative to eval(Meta.parse))
function parse_methods(s::String)
    if s == "HiGHS_Prob()"
        return HiGHS_Prob()
    elseif s == "HighsSimplexMethod()"
        return HighsSimplexMethod()
    elseif s == "HighsSimplexMethod(warmstart=false)"
        return HighsSimplexMethod(warmstart=false)
    elseif s == "HighsSimplexSIPMethod(warmstart=false)"
        return HighsSimplexSIPMethod(warmstart=false)
    elseif s == "KMeansAHMethod()"
        return KMeansAHMethod()
    end
end

# Prognosis util functions
function getrhsdata(rhsdata::Dict, datayear::Int64, scenarioyearstart::Int64, scenarioyearstop::Int64)
    method = rhsdata["function"]
    if method == "DynamicExogenPriceAHData"
        return DynamicExogenPriceAHData(Id("Balance", rhsdata["balance"])) # TODO: If dynamic use tphasein
    elseif method == "StaticRHSAHData"
        return StaticRHSAHData("Power", datayear, scenarioyearstart, scenarioyearstop)
    elseif method == "DynamicRHSAHData"
        return DynamicRHSAHData("Power")
    else
        error("$method not supported")
    end
end

function getscenmodmethod(problem::Dict, numscen::Int64)
    method = problem["function"]
    if method == "InflowClusteringMethod"
        parts = problem["parts"] # divide scendelta into this many parts, calculate sum inflow for each part of the inflow series, then use clustering algorithm
        return InflowClusteringMethod(numscen, parts)
    else
        error("$method not supported")
    end
end

function getstartstates!(startstates::Dict, problemsconfig::Dict, problem::String, dataset::Dict, objects::Dict, storages::Vector, tnormal::ProbTime)
    startstorages = problemsconfig[problem]["startstorages"]
    if startstorages["function"] == "percentages"
        shorttermstorages = getshorttermstorages(collect(values(objects)), Hour(problemsconfig["shorttermstoragecutoff_hours"]))
        longtermstorages = setdiff(storages, shorttermstorages)
        merge!(startstates, getstartstoragepercentage(shorttermstorages, tnormal, startstorages["shortpercentage"]))
        merge!(startstates, getstartstoragepercentage(longtermstorages, tnormal, startstorages["longpercentage"]))
    elseif startstorages["function"] == "percentage"
        merge!(startstates, getstartstoragepercentage(storages, tnormal, startstorages["percentage"]))
    elseif haskey(dataset, startstorages["function"])
        merge!(startstates, dataset[startstorages["function"]])
    end
end

# Get dictionary with each detailed reservoir and their water value for each scenario
# TODO: Detailed run-of-river reservoirs get water value from aggregated reservoir hydro
function getendvaluesdicts(endvaluesobjs::Any, detailedrescopl::Dict, enekvglobaldict::Dict)
    endvaluesdicts = Dict[];
    for endvaluesobj in endvaluesobjs
        instance = [getinstancename(getid(obj)) for obj in endvaluesobj.objects]
        endvalues = endvaluesobj.values
        endvaluesdict = Dict(instance .=> endvalues)

        for (k,v) in detailedrescopl
            endvaluesdict["Reservoir_" * k] = endvaluesdict["Reservoir_" * v * "_hydro_reservoir"] * enekvglobaldict[k]
        end
        push!(endvaluesdicts, endvaluesdict)
    end
    
    return endvaluesdicts
end

function getoutputindex(mainconfig::Dict, datayear::Int64, scenarioyear::Int64)
    if mainconfig["outputindex"] == "datayear"
        return datayear
    elseif mainconfig["outputindex"] == "scenarioyear"
        return scenarioyear
    end
end

# Find first exogen price in a vector of model objects
function findfirstprice(objects)
    for obj in objects
        if obj isa ExogenBalance
            return getprice(obj)
        end
    end
end;

# Get aggzone from config
function getaggzone(settings::Dict)
    if haskey(settings["problems"], "aggzone")
        return settings["problems"]["aggzone"]
    else
        return Dict()
    end
end

# Get if onlyagghydro
function getonlyagghydro(settings::Dict)
    if haskey(settings["problems"], "onlyagghydro")
        return settings["problems"]["onlyagghydro"]
    else
        return false
    end
end

# Get if statedependentprod
function getstatedependentprod(settings::Dict)
    if haskey(settings, "statedependentprod")
        return settings["statedependentprod"]
    else
        return false
    end
end

# Get if statedependentpump
function getstatedependentpump(settings::Dict)
    if haskey(settings, "statedependentpump")
        return settings["statedependentpump"]
    else
        return false
    end
end

# Get if statedependentpump
function getheadlosscost(settings::Dict)
    if haskey(settings, "headlosscost")
        return settings["headlosscost"]
    else
        return false
    end
end