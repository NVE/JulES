"""
Definition of default input and output types
"""

struct DefaultJulESInput <: AbstractJulESInput
    dataset::Dict
    mainconfig::Dict
    settings::Dict
    onlysubsystemmodel::Bool # TODO: can probably remove this

    steps::Int
    steplength::Millisecond
    simstarttime::ProbTime
    scenmod_data::Vector{Scenarios}

    tnormaltype::String
    tphaseintype::String
    phaseinoffset::Millisecond
    phaseindelta::Millisecond
    phaseinsteps::Int

    horizons::Dict{Tuple{ScenarioIx, TermName, CommodityName}, Horizon}

    function DefaultJulESInput(dataset, config)
        mainconfig = config["main"]
        settings = config[mainconfig["settings"]]

        onlysubsystemmodel = false
        if !haskey(settings["problems"], "prognosis") && !haskey(settings["problems"], "endvalue") && haskey(settings["problems"], "stochastic") && !haskey(settings["problems"], "clearing")
            onlysubsystemmodel = true
        end

        println("Time parameters")
        @time timeparams = get_timeparams(mainconfig, settings)
        steps, steplength, simstarttime, scenmod_data, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps = timeparams

        println("Handle elements")
        @time begin
            elements = dataset["elements"]
            addscenariotimeperiod_vector!(elements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
    
            if haskey(dataset, "progelements")
                progelements = dataset["progelements"]
                addscenariotimeperiod_vector!(progelements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
            end

            enekvglobaldict = Dict{String,Float64}()
            if !onlysubsystemmodel
                for element in elements
                    if element.typename == GLOBALENEQKEY
                        enekvglobaldict[split(element.instancename,"GlobalEneq_")[2]] = element.value["Value"]
                    end
                end
            end
            dataset["enekvglobaldict"] = enekvglobaldict
        end

        horizons = get_horizons(config)

        return new(dataset, mainconfig, settings, onlysubsystemmodels,
            steps, steplength, simstarttime, scenmod_data,
            tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps,
            horizons)
    end
end

get_dataset(input::DefaultJulESInput) = input.dataset
get_elements(input::DefaultJulESInput) = get_dataset(input)["elements"]
get_elements_ppp(input::DefaultJulESInput) = get_dataset(input)["elements_ppp"]

get_mainconfig(input::DefaultJulESInput) = input.mainconfig
get_settings(input::DefaultJulESInput) = input.settings
get_onlysubsystemmodel(input::DefaultJulESInput) = input.onlysubsystemmodel

get_steps(input::DefaultJulESInput) = input.steps
get_steplength(input::DefaultJulESInput) = input.steplength

get_simstarttime(input::DefaultJulESInput) = input.simstarttime
get_scenmod_data(input::DefaultJulESInput) = input.scenmod_data
get_numscen_data(input::DefaultJulESInput) = length(input.scenmod_data.scenarios)

get_tnormaltype(input::DefaultJulESInput) = input.tnormaltype
get_tphaseintype(input::DefaultJulESInput) = input.tphaseintype
get_phaseinoffset(input::DefaultJulESInput) = input.phaseinoffset
get_phaseindelta(input::DefaultJulESInput) = input.phaseindelta
get_phaseinsteps(input::DefaultJulESInput) = input.phaseinsteps

get_horizons(input::DefaultJulESInput) = input.horizons

function get_datascenarios(datayear::Int64, weatheryear::Int64, weekstart::Int64, datanumscen::Int64, simtimetype::String)
    # Standard time for market clearing - perfect information so simple time type
    datasimtime = getisoyearstart(datayear) + Week(weekstart-1)
    weathersimtime = getisoyearstart(weatheryear) + Week(weekstart-1)
    simtime = gettnormal(simtimetype, datasimtime, weathersimtime)

    # Make scenariooffset for all uncertainty scenarios
    datascenarios = Vector{Scenario}(undef, datanumscen)
    for scen in 1:datanumscen
        weatherscenariotime = getisoyearstart(weatheryear + scen - 1) + Week(weekstart-1)
        weatheroffset = weatherscenariotime - weathersimtime
        datascenarios[scen] = Scenarios(weatheroffset, 1/datanumscen, Dict{String, Tuple{Any,Float64}}, scen)
    end

    return (simtime, NoScenarioModellingMethod(datascenarios))
end

function get_timeparams(mainconfig::Dict, settings::Dict)
    weekstart = mainconfig["weekstart"]
    
    weatheryearstart = settings["time"]["weatheryearstart"]
    weatheryearstop = settings["time"]["weatheryearstop"]
    datanumscen = scenarioyearstop - scenarioyearstart # scenarios to consider uncertainty for
    
    simulationyears = mainconfig["simulationyears"]
    extrasteps = mainconfig["extrasteps"]
    steplength = Millisecond(Hour(settings["time"]["steplength_hours"]))
    steps = Int(ceil((getisoyearstart(datayear + simulationyears) - getisoyearstart(datayear)).value/steplength.value) + extrasteps);
    
    # Phasein settings
    phaseinoffset = steplength # phase in straight away from second stage scenarios
    phaseindelta = Millisecond(Day(settings["time"]["probtime"]["phaseindelta_days"])) # Phase in the second stage scenario over 5 weeks
    phaseinsteps = settings["time"]["probtime"]["phaseinsteps"] # Phase in second stage scenario in 5 steps

    # Make standard time and scenario uncertainty times
    tnormaltype = settings["time"]["probtime"]["normaltime"]
    tphaseintype = settings["time"]["probtime"]["phaseintime"]
    simstarttime, scenmod_data = get_datascenarios(datayear, weatheryear, weekstart, datanumscen, simtimetype)

    return (steps, steplength, simstarttime, scenmod_data, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps)
end

function _get_scenariotime(simtime::ProbTime, scenario::Scenario, input::AbstractJulESInput, inputtime::String)
    phaseinoffset = get_phaseinoffset(input)
    phaseindelta = get_phaseindelta(input)
    phaseinsteps = get_phaseinsteps(input)
    datasimtime = getdattime(simtime)
    weathersimtime = getscenariotime(simtime)
    weatherscenariotime = getscenariotime(simtime) + scenario.weatheroffset

    if timetype == "PrognosisTime"
        return PrognosisTime(datasimtime, datasimtime, weatherscenariotime)
    elseif timetype == "FixedDataTwoTime"
        return FixedDataTwoTime(datasimtime, weatherscenariotime)
    elseif timetype == "PhaseinPrognosisTime"
        return PhaseinPrognosisTime(datasimtime, datasimtime, weathersimtime, weatherscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
    elseif timetype == "PhaseinFixedDataTwoTime"
        return PhaseinFixedDataTwoTime(datasimtime, weathersimtime, weatherscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
    elseif timetype == "PrognosisTime"
        return PrognosisTime(datasimtime, datasimtime, weatherscenariotime)
    elseif timetype == "FixedDataTwoTime"
        return FixedDataTwoTime(datasimtime, weatherscenariotime)
    else
        error("$timetype not implementet in getscenariotime-function")
    end
end

function get_scentnormal(simtime::ProbTime, scenario::Scenario, input::AbstractJulESInput)
    timetype = get_tnormaltype(input)
    return _get_scentime(simtime, scenario, input, timetype)
end
function get_scentphasein(simtime::ProbTime, scenario::Scenario, input::AbstractJulESInput)
    timetype = get_tphaseintype(input)
    return _get_scentime(simtime, scenario, input, timetype)
end


function get_scenmod(problem::Dict, numscen::Int64, objects::Vector)
    method = problem["function"]
    if method == "InflowClusteringMethod"
        parts = problem["parts"] # divide scendelta into this many parts, calculate sum inflow for each part of the inflow series, then use clustering algorithm
        scendelta = MsTimeDelta(Day(problem["scendelta"]))
        return InflowClusteringMethod(numscen, objects, parts, scendelta)
    elseif method == "SumInflowQuantileMethod"
        a = problem["a"]
        b = problem["b"]
        c = problem["c"]
        maxquantile = problem["maxquantile"]
        scendelta = problem["scendelta"]
        usedensity = MsTimeDelta(Day(problem["usedensity"]))
        return SumInflowQuantileMethod(numscen, objects, maxquantile, a, b, c, scendelta, usedensity=usedensity)
    else
        error("$method not supported")
    end
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

function get_numscen_sim(input::AbstractJulESInput)
    settings = get_settings(input)
    if haskey(settings["scenariogeneration"], "simulation")
        return settings["scenariogeneration"]["simulation"]["numscen"]
    else
        return get_numscen_data(input)
    end
end
function get_numscen_ppp(input::AbstractJulESInput)
    settings = get_settings(input)
    if haskey(settings["scenariogeneration"], "prognosis")
        return settings["scenariogeneration"]["prognosis"]["numscen"]
    else
        return get_numscen_sim(input)
    end
end
function get_numscen_evp(input::AbstractJulESInput)
    settings = get_settings(input)
    if haskey(settings["scenariogeneration"], "endvalue")
        return settings["scenariogeneration"]["endvalue"]["numscen"]
    else
        return get_numscen_ppp(input)
    end
end
function get_numscen_sp(input::AbstractJulESInput)
    settings = get_settings(input)
    if haskey(settings["scenariogeneration"], "stochastic")
        return settings["scenariogeneration"]["stochastic"]["numscen"]
    else
        return get_numscen_evp(input)
    end
end

function get_simperiod(input::AbstractJulESInput)
    t = get_simstarttime(input)
    N = get_steps(input)
    delta = get_steplength(input)
    skipmed = Millisecond(Hour(0))
    skipmax = Millisecond(Hour(delta*(db.settings["time"]["skipmax"]-1)))

    return (t, N, delta, skipmed, skipmax)

function get_aggzone(settings::Dict)
    if haskey(settings["problems"], "aggzone")
        return settings["problems"]["aggzone"]
    else
        return Dict()
    end
end

# -------------------------------------------------------------------------------------
# Other inpututils not used yet
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

function getoutputindex(mainconfig::Dict, datayear::Int64, scenarioyear::Int64)
    if mainconfig["outputindex"] == "datayear"
        return datayear
    elseif mainconfig["outputindex"] == "scenarioyear"
        return scenarioyear
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

# -------------------------------------------------------------------------------------------

# TODO: complete
get_cores(input) = nothing   # should return non-empty CorId[]
get_horizons(input) = nothing # should return Dict{Tuple{TermName, CommodityName}, Horizon}

# Should live here and not in slot in PricePrognosisProblem?
function get_medendvaluesdict(input)
end


struct DefaultJulESOutput <: AbstractJulESOutput
    # TODO: complete
end
