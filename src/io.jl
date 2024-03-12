"""
Definition of default input and output types
"""


struct DefaultJulESInput <: AbstractJulESInput
    mainconfig::Dict
    settings::Dict
    onlysubsystemmodels::Bool

    steps::Int
    steplength::Millisecond
    simstarttime::ProbTime
    datascenarios::Vector{Scenarios}

    tnormaltype::String
    tphaseintype::String
    phaseinoffset::Millisecond
    phaseindelta::Millisecond
    phaseinsteps::Int

    elements::Vector
    progelements::Vector

    function DefaultJulESInput(dataset, config)
        mainconfig = config["main"]
        settings = config[mainconfig["settings"]]

        onlysubsystemmodel = false
        if !haskey(settings["problems"], "prognosis") && haskey(settings["problems"], "stochastic") && !haskey(settings["problems"], "clearing")
            onlysubsystemmodel = true
        end

        println("Time parameters")
        @time timeparams = gettimeparams(mainconfig, settings)
        steps, steplength, simstarttime, datascenarios, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps = timeparams

        println("Get data")
        @time begin
            elements = dataset["elements"]
            detailedrescopl = dataset["detailedrescopl"]
            addscenariotimeperiod_vector!(elements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
    
            if haskey(dataset, "progelements")
                progelements = dataset["progelements"]
                addscenariotimeperiod_vector!(progelements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
            else
                progelements = elements
            end
        end

        horizons = gethorizons(config)

        subsystems = getsubsystems(dataset, config, horizons)::Dict{String, Vector{DataElement}}

        # TODO: use throw-away-dummy-modelobjs for data proc

        return new(mainconfig, settings, onlysubsystemmodels,
            steps, steplength, simstarttime, datascenarios,
            tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps,
            elements, progelements)
    end
end

function getscenarios(datayear::Int64, weatheryear::Int64, weekstart::Int64, datanumscen::Int64, simtimetype::String)
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

    return (simtime, datascenarios)
end

function gettimeparams(mainconfig::Dict, settings::Dict)
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
    simstarttime, datascenarios = getscenarios(datayear, weatheryear, weekstart, datanumscen, simtimetype)

    return (steps, steplength, simstarttime, datascenarios, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps)
end

function _getscenariotime(simtime::ProbTime, scenario::Scenario, input::AbstractJulESInput, inputtime::String)
    scenarios = input.scenarios
    phaseinoffset = input.phaseinoffset
    phaseindelta = input.phaseindelta
    phaseinsteps = input.phaseinsteps
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

function getscenariotnormal(simtime::ProbTime, scenario::Scenario, input::AbstractJulESInput)
    timetype = input.tnormaltype
    return _getscenariotime(simtime, scenario, input, timetype)
end
function getscenariotphasein(simtime::ProbTime, scenario::Scenario, input::AbstractJulESInput)
    timetype = input.tphasein
    return _getscenariotime(simtime, scenario, input, timetype)
end


function getscenmodmethod(problem::Dict, numscen::Int64, objects::Vector)
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




# -------------------------------------------------------------------------------------------
struct DefaultJulESOutput <: AbstractJulESOutput
    # TODO: complete
end

# TODO: complete
get_cores(input) = nothing
get_horizons(input) = nothing
get_simulation_period(input) = nothing
get_startstates_ppp(input) = nothing



# Should live here and not in slot in PricePrognosisProblem?
function get_medendvaluesdict(input)
end

