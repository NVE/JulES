"""
Definition of default input and output types
"""

struct DefaultJulESInput <: AbstractJulESInput
    cores::Vector{CoreId}
    dataset::Dict
    mainconfig::Dict
    settings::Dict
    onlysubsystemmodel::Bool # TODO: can probably remove this

    steps::Int
    steplength::Millisecond
    simstarttime::ProbTime
    scenmod_data::AbstractScenarioModellingMethod

    tnormaltype::String
    tphaseintype::String
    phaseinoffset::Millisecond
    phaseindelta::Millisecond
    phaseinsteps::Int

    horizons::Dict{Tuple{TermName, CommodityName}, Horizon}

    function DefaultJulESInput(config, dataset, datayear, weatheryear)
        mainconfig = config["main"]
        settings = config[mainconfig["settings"]]
        numcores = mainconfig["numcores"]
        cores = collect(1:numcores)

        onlysubsystemmodel = false
        if !haskey(settings["problems"], "prognosis") && !haskey(settings["problems"], "endvalue") && haskey(settings["problems"], "stochastic") && !haskey(settings["problems"], "clearing")
            onlysubsystemmodel = true
        end

        println("Time parameters")
        @time timeparams = get_timeparams(mainconfig, settings, datayear, weatheryear)
        steps, steplength, simstarttime, scenmod_data, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps = timeparams

        println("Handle elements")
        @time begin
            elements = dataset["elements"]
            add_scenariotimeperiod_int!(elements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
    
            if haskey(dataset, "progelements")
                progelements = dataset["progelements"]
                add_scenariotimeperiod_int!(progelements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
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

        horizons = get_horizons(settings, datayear)

        return new(cores, dataset, mainconfig, settings, onlysubsystemmodel,
            steps, steplength, simstarttime, scenmod_data,
            tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps,
            horizons)
    end
end

get_cores(input::DefaultJulESInput) = input.cores
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
    simtime = get_tnormal(simtimetype, datasimtime, weathersimtime)

    # Make scenariooffset for all uncertainty scenarios
    datascenarios = Vector{AbstractScenario}(undef, datanumscen)
    for scen in 1:datanumscen
        weatherscenariotime = getisoyearstart(weatheryear + scen - 1) + Week(weekstart-1)
        weatheroffset = weatherscenariotime - weathersimtime
        datascenarios[scen] = WeatherScenario(weatheroffset, 1/datanumscen, scen)
    end

    return (simtime, NoScenarioModellingMethod(datascenarios))
end

function get_timeparams(mainconfig::Dict, settings::Dict, datayear::Int, weatheryear::Int)
    weekstart = mainconfig["weekstart"]
    
    weatheryearstart = settings["time"]["weatheryearstart"]
    weatheryearstop = settings["time"]["weatheryearstop"]
    datanumscen = weatheryearstop - weatheryearstart # scenarios to consider uncertainty for
    
    simulationyears = mainconfig["simulationyears"]
    extrasteps = mainconfig["extrasteps"]
    steplength = parse_duration(settings["horizons"]["clearing"], "termduration")
    steps = Int(ceil((getisoyearstart(datayear + simulationyears) - getisoyearstart(datayear)).value/steplength.value) + extrasteps);
    
    # Phasein settings
    phaseinoffset = steplength # phase in straight away from second stage scenarios
    phaseindelta = parse_duration(settings["time"]["probtime"], "phaseindelta") # Phase in the second stage scenario over 5 weeks
    phaseinsteps = settings["time"]["probtime"]["phaseinsteps"] # Phase in second stage scenario in 5 steps

    # Make standard time and scenario uncertainty times
    tnormaltype = settings["time"]["probtime"]["normaltime"]
    tphaseintype = settings["time"]["probtime"]["phaseintime"]
    simstarttime, scenmod_data = get_datascenarios(datayear, weatheryear, weekstart, datanumscen, tnormaltype)

    return (steps, steplength, simstarttime, scenmod_data, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps)
end

function get_tnormal(type::String, datatime::DateTime, scenariotime::DateTime)
    if type == "PrognosisTime"
        return PrognosisTime(datatime, datatime, scenariotime)
    elseif type == "FixedDataTwoTime"
        return FixedDataTwoTime(datatime, scenariotime)
    else
        error("$type not implementet in get_normal-function")
    end
end

function get_scenariotime(simtime::ProbTime, scenario::AbstractScenario, input::AbstractJulESInput, normal_phasein::String)
    phaseinoffset = get_phaseinoffset(input)
    phaseindelta = get_phaseindelta(input)
    phaseinsteps = get_phaseinsteps(input)
    datasimtime = getdattime(simtime)
    weathersimtime = getscenariotime(simtime)
    weatherscenariotime = getscenariotime(simtime) + scenario.weatheroffset

    if normal_phasein == "normaltime"
        timetype = get_settings(input)["time"]["probtime"]["normaltime"]
    else
        @assert normal_phasein == "phaseintime"
        timetype = get_settings(input)["time"]["probtime"]["phaseintime"]
    end

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

function add_scenariotimeperiod_int!(elements::Vector{DataElement}, start::Int, stop::Int)
    push!(elements, getelement(TIMEPERIOD_CONCEPT, "ScenarioTimePeriod", "ScenarioTimePeriod", 
            ("Start", getisoyearstart(start)), ("Stop", getisoyearstart(stop))))
end

function get_scentnormal(simtime::ProbTime, scenario::AbstractScenario, input::AbstractJulESInput)
    timetype = get_tnormaltype(input)
    return _get_scentime(simtime, scenario, input, timetype)
end
function get_scentphasein(simtime::ProbTime, scenario::AbstractScenario, input::AbstractJulESInput)
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
function get_numscen_stoch(input::AbstractJulESInput)
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
    skipmax = Millisecond(Hour(delta*(input.settings["time"]["skipmax"]-1)))

    return (t, N, delta, skipmed, skipmax)
end

function get_aggzone(settings::Dict)
    if haskey(settings["problems"], "aggzone")
        return settings["problems"]["aggzone"]
    else
        return Dict()
    end
end

function get_statedependentprod(settings::Dict)
    if haskey(settings, "statedependentprod")
        return settings["statedependentprod"]
    else
        return false
    end
end

function get_statedependentpump(settings::Dict)
    if haskey(settings, "statedependentpump")
        return settings["statedependentpump"]
    else
        return false
    end
end

function get_headlosscost(settings::Dict)
    if haskey(settings, "headlosscost")
        return settings["headlosscost"]
    else
        return false
    end
end

# -------------------------------------------------------------------------------------
# Other inpututils not used yet
function get_onlyagghydro(settings::Dict)
    if haskey(settings["problems"], "onlyagghydro")
        return settings["problems"]["onlyagghydro"]
    else
        return false
    end
end

function getoutputindex(mainconfig::Dict, datayear::Int64, weatheryear::Int64)
    if mainconfig["outputindex"] == "datayear"
        return datayear
    elseif mainconfig["outputindex"] == "weatheryear"
        return weatheryear
    end
end

# Prognosis util functions
function getrhsdata(rhsdata::Dict, datayear::Int64, weatheryearstart::Int64, weatheryearstop::Int64)
    method = rhsdata["function"]
    if method == "DynamicExogenPriceAHData"
        return DynamicExogenPriceAHData(Id("Balance", rhsdata["balance"])) # TODO: If dynamic use tphasein
    elseif method == "StaticRHSAHData"
        return StaticRHSAHData("Power", datayear, weatheryearstart, weatheryearstop)
    elseif method == "DynamicRHSAHData"
        return DynamicRHSAHData("Power")
    else
        error("$method not supported")
    end
end

# -------------------------------------------------------------------------------------------

function get_horizons(settings, datayear)
    horizons = Dict{Tuple{TermName, CommodityName}, Horizon}()
    commoditites = settings["horizons"]["commodities"]
    n_durations = Dict{Tuple{TermName, CommodityName}, Tuple{Int, Millisecond}}()

    for term in keys(settings["horizons"])
        if term != "commodities"
            for commodity in commoditites
                n_durations = get_n_durations(term, commodity, settings)

                method = settings["horizons"][term][commodity]["function"]
                if method == "SequentialHorizon"
                    horizons[(term, commodity)] = build_sequentialhorizon(term, commodity, settings, n_durations)
                elseif method == "AdaptiveHorizon"
                    horizons[(term, commodity)] = build_adaptivehorizon(term, commodity, settings, n_durations, datayear)
                end
            end
        end
    end

    return horizons
end

function build_sequentialhorizon(term, commodity, settings, n_durations)
    horizon = SequentialHorizon(n_durations...)

    # TODO: Shrinkable
    # if settings["horizons"][term]["shrinkable"]
    #     horizon = ShrinkableHorizon(horizon, )
    # end

    return horizon
end

function build_adaptivehorizon(term, commodity, settings, n_durations, datayear)
    weatheryearstart = settings["time"]["weatheryearstart"]
    weatheryearstop = settings["time"]["weatheryearstop"]
    rhsdata = getrhsdata(settings["horizons"][term][commodity]["rhsdata"], datayear, weatheryearstart, weatheryearstop)
    rhsmethod = parse_methods(settings["horizons"][term][commodity]["rhsmethod"])
    clusters = settings["horizons"][term][commodity]["clusters"]
    unitduration = Millisecond(Hour(settings["horizons"][term][commodity]["unitduration_hours"]))

    horizon = AdaptiveHorizon(clusters, unitduration, rhsdata, rhsmethod, n_durations...)

    # TODO: Shrinkable
    # if settings["horizons"][term]["shrinkable"]
    #     horizon = ShrinkableHorizon(horizon)
    # end

    return horizon
end

function get_shrinkable(settings::Dict)
    if haskey(settings, "shrinkable")
        return settings["shrinkable"]
    else
        return false
    end
end

function get_n_durations(term, commodity, settings)
    maintermconfig = settings["horizons"][term]
    maincommodityconfig = maintermconfig[commodity]

    method = maincommodityconfig["function"]
    if method == "AdaptiveHorizon"
        maincommodityconfig = maintermconfig[maincommodityconfig["macro"]]
    end

    mainperiodduration = parse_duration(maincommodityconfig, "periodduration")

    if term == LongTermName
        subterms = [ClearingTermName, ShortTermName, MedTermName, LongTermName]
    elseif term == MedTermName
        subterms = [ClearingTermName, ShortTermName, MedTermName]
    elseif term == ShortTermName
        subterms = [ClearingTermName, ShortTermName]
    elseif term == ClearingTermName
        subterms = [ClearingTermName]
    end

    n_durations = []
    for subterm in subterms
        subtermduration = parse_duration(settings["horizons"][subterm], "termduration")

        n = subtermduration.value / mainperiodduration.value
        if n < 1
            push!(n_durations, 1)
            push!(n_durations, subtermduration)
        else
            if !isinteger(n)
                error("Period $mainperiodduration don't fit termduration $subtermduration for term $term, subterm $subterm and commodity $commodity")
            end
            push!(n_durations, Int(n))
            push!(n_durations, mainperiodduration)
        end
    end

    return n_durations
end

function parse_duration(config, namestart)
    for key in keys(config)
        if startswith(key, namestart)
            res = split(key, "_")[2]
            if res == "hours"
                return Millisecond(Hour(config[key]))
            elseif res == "days"
                return Millisecond(Day(config[key]))
            elseif res == "weeks"
                return Millisecond(Week(config[key]))
            end
        end
    end
    println(config)
    error("Key $namestart not in config")
end


# -----------------------------------------------------------
struct DefaultJulESOutput <: AbstractJulESOutput
    cores::Vector{CoreId}
    function DefaultJulESOutput(input)
        cores = get_cores(input)
        return new(cores)
    end
end
get_cores(output::DefaultJulESOutput) = output.cores
cleanup_output(output) = nothing
