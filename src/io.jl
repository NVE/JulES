"""
Definition of default input and output types
"""

struct DefaultJulESInput <: AbstractJulESInput
    cores::Vector{CoreId}
    dataset::Dict
    mainconfig::Dict
    settings::Dict
    datayear::Int
    weatheryear::Int
    onlysubsystemmodel::Bool # TODO: can probably remove this
    
    steps::Int
    steplength::Millisecond
    simstarttime::TuLiPa.ProbTime
    scenmod_data::AbstractScenarioModellingMethod

    tnormaltype::String
    tphaseintype::String
    phaseinoffset::Millisecond
    phaseindelta::Millisecond
    phaseinsteps::Int

    horizons::Dict{Tuple{TermName, CommodityName}, TuLiPa.Horizon}

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
    
            if haskey(dataset, "elements_ppp")
                elements_ppp = dataset["elements_ppp"]
                add_scenariotimeperiod_int!(elements_ppp, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
            end

            iprogtype = get(dataset, "iprogtype", "direct")
            useifm = iprogtype != "direct"
            if useifm
                ifm_elements = dataset["ifm_elements"]
                add_scenariotimeperiod_int!(ifm_elements, settings["time"]["weatheryearstart"], settings["time"]["weatheryearstop"])
            end

            enekvglobaldict = Dict{String,Float64}()
            if !onlysubsystemmodel
                for element in elements
                    if element.typename == TuLiPa.GLOBALENEQKEY
                        enekvglobaldict[split(element.instancename,"GlobalEneq_")[2]] = element.value["Value"]
                    end
                end
            end
            dataset["enekvglobaldict"] = enekvglobaldict
        end

        horizons = get_horizons(settings, datayear)

        return new(cores, dataset, mainconfig, settings, datayear, weatheryear, onlysubsystemmodel,
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
get_datayear(input::DefaultJulESInput) = input.datayear
get_weatheryear(input::DefaultJulESInput) = input.weatheryear
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
"""

Returns the subsystem distribution method from the config file. If there is no distribution method, the default method "bysize" is returned. 
"Bysize" is the orginal method where the subsystems are sorted from biggest to smallest, and then distributed forward and then backward until all subsystems are distributed. 
"""
function get_distribution_method_mp(input::DefaultJulESInput, default::String="bysize")
    settings = get_settings(input)
    # Retrieve the distribution method value
    if !get_onlysubsystemmodel(input)
        method = get(settings["problems"]["stochastic"],"distribution_method_mp", default)
    else
        method = "core_main"
    end
    # Check if the method is not nothing and not an empty string
    if !isnothing(method) && !isempty(method)
        return method
    else
        return default
    end
end


"""
Returns the scenario problem distribution method from the config file. If there is no distribution method, the default method "withmp" is returned. 
Withmp is the original method where scenarios are distributed on the same core as the master problem. 
"""
function get_distribution_method_sp(input::DefaultJulESInput, default::String="withmp")
    settings = get_settings(input)
    
    # Retrieve the distribution method value
    if !get_onlysubsystemmodel(input)
        method = get(settings["problems"]["stochastic"],"distribution_method_sp", default)
    else
        method = "even"
    end

    # Check if the method is not nothing and not an empty string
    if !isnothing(method) && !isempty(method)
        return method
    else
        return default
    end
end

get_iprogtype(input::DefaultJulESInput) = get(input.dataset, "iprogtype", "direct")
has_ifm_results(input::DefaultJulESInput) = get_iprogtype(input) != "direct"
get_ifm_normfactors(input::DefaultJulESInput) = get(input.dataset, "ifm_normfactors", Dict{String, Float64}())
get_ifm_elements(input::DefaultJulESInput) = get(input.dataset, "ifm_elements", JulES.TuLiPa.DataElement[])

function get_ifm_names(input::DefaultJulESInput)
    if haskey(input.dataset, "ifm_names")
        s1 = Set(input.dataset["ifm_names"])
        s2 = Set([e.instancename for e in input.dataset["ifm_elements"] if e.conceptname == ABSTRACT_INFLOW_MODEL])
        return String[i for i in intersect(s1, s2)]
    else
        return String[]
    end
end

function get_ifm_weights(input::DefaultJulESInput)
    if haskey(input.dataset, "ifm_weights")
        w = input.dataset["ifm_weights"]
        # remove stations from w that does not exist
        # and update weights accordingly
        names = get_ifm_names(input)
        missings_dict = Dict()
        for k in keys(w)
            for station in keys(w[k])
                if !(station in names)
                    if haskey(missings_dict, k) == false
                        missings_dict[k] = Set()
                    end
                    push!(missings_dict[k], station)
                end
            end
        end
        for k in keys(missings_dict)
            sum_nonmissing = sum(weight for (station, weight) in w[k] if !(station in missings_dict[k]))
            for (station, weight) in w[k]
                if !(station in missings_dict[k])
                    w[k][station] = weight / sum_nonmissing
                end
            end
        end
        for k in keys(missings_dict)
            for station in missings_dict[k]
                delete!(w[k], station)
            end
        end
        for (__, weights) in w
            @assert isapprox(round(sum(values(weights)); digits=4), 1.0)
        end

        return w
    else
        return Dict{String, Dict{String, Float64}}()
    end
end 


function get_datascenarios(datayear::Int64, weatheryear::Int64, weekstart::Int64, datanumscen::Int64, simtimetype::String)
    # Standard time for market clearing - perfect information so simple time type
    datasimtime = TuLiPa.getisoyearstart(datayear) + Week(weekstart-1)
    weathersimtime = TuLiPa.getisoyearstart(weatheryear) + Week(weekstart-1)
    simtime = get_tnormal(simtimetype, datasimtime, weathersimtime)

    # Make scenariooffset for all uncertainty scenarios
    datascenarios = Vector{WeatherScenario}(undef, datanumscen)
    for scen in 1:datanumscen
        weatherscenariotime = TuLiPa.getisoyearstart(weatheryear + scen - 1) + Week(weekstart-1)
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
    steplength = get_steplength(settings)
    steps = Int(ceil((TuLiPa.getisoyearstart(datayear + simulationyears) - TuLiPa.getisoyearstart(datayear)).value/steplength.value) + extrasteps);
    
    # Phasein settings
    phaseinoffset = steplength # phase in straight away from second stage scenarios
    if haskey(settings["time"]["probtime"], "phaseintime")
        phaseindelta = parse_duration(settings["time"]["probtime"], "phaseindelta") # Phase in the second stage scenario over 5 weeks
    else
        phaseindelta = Millisecond(0)
    end
    phaseinsteps = get(settings["time"]["probtime"], "phaseinsteps", 0) # Phase in second stage scenario in 5 steps

    # Make standard time and scenario uncertainty times
    tnormaltype = settings["time"]["probtime"]["normaltime"]
    if haskey(settings["time"]["probtime"], "phaseintime")
        tphaseintype = settings["time"]["probtime"]["phaseintime"]
    else
        tphaseintype = tnormaltype
    end
    simstarttime, scenmod_data = get_datascenarios(datayear, weatheryear, weekstart, datanumscen, tnormaltype)

    return (steps, steplength, simstarttime, scenmod_data, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps)
end

function get_steplength(settings)
    if haskey(settings["horizons"], "clearing")
        return steplength = parse_duration(settings["horizons"]["clearing"], "termduration")
    elseif haskey(settings["horizons"], "master")
        return steplength = parse_duration(settings["horizons"]["master"], "termduration")
    else
        error("No key clearing or master for settings[horizons]")
    end
end

function get_tnormal(type::String, datatime::DateTime, scenariotime::DateTime)
    if type == "PrognosisTime"
        return TuLiPa.PrognosisTime(datatime, datatime, scenariotime)
    elseif type == "FixedDataTwoTime"
        return TuLiPa.FixedDataTwoTime(datatime, scenariotime)
    else
        error("$type not implementet in get_normal-function")
    end
end

function get_scentime(simtime::TuLiPa.ProbTime, scenario::AbstractScenario, input::AbstractJulESInput, timetype::String)
    phaseinoffset = get_phaseinoffset(input)
    phaseindelta = get_phaseindelta(input)
    phaseinsteps = get_phaseinsteps(input)
    datasimtime = TuLiPa.getdatatime(simtime)
    weathersimtime = TuLiPa.getscenariotime(simtime)
    weatherscenariotime = TuLiPa.getscenariotime(simtime) + scenario.weatheroffset

    if timetype == "PrognosisTime"
        return TuLiPa.PrognosisTime(datasimtime, datasimtime, weatherscenariotime)
    elseif timetype == "FixedDataTwoTime"
        return TuLiPa.FixedDataTwoTime(datasimtime, weatherscenariotime)
    elseif timetype == "PhaseinPrognosisTime"
        return TuLiPa.PhaseinPrognosisTime(datasimtime, datasimtime, weathersimtime, weatherscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
    elseif timetype == "PhaseinFixedDataTwoTime"
        return TuLiPa.PhaseinFixedDataTwoTime(datasimtime, weathersimtime, weatherscenariotime, phaseinoffset, phaseindelta, phaseinsteps)
    else
        error("$timetype not implementet in getscenariotime-function")
    end
end

function add_scenariotimeperiod_int!(elements::Vector{TuLiPa.DataElement}, start::Int, stop::Int)
    push!(elements, TuLiPa.getelement(TuLiPa.TIMEPERIOD_CONCEPT, "ScenarioTimePeriod", "ScenarioTimePeriod", 
            ("Start", TuLiPa.getisoyearstart(start)), ("Stop", TuLiPa.getisoyearstart(stop))))
end

function get_scentnormal(simtime::TuLiPa.ProbTime, scenario::AbstractScenario, input::AbstractJulESInput)
    timetype = get_tnormaltype(input)
    return get_scentime(simtime, scenario, input, timetype)
end
function get_scentphasein(simtime::TuLiPa.ProbTime, scenario::AbstractScenario, input::AbstractJulESInput)
    timetype = get_tphaseintype(input)
    return get_scentime(simtime, scenario, input, timetype)
end


function get_scenmod(allscenarios::Vector, problem::Dict, numscen::Int64, objects::Vector)
    method = problem["function"]
    if method == "InflowClusteringMethod"
        scenarios = Vector{eltype(allscenarios)}(undef, numscen)
        parts = problem["parts"] # divide scendelta into this many parts, calculate sum inflow for each part of the inflow series, then use clustering algorithm
        scendelta = TuLiPa.MsTimeDelta(parse_duration(problem, "scendelta"))
        return InflowClusteringMethod(scenarios, objects, parts, scendelta)
    elseif method == "SumInflowQuantileMethod"
        scenarios = Vector{eltype(allscenarios)}(undef, numscen)
        a = problem["a"]
        b = problem["b"]
        c = problem["c"]
        maxquantile = problem["maxquantile"]
        scendelta = TuLiPa.MsTimeDelta(parse_duration(problem, "scendelta"))
        usedensity = problem["usedensity"]
        return SumInflowQuantileMethod(scenarios, objects, maxquantile, a, b, c, scendelta, usedensity=usedensity)
    else
        error("$method not supported")
    end
end

# Parse methods (alternative to eval(Meta.parse))
function parse_methods(s::String)
    if s == "HiGHS_Prob()"
        return TuLiPa.HiGHS_Prob()
    elseif s == "HighsSimplexMethod()"
        return TuLiPa.HighsSimplexMethod()
    elseif s == "HighsSimplexMethod(warmstart=false)"
        return TuLiPa.HighsSimplexMethod(warmstart=false)
    elseif s == "HighsSimplexSIPMethod(warmstart=false)"
        return TuLiPa.HighsSimplexSIPMethod(warmstart=false)
    elseif s == "JuMPHiGHSMethod()"
        return TuLiPa.JuMPHiGHSMethod()
    elseif s == "KMeansAHMethod()"
        return TuLiPa.KMeansAHMethod()
    elseif s == "PercentilesAHMethod()"
        return TuLiPa.PercentilesAHMethod()
    end
end

function get_numscen_sim(input::AbstractJulESInput)
    settings = get_settings(input)
    if !isnothing(settings["scenariogeneration"])
        if haskey(settings["scenariogeneration"], "simulation")
            return settings["scenariogeneration"]["simulation"]["numscen"]
        end
    end
    return get_numscen_data(input)
end
function get_numscen_stoch(input::AbstractJulESInput)
    settings = get_settings(input)
    if !isnothing(settings["scenariogeneration"])
        if haskey(settings["scenariogeneration"], "stochastic")
            return settings["scenariogeneration"]["stochastic"]["numscen"]
        end
    end
    return get_numscen_sim(input)
end

function get_simperiod(input::AbstractJulESInput)
    t = get_simstarttime(input)
    N = get_steps(input)
    steplength = get_steplength(input)
    skipmed = Millisecond(Hour(0))
    skipmax = Millisecond(Hour(steplength*(input.settings["time"]["skipmax"]-1)))

    return (t, N, steplength, skipmed, skipmax)
end

function get_aggzone(settings::Dict)
    if haskey(settings["problems"], "aggzone")
        return settings["problems"]["aggzone"]
    else
        return Dict()
    end
end

function get_aggzonecopl(aggzone::Dict)
    aggzonecopl = Dict()
    for (k,v) in aggzone
        for vv in v
            aggzonecopl["PowerBalance_" * vv] = "PowerBalance_" * k
        end
    end

    return aggzonecopl
end

get_result_prices_ppp(settings::Dict)::Int = get(settings["results"], "prices_ppp", 0)

has_statedependentprod(settings::Dict)::Bool = get(settings, "statedependentprod", false)
has_statedependentpump(settings::Dict)::Bool = get(settings, "statedependentpump", false)
has_headlosscost(settings::Dict)::Bool = get(settings, "statedependentpump", false)
has_onlyagghydro(settings::Dict)::Bool = get(settings["problems"], "onlyagghydro", false)

has_result_times(settings::Dict)::Bool = get(settings["results"], "times", false)
has_result_scenarios(settings::Dict)::Bool = get(settings["results"], "scenarios", false)
has_result_memory(settings::Dict)::Bool = get(settings["results"], "memory", false)
has_result_storagevalues(settings::Dict)::Bool = get(settings["results"], "storagevalues", false)
has_result_hydrolevels_water(settings::Dict)::Bool = get(settings["results"], "hydrolevels_water", false)

function get_outputindex(mainconfig::Dict, datayear::Int64, weatheryear::Int64)
    if mainconfig["outputindex"] == "datayear"
        return datayear
    elseif mainconfig["outputindex"] == "weatheryear"
        return weatheryear
    end
end

# Prognosis util functions
function get_rhsdata(rhsdata::Dict, datayear::Int64, weatheryearstart::Int64, weatheryearstop::Int64)
    method = rhsdata["function"]
    if method == "DynamicExogenPriceAHData"
        return TuLiPa.DynamicExogenPriceAHData(TuLiPa.Id("Balance", rhsdata["balance"])) # TODO: If dynamic use tphasein
    elseif method == "FindFirstDynamicExogenPriceAHData"
        return TuLiPa.FindFirstDynamicExogenPriceAHData() # TODO: If dynamic use tphasein
    elseif method == "StaticRHSAHData"
        return TuLiPa.StaticRHSAHData("Power", datayear, weatheryearstart, weatheryearstop)
    elseif method == "DynamicRHSAHData"
        return TuLiPa.DynamicRHSAHData("Power")
    else
        error("$method not supported")
    end
end

# -------------------------------------------------------------------------------------------

function get_horizons(settings, datayear)
    horizons = Dict{Tuple{TermName, CommodityName}, TuLiPa.Horizon}()
    commoditites = settings["horizons"]["commodities"]
    n_durations = Dict{Tuple{TermName, CommodityName}, Tuple{Int, Millisecond}}()

    for term in keys(settings["horizons"])
        if !(term in ["commodities", "shrinkable"])
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
    horizon = TuLiPa.SequentialHorizon(n_durations...)

    if get_shrinkable(settings["horizons"][term])
        startafter = parse_duration(settings["horizons"]["shrinkable"], "startafter")
        shrinkatleast = parse_duration(settings["horizons"]["shrinkable"], "shrinkatleast")
        minperiod = get_steplength(settings)
        horizon = TuLiPa.ShrinkableHorizon(horizon, startafter, shrinkatleast, minperiod)
    end

    return horizon
end

function build_adaptivehorizon(term, commodity, settings, n_durations, datayear)
    weatheryearstart = settings["time"]["weatheryearstart"]
    weatheryearstop = settings["time"]["weatheryearstop"]
    rhsdata = get_rhsdata(settings["horizons"][term][commodity]["rhsdata"], datayear, weatheryearstart, weatheryearstop)
    rhsmethod = parse_methods(settings["horizons"][term][commodity]["rhsmethod"])
    clusters = settings["horizons"][term][commodity]["clusters"]
    unitduration = Millisecond(Hour(settings["horizons"][term][commodity]["unitduration_hours"]))

    horizon = TuLiPa.AdaptiveHorizon(clusters, unitduration, rhsdata, rhsmethod, n_durations...)

    if get_shrinkable(settings["horizons"][term])
        startafter = parse_duration(settings["horizons"]["shrinkable"], "startafter")
        shrinkatleast = parse_duration(settings["horizons"]["shrinkable"], "shrinkatleast")
        minperiod = get_steplength(settings)
        horizon = TuLiPa.ShrinkableHorizon(horizon, startafter, shrinkatleast, minperiod)
    end

    return horizon
end

function get_shrinkable(horizonterm::Dict)
    if haskey(horizonterm, "shrinkable")
        return horizonterm["shrinkable"]
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
    elseif term == MasterTermName
        subterms = [MasterTermName]
    elseif term == SubTermName
        subterms = [MasterTermName, SubTermName]
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
mutable struct DefaultJulESOutput <: AbstractJulESOutput
    timing_ppp::Dict
    timing_evp::Dict
    timing_mp::Dict
    timing_sp::Dict
    timing_cp::Array

    storagevalues::Dict

    prices_balances::Vector{TuLiPa.Id}
    prices_long::Array{Float64}
    deltas_long::Array{Float64}
    prices_med::Array{Float64}
    deltas_med::Array{Float64}
    prices_short::Array{Float64}
    deltas_short::Array{Float64}

    scenweights_sim::Array{Float64}
    scenweights_stoch::Array{Float64}

    prices::Array{Float64}
    rhstermvalues::Array{Float64}
    production::Array{Float64}
    consumption::Array{Float64}
    hydrolevels::Array{Float64}
    batterylevels::Array{Float64}
    othervalues::Dict
    
    modelobjects::Dict
    powerbalances::Vector
    rhsterms::Vector
    rhstermbalances::Vector
    plants::Vector
    plantbalances::Vector
    plantarrows::Dict
    demands::Vector
    demandbalances::Vector
    demandarrows::Dict
    hydrostorages::Vector
    batterystorages::Vector
    otherobjects::Dict
    otherbalances::Dict

    statenames::Vector{String}
    statematrix::Array{Float64} # end states after each step

    ifm_stations::Vector{String}
    ifm_statenames::Vector{String}
    ifm_u0::Vector{Matrix{Float64}}
    ifm_Q::Matrix{Float64}

    function DefaultJulESOutput(input)
        return new(Dict(),Dict(),Dict(),Dict(),[],
        Dict(),
        [],[],[],[],[],[],[],
        [],[],
        [],[],[],[],[],[],Dict(),
        Dict(),[],[],[],[],[],Dict(),[],[],Dict(),[],[],Dict(),Dict(),
        [],[],
        [],[],[],Matrix{Float64}(undef, (0,0)))
    end
end

get_cores(output::DefaultJulESOutput) = output.cores

function init_local_output()
    db = get_local_db()
    settings = get_settings(db)
    db.output = get_output_from_input(db.input)

    steps = get_steps(db)

    if has_result_times(settings)
        for (scenix, core) in db.dist_ppp
            db.output.timing_ppp[scenix] = zeros(steps, 3, 3) # TODO: Matrix more flexible long term, and array instead of dict?
        end

        for (scenix, subix, core) in db.dist_evp
            db.output.timing_evp[(scenix, subix)] = zeros(steps, 3)
        end

        for (subix, core) in db.dist_mp
            db.output.timing_mp[subix] = zeros(steps, 5)
        end

        for (scenix, subix, core) in db.dist_sp
            db.output.timing_sp[(scenix, subix)] = zeros(steps, 3)
        end

        db.output.timing_cp = zeros(steps, 3)
    end

    if has_result_scenarios(settings)
        db.output.scenweights_sim = zeros(steps, get_numscen_sim(db.input))
        db.output.scenweights_stoch = zeros(steps, get_numscen_stoch(db.input))
    end

    if has_result_storagevalues(settings)
        for (subix, core) in db.dist_mp
            if settings["results"]["storagevalues"]
                if has_headlosscost(settings["problems"]["stochastic"]["master"])
                    num_storagevalues = get_numscen_stoch(db.input)*2 + 2 # scenarios + master operative + master operative after headlosscost adjustment
                else
                    num_storagevalues = get_numscen_stoch(db.input)*2 + 1 # scenarios + master operative 
                end
                if haskey(settings["problems"], "clearing")
                    num_storagevalues += 1
                end
                if haskey(settings["problems"], "stochastic")
                    num_storagevalues += get_numscen_stoch(db.input)
                end
                if haskey(settings["problems"], "endvalue")
                    num_storagevalues += get_numscen_stoch(db.input)
                end
                f = @spawnat core get_numstates(subix)
                db.output.storagevalues[subix] = zeros(steps, num_storagevalues, fetch(f))
            end
        end
    end

    collect_interval = get_result_prices_ppp(settings)
    if collect_interval != 0
        collect_steps = div(steps, collect_interval) + 1

        for (id, obj) in first(get_dummyobjects_ppp(db))
            if obj isa TuLiPa.BaseBalance
                if TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(obj))) == "Power"
                    push!(db.output.prices_balances, id)
                end
            end
        end
        num_balances = length(db.output.prices_balances)

        numperiods_long = TuLiPa.getnumperiods(get_horizons(db.input)[(LongTermName, "Power")])
        db.output.prices_long = zeros(collect_steps, num_balances, get_numscen_stoch(db.input), numperiods_long)
        db.output.deltas_long = zeros(collect_steps, numperiods_long)

        numperiods_med = TuLiPa.getnumperiods(get_horizons(db.input)[(MedTermName, "Power")])
        db.output.prices_med = zeros(collect_steps, num_balances, get_numscen_stoch(db.input), numperiods_med)
        db.output.deltas_med = zeros(collect_steps, numperiods_med)

        numperiods_short = TuLiPa.getnumperiods(get_horizons(db.input)[(ShortTermName, "Power")])
        db.output.prices_short = zeros(collect_steps, num_balances, get_numscen_stoch(db.input), numperiods_short)
        db.output.deltas_short = zeros(collect_steps, numperiods_short)

        @sync for (scenix, core) in db.dist_ppp
            @spawnat core init_prices_ppp(scenix, num_balances, numperiods_long, numperiods_med, numperiods_short)
        end
    end

    if has_ifm_results(db.input)
        db.output.ifm_stations = collect(get_ifm_names(db.input))
        db.output.ifm_statenames = collect(get_ifm_statenames())
        num_states = get_ifm_numstates()
        num_stations = length(db.output.ifm_stations)
        @assert length(db.output.ifm_u0) == 0
        for __ in 1:num_states
            push!(db.output.ifm_u0, zeros(Float64, (num_stations, steps)))
        end
        db.output.ifm_Q = zeros(Float64, (num_stations, steps))
    end
end

function get_ifm_numstates()
    # TODO: find common numstates by calling get_numstates on each ifm and verify all ifm of same type
    return 2
end

function get_ifm_statenames()
    # TODO: find common statenames by calling get_statenames on each ifm and verify all ifm of same type
    return ["snow", "ground"]
end

function init_prices_ppp(scenix, num_balances, numperiods_long, numperiods_med, numperiods_short)
    db = get_local_db()

    div = db.ppp[scenix].div
    div["prices_long"] = zeros(num_balances, numperiods_long)
    div["prices_med"] = zeros(num_balances, numperiods_med)
    div["prices_short"] = zeros(num_balances, numperiods_short)
    div["deltas_long"] = zeros(numperiods_long)
    div["deltas_med"] = zeros(numperiods_med)
    div["deltas_short"] = zeros(numperiods_short)
    return
end

function get_numstates(subix)
    db = get_local_db()
    return length(db.mp[subix].states)
end

function collect_ifm_u0(stepnr)
    db = get_local_db()
    d = Dict{String, Vector{Float64}}()
    for core in get_cores(db.input)
        fetched = fetch(@spawnat core local_collect_ifm_u0(stepnr))
        if fetched isa RemoteException
            throw(fetched)
        end
        for (name, x) in fetched
            @assert !haskey(d, name)
            d[name] = x
        end
    end
    return d
end

function local_collect_ifm_u0(stepnr)
    db = get_local_db()
    d = Dict{String, Vector{Float64}}()
    for (name, core) in db.dist_ifm
        if core == db.core
            (stored_stepnr, u0) = db.div[IFM_DB_STATE_KEY][name]
            @assert stored_stepnr == stepnr
            d[name] = u0
        end
    end
    return d
end

function collect_ifm_Q(stepnr)
    db = get_local_db()
    d = Dict{String, Float64}()
    for core in get_cores(db.input)
        fetched = fetch(@spawnat core local_collect_ifm_Q(stepnr))
        if fetched isa RemoteException
            throw(fetched)
        end
        for (name, x) in fetched
            @assert !haskey(d, name)
            d[name] = x
        end
    end
    return d
end

function local_collect_ifm_Q(stepnr)
    db = get_local_db()
    d = Dict{String, Float64}()
    for (name, core) in db.dist_ifm
        if core == db.core
            (stored_stepnr, Q) = db.div[IFM_DB_FLOW_KEY][name]
            @assert stored_stepnr == stepnr
            d[name] = Q
        end
    end
    return d
end

function update_output(t::TuLiPa.ProbTime, stepnr::Int)
    db = get_local_db()
    settings = get_settings(db)
    steps = get_steps(db)

    if has_ifm_results(db.input)
        u0 = collect_ifm_u0(stepnr)
        for (i, station) in enumerate(db.output.ifm_stations)
            for (j, v) in enumerate(u0[station])
                db.output.ifm_u0[j][i, stepnr] = v
            end
        end
        ifm_Q = collect_ifm_Q(stepnr)
        for (i, station) in enumerate(db.output.ifm_stations)
            db.output.ifm_Q[i, stepnr] = ifm_Q[station]
        end
    end

    if has_result_times(settings)
        for (scenix, core) in db.dist_ppp
            f = @spawnat core get_maintiming_ppp(scenix)
            db.output.timing_ppp[scenix][stepnr, :, :] .= fetch(f)
            @spawnat core reset_maintiming_ppp(scenix)
        end

        for (scenix, subix, core) in db.dist_evp
            f = @spawnat core get_maintiming_evp(scenix, subix)
            db.output.timing_evp[(scenix, subix)][stepnr, :] .= fetch(f)
            @spawnat core reset_maintiming_evp(scenix, subix)
        end

        for (subix, core) in db.dist_mp
            f = @spawnat core get_maintiming_mp(subix)
            db.output.timing_mp[subix][stepnr, :] .= fetch(f)
            @spawnat core reset_maintiming_mp(subix)
        end

        for (scenix, subix, core) in db.dist_sp
            f = @spawnat core get_maintiming_sp(scenix, subix)
            db.output.timing_sp[(scenix, subix)][stepnr, :] .= fetch(f)
            @spawnat core reset_maintiming_sp(scenix, subix)
        end
    end    

    if has_result_scenarios(settings)
        db.output.scenweights_sim[stepnr, :] .= [get_probability(scen) for scen in get_scenarios(db.scenmod_sim)]
        db.output.scenweights_stoch[stepnr, :] .= [get_probability(scen) for scen in get_scenarios(db.scenmod_stoch)]
    end

    if has_result_storagevalues(settings)
        for (subix, core) in db.dist_mp
            f = @spawnat core get_storagevalues_stoch(subix)
            storagevalues_stoch = fetch(f)
            dim = (size(storagevalues_stoch, 1))
            db.output.storagevalues[subix][stepnr, 1:dim, :] .= storagevalues_stoch

            if haskey(settings["problems"], "clearing")
                cutid = fetch(@spawnat core get_cutsid(subix))
                cuts = get_obj_from_id(TuLiPa.getobjects(db.cp.prob), cutid)
                for (j, statevar) in enumerate(cuts.statevars) # master / operative water values after headlosscost
                    obj = get_obj_from_id(TuLiPa.getobjects(db.cp.prob), first(TuLiPa.getvarout(statevar))) # TODO: OK to assume objid = varoutid?
                    balance = TuLiPa.getbalance(obj)

                    db.output.storagevalues[subix][stepnr, dim+1, j] = TuLiPa.getcondual(db.cp.prob, TuLiPa.getid(balance), TuLiPa.getnumperiods(TuLiPa.gethorizon(balance)))
                    if haskey(balance.metadata, TuLiPa.GLOBALENEQKEY)
                        db.output.storagevalues[subix][stepnr, dim+1, j] = db.output.storagevalues[subix][stepnr, dim+1, j] / balance.metadata[TuLiPa.GLOBALENEQKEY]
                    end
                    if haskey(settings["problems"], "stochastic")
                        for scenix in 1:get_numscen_stoch(db.input)
                            core_stoch = get_core_sp(db, scenix, subix)
                            f = @spawnat core_stoch get_enddual_stoch(scenix, subix, first(TuLiPa.getvarout(statevar)))
                            db.output.storagevalues[subix][stepnr, dim+1+scenix, j] = fetch(f)
                        end
                    end
                    if haskey(settings["problems"], "endvalue") && is_subsystem_evp(db.subsystems[subix])
                        for scenix in 1:get_numscen_stoch(db.input)
                            core_evp = get_core_evp(db, scenix, subix)
                            f = @spawnat core_evp get_enddual_evp(scenix, subix, first(TuLiPa.getvarout(statevar)))
                            db.output.storagevalues[subix][stepnr, dim+1+get_numscen_stoch(db.input)+scenix, j] = fetch(f)
                        end
                    end
                end
            elseif haskey(settings["problems"], "stochastic")
                statevars = db.mp[subix].cuts.statevars
                for scenix in 1:get_numscen_stoch(db.input)
                    for (j, statevar) in enumerate(statevars)
                        core_stoch = get_core_sp(db, scenix, subix)
                        f = @spawnat core_stoch get_enddual_stoch(scenix, subix, first(TuLiPa.getvarout(statevar)))
                        db.output.storagevalues[subix][stepnr, dim+scenix, j] = fetch(f)
                    end
                end
            end
        end
    end

    if haskey(settings["results"], "mainresults")
        termduration = get_steplength(db.input)
        if get_onlysubsystemmodel(db.input)
            periodduration_power = parse_duration(settings["horizons"]["master"]["Power"], "periodduration")
            if haskey(settings["horizons"]["master"], "Hydro")
                periodduration_hydro = parse_duration(settings["horizons"]["master"]["Hydro"], "periodduration")
            else
                periodduration_hydro = periodduration_power
            end
            prob_results = db.mp[first(db.dist_mp[1])].prob
        else
            periodduration_power = parse_duration(settings["horizons"]["clearing"]["Power"], "periodduration")
            periodduration_hydro = parse_duration(settings["horizons"]["clearing"]["Hydro"], "periodduration")
            prob_results = db.cp.prob
        end
        numperiods_powerhorizon = Int(termduration.value / periodduration_power.value)
        numperiods_hydrohorizon = Int(termduration.value / periodduration_hydro.value)

        if stepnr == 1 # TODO: move to init
            db.output.modelobjects = Dict(zip([TuLiPa.getid(obj) for obj in TuLiPa.getobjects(prob_results)], TuLiPa.getobjects(prob_results)))
            if settings["results"]["mainresults"] == "all"
                resultobjects = TuLiPa.getobjects(prob_results) # collect results for all areas
            else
                resultobjects = TuLiPa.getpowerobjects(db.output.modelobjects, settings["results"]["mainresults"]); # only collect results for one area
            end

            powerbalances, rhsterms, rhstermbalances, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages = TuLiPa.order_result_objects(resultobjects, true)
            db.output.powerbalances = powerbalances
            db.output.rhsterms = rhsterms
            db.output.rhstermbalances = rhstermbalances
            db.output.plants = plants
            db.output.plantbalances = plantbalances
            db.output.plantarrows = plantarrows
            db.output.demands = demands
            db.output.demandbalances = demandbalances
            db.output.demandarrows = demandarrows
            db.output.hydrostorages = hydrostorages
            db.output.batterystorages = batterystorages

            db.output.prices = zeros(Int(numperiods_powerhorizon*steps), length(db.output.powerbalances))
            db.output.rhstermvalues = zeros(Int(numperiods_powerhorizon*steps), length(db.output.rhsterms))
            db.output.production = zeros(Int(numperiods_powerhorizon*steps), length(db.output.plants))
            db.output.consumption = zeros(Int(numperiods_powerhorizon*steps), length(db.output.demands))
            db.output.hydrolevels = zeros(Int(numperiods_hydrohorizon*steps), length(db.output.hydrostorages))
            db.output.batterylevels = zeros(Int(numperiods_powerhorizon*steps), length(db.output.batterystorages))

            if haskey(settings["results"], "otherterms")
                otherinfo = settings["results"]["otherterms"]
                otherobjects, otherbalances = TuLiPa.order_result_objects_other(resultobjects, otherinfo)
                db.output.otherobjects = otherobjects
                db.output.otherbalances = otherbalances
    
                for key in keys(otherinfo)
                    db.output.othervalues[key] = Dict()
            
                    for commodity in keys(otherinfo[key])
                        horizon = TuLiPa.get_horizon_commodity(resultobjects, commodity)
                        if key == "RHSTerms"
                            db.output.othervalues[key][commodity] = zeros(TuLiPa.getnumperiods(horizon)*steps, length(otherobjects[key][commodity]))
                        elseif key == "Vars"
                            db.output.othervalues[key][commodity] = zeros(TuLiPa.getnumperiods(horizon)*steps, length(otherobjects[key][commodity]))
                        end
                    end
                end
            end
        end

        if stepnr == 2 
            db.output.statenames = collect(keys(db.startstates))
            db.output.statematrix = zeros(length(values(db.startstates)), Int(steps))
        end
        if stepnr != 1
            db.output.statematrix[:,stepnr-1] .= collect(values(db.startstates))
        end

        powerrange = Int(numperiods_powerhorizon*(stepnr-1)+1):Int(numperiods_powerhorizon*(stepnr))
        hydrorange = Int(numperiods_hydrohorizon*(stepnr-1)+1):Int(numperiods_hydrohorizon*(stepnr))
        TuLiPa.get_results!(prob_results, db.output.prices, db.output.rhstermvalues, db.output.production, db.output.consumption, db.output.hydrolevels, db.output.batterylevels, db.output.powerbalances, db.output.rhsterms, db.output.plants, db.output.plantbalances, db.output.plantarrows, db.output.demands, db.output.demandbalances, db.output.demandarrows, db.output.hydrostorages, db.output.batterystorages, db.output.modelobjects, powerrange, hydrorange, periodduration_power, t)

        if haskey(settings["results"], "otherterms")
            TuLiPa.get_results!(stepnr, prob_results, db.output.otherobjects, db.output.otherbalances, db.output.othervalues, db.output.modelobjects, t)
        end
    end

    collect_interval = get_result_prices_ppp(settings)
    if collect_interval != 0
        if (stepnr-1) % collect_interval == 0
            collect_step = div((stepnr-1), collect_interval) + 1

            futures = []
            stoch_scenixs = [scenario.parentscenario for scenario in db.scenmod_stoch.scenarios]
            @sync for (scenix, core) in get_dist_ppp(db)
                if scenix in stoch_scenixs
                    f = @spawnat core get_ppp_prices(scenix, db.output.prices_balances)
                    push!(futures, f)
                end
            end

            for (i, f) in enumerate(futures)
                pl, dl, pm, dm, ps, ds = fetch(f)

                db.output.prices_long[collect_step, :, i, :] .= pl
                db.output.deltas_long[collect_step, :] .= dl

                db.output.prices_med[collect_step, :, i, :] .= pm
                db.output.deltas_med[collect_step, :] .= dm

                db.output.prices_short[collect_step, :, i, :] .= ps
                db.output.deltas_short[collect_step, :] .= ds
            end

            @sync for (scenix, core) in get_dist_ppp(db)
                @spawnat core reset_ppp_prices(scenix)
            end
        end
    end
end

function get_enddual_stoch(scenix, subix, objid)
    db = get_local_db()
    sp = db.sp[(scenix, subix)]

    obj = get_obj_from_id(TuLiPa.getobjects(sp.prob), objid) # TODO: OK to assume objid = varoutid?
    balance = TuLiPa.getbalance(obj)
    dual = TuLiPa.getcondual(sp.prob, TuLiPa.getid(balance), TuLiPa.getnumperiods(TuLiPa.gethorizon(balance)))
    if haskey(balance.metadata, TuLiPa.GLOBALENEQKEY)
        dual /= balance.metadata[TuLiPa.GLOBALENEQKEY]
    end

    return dual
end

function get_enddual_evp(scenix, subix, objid)
    db = get_local_db()
    evp = db.evp[(scenix, subix)]

    obj = get_obj_from_id(TuLiPa.getobjects(evp.prob), objid) # TODO: OK to assume objid = varoutid?
    balance = TuLiPa.getbalance(obj)
    dual = TuLiPa.getcondual(evp.prob, TuLiPa.getid(balance), TuLiPa.getnumperiods(TuLiPa.gethorizon(balance)))
    if haskey(balance.metadata, TuLiPa.GLOBALENEQKEY)
        dual /= balance.metadata[TuLiPa.GLOBALENEQKEY]
    end

    return dual
end

function reset_ppp_prices(scenix)
    db = get_local_db()
    ppp = db.ppp[scenix]
    fill!(ppp.div["prices_long"], 0.0)
    fill!(ppp.div["deltas_long"], 0.0)
    fill!(ppp.div["prices_med"], 0.0)
    fill!(ppp.div["deltas_med"], 0.0)
    fill!(ppp.div["prices_short"], 0.0)
    fill!(ppp.div["deltas_short"], 0.0)
    return
end

function get_ppp_prices(scenix, bids)
    db = get_local_db()
    ppp = db.ppp[scenix]

    balance_long = get_obj_from_id(TuLiPa.getobjects(ppp.longprob), bids[1])
    horizon_long = TuLiPa.gethorizon(balance_long)
    numperiods_long = TuLiPa.getnumperiods(horizon_long)

    delta_long = 0
    for t in 1:numperiods_long
        delta_long += TuLiPa.getduration(TuLiPa.gettimedelta(horizon_long, t)).value
        ppp.div["deltas_long"][t] = delta_long
        for (i, bid) in enumerate(bids)
            ppp.div["prices_long"][i,t]  = TuLiPa.getcondual(ppp.longprob, bid, t)
        end
    end

    balance_med = get_obj_from_id(TuLiPa.getobjects(ppp.medprob), bids[1])
    horizon_med = TuLiPa.gethorizon(balance_med)
    numperiods_med = TuLiPa.getnumperiods(horizon_med)

    delta_med = 0
    for t in 1:numperiods_med
        delta_med += TuLiPa.getduration(TuLiPa.gettimedelta(horizon_med, t)).value
        ppp.div["deltas_med"][t] = delta_med
        for (i, bid) in enumerate(bids)
            ppp.div["prices_med"][i,t]  = TuLiPa.getcondual(ppp.medprob, bid, t)
        end
    end

    balance_short = get_obj_from_id(TuLiPa.getobjects(ppp.shortprob), bids[1])
    horizon_short = TuLiPa.gethorizon(balance_short)
    numperiods_short = TuLiPa.getnumperiods(horizon_short)

    delta_short = 0
    for t in 1:numperiods_short
        delta_short += TuLiPa.getduration(TuLiPa.gettimedelta(horizon_short, t)).value
        ppp.div["deltas_short"][t] = delta_short
        for (i, bid) in enumerate(bids)
            ppp.div["prices_short"][i,t]  = TuLiPa.getcondual(ppp.shortprob, bid, t)
        end
    end

    return ppp.div["prices_long"], ppp.div["deltas_long"], ppp.div["prices_med"], ppp.div["deltas_med"], ppp.div["prices_short"], ppp.div["deltas_short"]
end

get_output_from_input(input::DefaultJulESInput) = DefaultJulESOutput(input)

get_maintiming_ppp(scenix) = get_local_db().ppp[scenix].div[MainTiming]
get_maintiming_evp(scenix, subix) = get_local_db().evp[(scenix, subix)].div[MainTiming]
get_maintiming_mp(subix) = get_local_db().mp[subix].div[MainTiming]
get_maintiming_sp(scenix, subix) = get_local_db().sp[(scenix, subix)].div[MainTiming]

reset_maintiming_ppp(scenix) = fill!(get_local_db().ppp[scenix].div[MainTiming], 0.0)
reset_maintiming_evp(scenix, subix) = fill!(get_local_db().evp[(scenix, subix)].div[MainTiming], 0.0)
reset_maintiming_mp(subix) = fill!(get_local_db().mp[subix].div[MainTiming], 0.0)
reset_maintiming_sp(scenix, subix) = fill!(get_local_db().sp[(scenix, subix)].div[MainTiming], 0.0)

get_storagevalues_stoch(subix) = get_local_db().mp[subix].div[StorageValues]

function get_output_final(steplength, skipmax)
    output = get_output_main()

    get_output_timing(output, steplength, skipmax)

    get_output_scenarios(output)

    get_output_storagevalues(output, steplength, skipmax)

    get_output_ppp_prices(output)

    get_output_memory(output) # TODO: Find problem

    return output
end

function get_output_storagevalues(output, steplength, skipmax)
    db = get_local_db()
    settings = get_settings(db)
    
    if has_result_storagevalues(settings)
        f = @spawnat db.core_main get_output_storagevalues_local(output, steplength, skipmax)
        storagenames, storagevalues, shorts, scenarionames, skipfactor = fetch(f)
        
        output[StorageValues] = cat(storagevalues..., dims=3)
        output["storagenames"] = storagenames
        output["shorts"] = shorts
        output["scenarionames"] = scenarionames
        output["skipfactor"] = skipfactor
    end
end

function get_output_storagevalues_local(output, steplength, skipmax)
    db = get_local_db()
    settings = get_settings(db)

    storagevalues = []
    storagenames = String[]
    shorts = Bool[]
    for (subix, core) in get_dist_mp(db)
        push!(storagevalues, db.output.storagevalues[subix])
        substoragenames = fetch(@spawnat core get_storagenames_from_subix(subix))
        storagenames = vcat(storagenames, substoragenames)
        short = !get_skipmed_impact(db.subsystems[subix])
        for substoragename in substoragenames
            push!(shorts, short)
        end
    end

    scenarionames = String[]
    for i in 1:get_numscen_stoch(db.input)
        push!(scenarionames, string(i) * " min")
        push!(scenarionames, string(i) * " max")
    end
    push!(scenarionames, "Operative master")
    if has_headlosscost(settings["problems"]["stochastic"]["master"])
        push!(scenarionames, "Operative master after")
    end
    if haskey(settings["problems"], "clearing")
        push!(scenarionames, "Operative clearing")
    end
    if haskey(settings["problems"], "stochastic")
        for i in 1:get_numscen_stoch(db.input)
            push!(scenarionames, string(i) * " stochend")
        end
    end
    if haskey(settings["problems"], "endvalue")
        for i in 1:get_numscen_stoch(db.input)
            push!(scenarionames, string(i) * " evpend")
        end
    end

    skipfactor = (skipmax+Millisecond(steplength))/Millisecond(steplength)

    return (storagenames, storagevalues, shorts, scenarionames, skipfactor)
end

function get_storagenames_from_subix(subix)
    db = get_local_db()
    
    storagenames = String[]
    for (j, statevar) in enumerate(db.mp[subix].cuts.statevars)
        push!(storagenames, TuLiPa.getinstancename(first(TuLiPa.getvarout(statevar))))
    end
    return storagenames
end

function get_output_ppp_prices(output)
    db = get_local_db()
    settings = get_settings(db)
    
    if get_result_prices_ppp(settings) != 0
        f = @spawnat db.core_main get_output_prices_ppp_local()
        ret = fetch(f)
        if ret isa RemoteException
            throw(ret)
        end
        balancenames, pl, dl, pm, dm, ps, ds = ret
        
        output["balancenames_ppp"] = balancenames
        output["prices_long"] = pl
        output["deltas_long"] = dl
        output["prices_med"] = pm
        output["deltas_med"] = dm
        output["prices_short"] = ps
        output["deltas_short"] = ds
    end
end

function get_output_prices_ppp_local()
    db = get_local_db()

    balancenames = [TuLiPa.getinstancename(bid) for bid in db.output.prices_balances]
    pl = db.output.prices_long
    dl = db.output.deltas_long
    pm = db.output.prices_med
    dm = db.output.deltas_med
    ps = db.output.prices_short
    ds = db.output.deltas_short

    return balancenames, pl, dl, pm, dm, ps, ds
end

function get_output_memory(output)
    db = get_local_db()
    settings = get_settings(db)

    if has_result_memory(settings)
        names = ["coreid", "sum_unique", "core", "input", "output", "horizons", "dummyobjects", "dummyobjects_ppp", "startstates", "subsystems", "subsystems_evp", "subsystems_stoch", "scenmod_sim", "scenmod_stoch", "ifm", "ppp", "prices_ppp", "evp", "mp", "sp", "cp", "dist_ifm", "dist_ppp", "dist_evp", "dist_mp", "dist_sp", "core_main", "div", "ifm_output", "ifm_derived"]
        df = DataFrame(DataFrame([[] for _ = names] , names))
        cores = get_cores(db)
        for core in cores # TODO: Do sync
            f = @spawnat core get_output_memory_local()
            push!(df, fetch(f))
        end
        df[!, :sum] = sum(eachcol(select(df, Not(:coreid, :sum_unique))))
        df = permutedims(df, "coreid")
        if length(cores) > 1
            df[!, :sum] = sum(eachcol(select(df, Not(:coreid))))
        end
        println(df)
    end
end

function get_output_memory_local()
    db = get_local_db()

    values = Any[string(db.core), Base.summarysize(db)/1e6]

    for field in fieldnames(typeof(db))
        field_value = getfield(db, field)
        field_memory_size = Base.summarysize(field_value)/1e6
        push!(values, field_memory_size)
    end
    return values
end

function get_output_scenarios(output)
    db = get_local_db()

    wait(@spawnat db.core_main get_output_scenarios_local(output))
end

function get_output_scenarios_local(data)
    db = get_local_db()

    data["scenweights_sim"] = db.output.scenweights_sim
    data["scenweights_stoch"] = db.output.scenweights_stoch
end

function get_output_timing(output, steplength, skipmax)
    db = get_local_db()

    wait(@spawnat db.core_main get_output_timing_local(output, steplength, skipmax))
end

function get_output_timing_local(data, steplength, skipmax)
    db = get_local_db()
    settings = get_settings(db)

    if has_result_times(settings)
        skipfactor = (skipmax+Millisecond(steplength))/Millisecond(steplength)
        
        timing_cp = get_timing_cp_local()

        timings_ppp = []
        for (scenix, values) in db.output.timing_ppp
            push!(timings_ppp, values)
        end

        # TODO: Add subix name
        df_evp = DataFrame([name => [] for name in ["scenix", "subix", "update", "solve", "total", "core", "skipmed"]])
        for (scenix, subix, core) in db.dist_evp
            values = dropdims(mean(db.output.timing_evp[(scenix, subix)], dims=1), dims=1)
            f = @spawnat core get_skipmed_impact(subix)
            push!(df_evp, [scenix, subix, values[1], values[2], values[3], core, fetch(f)])
        end
        df_evp[!, :other] = df_evp[!, :total] - df_evp[!, :solve] - df_evp[!, :update]
        df_evp[df_evp.skipmed .== true, [:update, :solve, :total]] .= df_evp[df_evp.skipmed .== true, [:update, :solve, :total]] .* skipfactor
        if nrow(df_evp) != 0
            df_evp_subix = combine(groupby(df_evp, [:subix]), 
            :update => sum => :evp_u, 
            :solve => sum => :evp_s, 
            :other => sum => :evp_o,
            :total => sum => :evp_tot)
            timings_evp = mean.(eachcol(select(df_evp_subix, Not([:subix, :evp_o]))))
            df_evp_core = combine(groupby(df_evp, [:core]), 
            :update => sum => :evp_u, 
            :solve => sum => :evp_s, 
            :other => sum => :evp_o,
            :total => sum => :evp_tot)
        else
            df_evp_subix = DataFrame([name => [] for name in ["subix", "evp_u", "evp_s", "evp_o", "evp_tot"]])
            timings_evp = [0.0, 0.0, 0.0]
            df_evp_core = DataFrame([name => [] for name in ["core", "evp_u", "evp_s", "evp_o", "evp_tot"]])
        end
        # TODO: df_evp_scen

        df_mp = DataFrame([name => [] for name in ["subix", "mp_u", "mp_s", "mp_fin", "mp_o", "bend_it", "core", "skipmed"]])
        for (subix, core) in db.dist_mp
            values = dropdims(mean(db.output.timing_mp[(subix)], dims=1), dims=1)
            f = @spawnat core get_skipmed_impact(subix)
            push!(df_mp, [subix, values[1], values[2], values[3], values[4], values[5], core, fetch(f)])
        end
        df_mp[!, :mp_tot] = df_mp[!, :mp_s] + df_mp[!, :mp_u] + df_mp[!, :mp_fin] + df_mp[!, :mp_o]
        df_mp[df_mp.skipmed .== true, [:mp_u, :mp_s, :mp_fin, :mp_o, :mp_tot, :bend_it]] .= df_mp[df_mp.skipmed .== true, [:mp_u, :mp_s, :mp_fin, :mp_o, :mp_tot, :bend_it]] .* skipfactor
        timings_mp = mean.(eachcol(select(df_mp, Not([:subix, :core, :skipmed, :mp_fin, :mp_o, :bend_it]))))
        df_mp_core = combine(groupby(df_mp, [:core]), 
        :mp_u => sum => :mp_u, 
        :mp_s => sum => :mp_s, 
        :mp_fin => sum => :mp_fin,
        :mp_o => sum => :mp_o,
        :mp_tot => sum => :mp_tot)

        df_sp = DataFrame([name => [] for name in ["scenix", "subix", "update", "solve", "other", "core", "skipmed"]])
        for (scenix, subix, core) in db.dist_sp
            values = dropdims(mean(db.output.timing_sp[(scenix, subix)], dims=1), dims=1)
            f = @spawnat core get_skipmed_impact(subix)
            push!(df_sp, [scenix, subix, values[1], values[2], values[3], core, fetch(f)])
        end
        df_sp[!, :total] = df_sp[!, :solve] + df_sp[!, :update] + df_sp[!, :other]
        df_sp[df_sp.skipmed .== true, [:update, :solve, :other, :total]] .= df_sp[df_sp.skipmed .== true, [:update, :solve, :other, :total]] .* skipfactor
        df_sp_subix = combine(groupby(df_sp, [:subix]), 
        :update => sum => :sp_u, 
        :solve => sum => :sp_s, 
        :other => sum => :sp_o,
        :total => sum => :sp_tot)
        timings_sp = mean.(eachcol(select(df_sp_subix, Not([:subix, :sp_o]))))
        df_sp_core = combine(groupby(df_sp, [:core]), 
        :update => sum => :sp_u, 
        :solve => sum => :sp_s, 
        :other => sum => :sp_o,
        :total => sum => :sp_tot)
        # TODO: df_sp_scen

        df_subix = outerjoin(df_evp_subix, df_mp, df_sp_subix, on = :subix)
        df_subix = coalesce.(df_subix, 0.0)
        df_subix[!, :tot] = df_subix[!, :evp_tot] + df_subix[!, :mp_tot] + df_subix[!, :sp_tot]
        df_subix = sort(df_subix, :tot, rev=true)
        df_subix = df_subix[!, [:subix, :tot, :evp_tot, :mp_tot, :sp_tot, :evp_u, :evp_s, :evp_o, :bend_it, :mp_u, :mp_s, :mp_fin, :mp_o, :sp_u, :sp_s, :sp_o]]

        df_core = outerjoin(df_evp_core, df_mp_core, df_sp_core, on = :core)
        df_core = coalesce.(df_core, 0.0)
        df_core[!, :tot] = df_core[!, :evp_tot] + df_core[!, :mp_tot] + df_core[!, :sp_tot]
        df_core = sort(df_core, :tot, rev=true)
        df_core = df_core[!, [:core, :tot, :evp_tot, :mp_tot, :sp_tot, :evp_u, :evp_s, :evp_o, :mp_u, :mp_s, :mp_fin, :mp_o, :sp_u, :sp_s, :sp_o]]

        if haskey(settings["problems"], "prognosis") && haskey(settings["problems"], "clearing")
            factors = [skipfactor,skipfactor,1]
            dims = size(timings_ppp[1])
            dims = (dims..., length(timings_ppp))
            timings_ppp1 = reshape(cat(timings_ppp..., dims=4), dims)
            timings_ppp2 = transpose(dropdims(mean(timings_ppp1,dims=(1,4)),dims=(1,4))).*factors
            all = vcat(timings_ppp2, reshape(timings_evp,1,3), reshape(timings_mp,1,3), reshape(timings_sp,1,3), mean(timing_cp, dims=1))
            df = DataFrame(model=["long","med","short","evp","mp","sp","clearing"], update=all[:,1], solve=all[:,2], total=all[:,3])
            df[!, :other] = df[!, :total] - df[!, :solve] - df[!, :update]
            display(df[!, [1, 2, 3, 5, 4]])
        end

        display(df_core)
        display(df_subix)
    end
    # # display number of elements of each type per subsystem
    # unique_types = ["subix", "name_first_element", "name_second_element", "total_count"] #  
    # # loop gjennom get_elements(db.input) og legg til push!(unique_type, )
    # df_sub_element_type = DataFrame([name => [] for name in unique_type])
    # for (subix, subsystem) in enumerate(get_subsystems(db))
    #     type_count = []
    #     total_count = # sum av type_count
    #     name_first_element = # element.instancename
    #     name_second_element = 
    #     push!(df_sub_element_type, [subix, name_first_element, name_second_element, total_count, type_count...])
    # end
    # # sorter df_sub_element_type etter total_count
    # display(df_sub_element_type)

    if settings["results"]["times"]
        if haskey(settings["problems"], "prognosis") 
            data["prognosistimes"] = timings_ppp1
        end
        if haskey(settings["problems"], "endvalue") 
            data["endvaluetimes"] = Matrix{Float64}(df_evp)
        end
        if haskey(settings["problems"], "stochastic") 
            data["mptimes"] = Matrix{Float64}(df_mp)
            data["sptimes"] = Matrix{Float64}(df_sp)
        end
        if haskey(settings["problems"], "clearing")
            data["clearingtimes"] = timing_cp
        end
    end
end

get_skipmed_impact(subix) = get_skipmed_impact(get_local_db().subsystems[subix])

function get_output_main()
    db = get_local_db()
    future = @spawnat db.core_main get_output_main_local()
    ret = fetch(future)
    if ret isa RemoteException
        throw(ret)
    end
    return ret
end

function get_timing_cp_local()
    db = get_local_db()
    return db.output.timing_cp
end

function get_output_main_local()
    db = get_local_db()
    settings = get_settings(db)
    mainconfig = get_mainconfig(db)

    data = Dict()

    # Final update of statevariables
    if get_onlysubsystemmodel(db.input)
        startstates_main = get_startstates_from_mp()
    else
        startstates_main =  get_startstates_from_cp()
    end
    for (k, v) in startstates_main
        db.startstates[k] = v
    end
    steps = get_steps(db)
    if steps == 1
        db.output.statenames = collect(keys(db.startstates))
        db.output.statematrix = collect(values(db.startstates))
    else
        db.output.statematrix[:,steps] .= collect(values(db.startstates))
    end

    if haskey(settings["results"], "mainresults")
        steplength = get_steplength(db.input)
        if get_onlysubsystemmodel(db.input)
            periodduration_power = parse_duration(settings["horizons"]["master"]["Power"], "periodduration")
            if haskey(settings["horizons"]["master"], "Hydro")
                periodduration_hydro = parse_duration(settings["horizons"]["master"]["Hydro"], "periodduration")
            else
                periodduration_hydro = periodduration_power
            end
        else
            periodduration_power = parse_duration(settings["horizons"]["clearing"]["Power"], "periodduration")
            periodduration_hydro = parse_duration(settings["horizons"]["clearing"]["Hydro"], "periodduration")
        end

        # Only keep rhsterms that have at least one value (TODO: Do the same for sypply and demands)
        rhstermtotals = dropdims(sum(db.output.rhstermvalues,dims=1),dims=1)
        rhstermsupplyidx = []
        rhstermdemandidx = []

        for k in 1:length(db.output.rhsterms)
            if rhstermtotals[k] > 0
                push!(rhstermsupplyidx, k)
            elseif rhstermtotals[k] < 0
                push!(rhstermdemandidx, k)
            end
        end

        # Put rhsterms together with supplies and demands
        rhstermsupplyvalues = db.output.rhstermvalues[:,rhstermsupplyidx]
        rhstermdemandvalues = db.output.rhstermvalues[:,rhstermdemandidx]*-1

        rhstermsupplynames = [TuLiPa.getinstancename(rhsterm) for rhsterm in db.output.rhsterms[rhstermsupplyidx]]
        rhstermsupplybalancenames = [split(TuLiPa.getinstancename(r), "PowerBalance_")[2] for r in db.output.rhstermbalances[rhstermsupplyidx]]
        rhstermdemandnames = [TuLiPa.getinstancename(rhsterm) for rhsterm in db.output.rhsterms[rhstermdemandidx]]
        rhstermdemandbalancenames = [split(TuLiPa.getinstancename(r), "PowerBalance_")[2] for r in db.output.rhstermbalances[rhstermdemandidx]]

        supplynames = [[TuLiPa.getinstancename(plant) for plant in db.output.plants];rhstermsupplynames]
        supplybalancenames = [[split(TuLiPa.getinstancename(p), "PowerBalance_")[2] for p in db.output.plantbalances];rhstermsupplybalancenames]
        supplyvalues = hcat(db.output.production,rhstermsupplyvalues)

        demandnames = [[TuLiPa.getinstancename(demand) for demand in db.output.demands];rhstermdemandnames]
        demandbalancenames = [[split(TuLiPa.getinstancename(p), "PowerBalance_")[2] for p in db.output.demandbalances];rhstermdemandbalancenames]
        demandvalues = hcat(db.output.consumption, rhstermdemandvalues)

        # Prepare for plotting results
        hydronames = [TuLiPa.getinstancename(hydro) for hydro in db.output.hydrostorages]
        batterynames = [TuLiPa.getinstancename(battery) for battery in db.output.batterystorages]
        powerbalancenames = [split(TuLiPa.getinstancename(TuLiPa.getid(powerbalance)), "PowerBalance_")[2] for powerbalance in db.output.powerbalances]

        # Convert reservoir filling to TWh
        hydrolevels1 = copy(db.output.hydrolevels)
        for (i,hydroname) in enumerate(hydronames)
            if haskey(TuLiPa.getbalance(db.output.modelobjects[db.output.hydrostorages[i]]).metadata, TuLiPa.GLOBALENEQKEY)
                hydrolevels1[:,i] .= hydrolevels1[:,i]*TuLiPa.getbalance(db.output.modelobjects[db.output.hydrostorages[i]]).metadata[TuLiPa.GLOBALENEQKEY]
            end
        end

        # Indexes TODO: Replace with generic (for each commodity, and for state)
        dim = get_outputindex(mainconfig, get_datayear(db), get_weatheryear(db))
        x1 = [TuLiPa.getisoyearstart(dim) + Week(mainconfig["weekstart"]-1) + periodduration_power*(t-1) for t in 1:first(size(supplyvalues))] # power/load resolution
        x2 = [TuLiPa.getisoyearstart(dim) + Week(mainconfig["weekstart"]-1) + periodduration_hydro*(t-1) for t in 1:first(size(hydrolevels1))]; # reservoir resolution
        x3 = [TuLiPa.getisoyearstart(dim) + Week(mainconfig["weekstart"]-1) + steplength*(t-1) for t in 1:steps]; # state resolution

        outputformat = mainconfig["outputformat"]
        if outputformat != "juliadict"
            datetimeformat = mainconfig["datetimeformat"]
            x1 = Dates.format.(x1, datetimeformat)
            x2 = Dates.format.(x2, datetimeformat)
            x3 = Dates.format.(x3, datetimeformat)
        end

        data["areanames"] = powerbalancenames |> Vector{String}
        data["pricematrix"] = db.output.prices
        data["priceindex"] = x1

        data["resnames"] = hydronames
        data["resmatrix"] = hydrolevels1
        data["resindex"] =  x2
        if has_result_hydrolevels_water(settings)
            data["resmatrix_water"] = db.output.hydrolevels
        end

        data["batnames"] = batterynames
        data["batmatrix"] = db.output.batterylevels
        data["batindex"] =  x2

        data["statenames"] = db.output.statenames
        data["statematrix"] = permutedims(db.output.statematrix)
        data["stateindex"] =  x3

        data["supplyvalues"] = supplyvalues
        data["supplynames"] = supplynames
        data["supplybalancenames"] = supplybalancenames

        data["demandvalues"] = demandvalues
        data["demandnames"] = demandnames
        data["demandbalancenames"] = demandbalancenames

        if haskey(settings["results"], "otherterms")
            for key in keys(db.output.othervalues)
                for commodity in keys(db.output.othervalues[key])
                    data["othernames_" * key * "_" * commodity] = [TuLiPa.getinstancename(id) for id in db.output.otherobjects[key][commodity]] |> Vector{String}
                    data["othervalues_" * key * "_" * commodity] = db.output.othervalues[key][commodity]
                end
            end
        end

        if has_ifm_results(db.input)
            data["ifm_names"] = db.output.ifm_stations
            data["ifm_index"] =  x3
            for stateix in eachindex(db.output.ifm_u0)
                statename = db.output.ifm_statenames[stateix]
                data["ifm_startstates_$(statename)_$(stateix)"] = db.output.ifm_u0[stateix]
            end
            data["ifm_steamflow"] =  db.output.ifm_Q
        end
    end

    return data
end