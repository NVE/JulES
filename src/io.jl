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

get_iprogtype(input::DefaultJulESInput) = input.dataset["iprogtype"]
get_ifm_normfactors(input::DefaultJulESInput) = input.dataset["ifm_normfactors"]
get_ifm_elements(input::DefaultJulESInput) = input.dataset["ifm_elements"]

function get_ifm_names(input::DefaultJulESInput)
    s1 = Set(input.dataset["ifm_names"])
    s2 = Set([e.instancename for e in input.dataset["ifm_elements"] if e.conceptname == ABSTRACT_INFLOW_MODEL])
    return String[i for i in intersect(s1, s2)]
end

function get_ifm_replacemap(input::DefaultJulESInput) 
    names = Set(get_ifm_names(input))
    aggnames = Set(keys(get_ifm_weights(input)))
    d = Dict{String, String}()
    for (k, v) in input.dataset["ifm_replacemap"]
        if (v in names) || (v in aggnames)
            d[k] = v
        end
    end
    return d
end

function get_ifm_weights(input::DefaultJulESInput)
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
    steplength = parse_duration(settings["horizons"]["clearing"], "termduration")
    steps = Int(ceil((TuLiPa.getisoyearstart(datayear + simulationyears) - TuLiPa.getisoyearstart(datayear)).value/steplength.value) + extrasteps);
    
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
    elseif timetype == "PrognosisTime"
        return TuLiPa.PrognosisTime(datasimtime, datasimtime, weatherscenariotime)
    elseif timetype == "FixedDataTwoTime"
        return TuLiPa.FixedDataTwoTime(datasimtime, weatherscenariotime)
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
    elseif s == "KMeansAHMethod()"
        return TuLiPa.KMeansAHMethod()
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
function get_numscen_stoch(input::AbstractJulESInput)
    settings = get_settings(input)
    if haskey(settings["scenariogeneration"], "stochastic")
        return settings["scenariogeneration"]["stochastic"]["numscen"]
    else
        return get_numscen_sim(input)
    end
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
        return TuLiPa.DynamicExogenPriceAHData(TuLiPa.Id("Balance", rhsdata["balance"])) # TODO: If dynamic use tphasein
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
        minperiod = parse_duration(settings["horizons"]["clearing"], "termduration")
        horizon = TuLiPa.ShrinkableHorizon(horizon, startafter, shrinkatleast, minperiod)
    end

    return horizon
end

function build_adaptivehorizon(term, commodity, settings, n_durations, datayear)
    weatheryearstart = settings["time"]["weatheryearstart"]
    weatheryearstop = settings["time"]["weatheryearstop"]
    rhsdata = getrhsdata(settings["horizons"][term][commodity]["rhsdata"], datayear, weatheryearstart, weatheryearstop)
    rhsmethod = parse_methods(settings["horizons"][term][commodity]["rhsmethod"])
    clusters = settings["horizons"][term][commodity]["clusters"]
    unitduration = Millisecond(Hour(settings["horizons"][term][commodity]["unitduration_hours"]))

    horizon = TuLiPa.AdaptiveHorizon(clusters, unitduration, rhsdata, rhsmethod, n_durations...)

    if get_shrinkable(settings["horizons"][term])
        startafter = parse_duration(settings["horizons"]["shrinkable"], "startafter")
        shrinkatleast = parse_duration(settings["horizons"]["shrinkable"], "shrinkatleast")
        minperiod = parse_duration(settings["horizons"]["clearing"], "termduration")
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

    # TODO: info on scenix -> scenario for each step

    prices::Array{Float64}
    rhstermvalues::Array{Float64}
    production::Array{Float64}
    consumption::Array{Float64}
    hydrolevels::Array{Float64}
    batterylevels::Array{Float64}
    
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

    statenames::Vector{String}
    statematrix::Array{Float64} # end states after each step

    function DefaultJulESOutput(input)
        return new(Dict(),Dict(),Dict(),Dict(),[],
        Dict(),
        [],[],[],[],[],[],
        Dict(),[],[],[],[],[],Dict(),[],[],Dict(),[],[],[],[])
    end
end

get_cores(output::DefaultJulESOutput) = output.cores

function init_local_output()
    db = get_local_db()
    settings = get_settings(db)
    db.output = get_output_from_input(db.input)

    steps = get_steps(db)

    for (scenix, core) in db.dist_ppp
        db.output.timing_ppp[scenix] = zeros(steps, 3, 3) # TODO: Matrix more flexible long term, and array instead of dict?
    end

    for (scenix, subix, core) in db.dist_evp
        db.output.timing_evp[(scenix, subix)] = zeros(steps, 3)
    end

    for (subix, core) in db.dist_mp
        db.output.timing_mp[subix] = zeros(steps, 4)
        if settings["results"]["storagevalues"]
            if get_headlosscost(settings["problems"]["stochastic"]["master"])
                num_storagevalues = get_numscen_stoch(db.input)*2 + 2 # scenarios + master operative + master operative after headlosscost adjustment
            else
                num_storagevalues = get_numscen_stoch(db.input)*2 + 1 # scenarios + master operative 
            end
            if haskey(settings["problems"], "clearing")
                num_storagevalues += 1
            end
            f = @spawnat core get_numstates(subix)
            db.output.storagevalues[subix] = zeros(steps, num_storagevalues, fetch(f))
        end
    end

    for (scenix, subix, core) in db.dist_sp
        db.output.timing_sp[(scenix, subix)] = zeros(steps, 3)
    end

    db.output.timing_cp = zeros(steps, 3)

end

function get_numstates(subix)
    db = get_local_db()
    return length(db.mp[subix].states)
end

function update_output(t::TuLiPa.ProbTime, stepnr::Int)
    db = get_local_db()
    settings = get_settings(db)
    steps = get_steps(db)

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

        if settings["results"]["storagevalues"]
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
                end
            end
        end
    end

    for (scenix, subix, core) in db.dist_sp
        f = @spawnat core get_maintiming_sp(scenix, subix)
        db.output.timing_sp[(scenix, subix)][stepnr, :] .= fetch(f)
        @spawnat core reset_maintiming_sp(scenix, subix)
    end

    if haskey(settings["results"], "mainresults")
        termduration = parse_duration(settings["horizons"]["clearing"], "termduration")
        periodduration_power = parse_duration(settings["horizons"]["clearing"]["Power"], "periodduration")
        periodduration_hydro = parse_duration(settings["horizons"]["clearing"]["Hydro"], "periodduration")
        numperiods_powerhorizon = Int(termduration.value / periodduration_power.value)
        numperiods_hydrohorizon = Int(termduration.value / periodduration_hydro.value)

        if stepnr == 1
            db.output.modelobjects = Dict(zip([TuLiPa.getid(obj) for obj in TuLiPa.getobjects(db.cp.prob)], TuLiPa.getobjects(db.cp.prob)))
            if settings["results"]["mainresults"] == "all"
                resultobjects = TuLiPa.getobjects(db.cp.prob) # collect results for all areas
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
        TuLiPa.get_results!(db.cp.prob, db.output.prices, db.output.rhstermvalues, db.output.production, db.output.consumption, db.output.hydrolevels, db.output.batterylevels, db.output.powerbalances, db.output.rhsterms, db.output.plants, db.output.plantbalances, db.output.plantarrows, db.output.demands, db.output.demandbalances, db.output.demandarrows, db.output.hydrostorages, db.output.batterystorages, db.output.modelobjects, powerrange, hydrorange, periodduration_power, t)
    end
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

    get_output_storagevalues(output, steplength, skipmax)

    get_output_memory(output) # TODO: Find problem

    return output
end

function get_output_storagevalues(output, steplength, skipmax)
    db = get_local_db()
    settings = get_settings(db)
    
    if settings["results"]["storagevalues"]
        f = @spawnat db.core_cp get_output_storagevalues_local(output, steplength, skipmax)
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
    storagenames = []
    shorts = []
    for (subix, core) in get_dist_mp(db)
        push!(storagevalues, db.output.storagevalues[subix])
        substoragenames = fetch(@spawnat core get_storagenames_from_subix(subix))
        storagenames = vcat(storagenames, substoragenames)
        short = !get_skipmed_impact(db.subsystems[subix])
        for substoragename in storagenames
            push!(shorts, short)
        end
    end

    scenarionames = String[]
    for i in 1:get_numscen_stoch(db.input)
        push!(scenarionames, string(i) * " min")
        push!(scenarionames, string(i) * " max")
    end
    push!(scenarionames, "Operative master")
    if get_headlosscost(settings["problems"]["stochastic"]["master"])
        push!(scenarionames, "Operative master after")
    end
    if haskey(settings["problems"], "clearing")
        push!(scenarionames, "Operative clearing")
    end

    skipfactor = (skipmax+Millisecond(steplength))/Millisecond(steplength)

    return (storagenames, storagevalues, shorts, scenarionames, skipfactor)
end

function get_storagenames_from_subix(subix)
    db = get_local_db()
    
    storagenames = []
    for (j, statevar) in enumerate(db.mp[subix].cuts.statevars)
        push!(storagenames, getinstancename(first(getvarout(statevar))))
    end
    return storagenames
end

function get_output_memory(output)
    db = get_local_db()
    settings = get_settings(db)

    if settings["results"]["memory"]
        names = ["coreid", "sum_unique", "core", "input", "output", "horizons", "dummyobjects", "dummyobjects_ppp", "startstates", "subsystems", "subsystems_evp", "subsystems_stoch", "scenmod_sim", "scenmod_stoch", "ppp", "prices_ppp", "evp", "mp", "sp", "cp", "dist_ppp", "dist_evp", "dist_mp", "dist_sp", "core_cp", "div"]
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

function get_output_timing(output, steplength, skipmax)
    db = get_local_db()

    wait(@spawnat db.core_cp get_output_timing_local(output, steplength, skipmax))
end

function get_output_timing_local(data, steplength, skipmax)
    db = get_local_db()
    settings = get_settings(db)

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
    # TODO: df_evp_scen

    df_mp = DataFrame([name => [] for name in ["subix", "mp_u", "mp_s", "mp_fin", "mp_o", "core", "skipmed"]])
    for (subix, core) in db.dist_mp
        values = dropdims(mean(db.output.timing_mp[(subix)], dims=1), dims=1)
        f = @spawnat core get_skipmed_impact(subix)
        push!(df_mp, [subix, values[1], values[2], values[3], values[4], core, fetch(f)])
    end
    df_mp[!, :mp_tot] = df_mp[!, :mp_s] + df_mp[!, :mp_u] + df_mp[!, :mp_fin] + df_mp[!, :mp_o]
    df_mp[df_mp.skipmed .== true, [:mp_u, :mp_s, :mp_fin, :mp_o, :mp_tot]] .= df_mp[df_mp.skipmed .== true, [:mp_u, :mp_s, :mp_fin, :mp_o, :mp_tot]] .* skipfactor
    timings_mp = mean.(eachcol(select(df_mp, Not([:subix, :core, :skipmed, :mp_fin, :mp_o]))))
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
    df_subix = df_subix[!, [:subix, :tot, :evp_tot, :mp_tot, :sp_tot, :evp_u, :evp_s, :evp_o, :mp_u, :mp_s, :mp_fin, :mp_o, :sp_u, :sp_s, :sp_o]]

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
    future = @spawnat db.core_cp get_output_cp_local()
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

function get_output_cp_local()
    db = get_local_db()
    settings = get_settings(db)
    mainconfig = get_mainconfig(db)

    data = Dict()

    # Final update of statevariables
    startstates_cp = get_startstates_from_cp()
    for (k, v) in startstates_cp
        db.startstates[k] = v
    end
    db.output.statematrix[:,get_steps(db)] .= collect(values(db.startstates))

    if haskey(settings["results"], "mainresults")
        steps = get_steps(db)
        steplength = parse_duration(settings["horizons"]["clearing"], "termduration")
        periodduration_power = parse_duration(settings["horizons"]["clearing"]["Power"], "periodduration")
        periodduration_hydro = parse_duration(settings["horizons"]["clearing"]["Hydro"], "periodduration")

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

        # Indexes
        dim = getoutputindex(mainconfig, get_datayear(db), get_weatheryear(db))
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
    end

    return data
end