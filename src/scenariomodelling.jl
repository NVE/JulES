struct NothingScenarioModellingMethod <: AbstractScenarioModellingMethod end
mutable struct NoScenarioModellingMethod{T <: AbstractScenario} <: AbstractScenarioModellingMethod
    scenarios::Vector{T}
end
# mutable struct ResidualLoadMethod <: ScenarioModellingMethod # choose scenario based on residual load (also energy inflow)
# end
mutable struct InflowClusteringMethod{T <: AbstractScenario} <: AbstractScenarioModellingMethod
    scenarios::Vector{T}
    inflowfactors::Vector{Float64}
    objects::Vector
    parts::Int
    scendelta::TuLiPa.MsTimeDelta

    function InflowClusteringMethod(scenarios, objects, parts, scendelta)
        inflowfactors = Vector{Float64}(undef, length(scenarios))
        return new{eltype(scenarios)}(scenarios, inflowfactors, objects, parts, scendelta)
    end
end
get_parts(method::InflowClusteringMethod{WeatherScenario}) = method.parts

mutable struct SumInflowQuantileMethod{T <: AbstractScenario} <: AbstractScenarioModellingMethod
    scenarios::Vector{T}
    inflowfactors::Vector{Float64}
    objects::Vector
    maxquantile::Float64 # parameter
    a::Float64 # parameter
    b::Float64 # parameter
    c::Float64 # parameter
    usedensity::Bool # parameter
    scendelta::TuLiPa.MsTimeDelta # parameter
    function SumInflowQuantileMethod(scenarios, objects, maxquantile, a, b, c, scendelta; usedensity=false)
        inflowfactors = Vector{Float64}(undef, length(scenarios))
        return new{eltype(scenarios)}(scenarios, inflowfactors, objects, maxquantile, a, b, c, scendelta, usedensity)
    end
end
get_maxquantile(method::SumInflowQuantileMethod) = method.maxquantile
get_a(method::SumInflowQuantileMethod) = method.a
get_b(method::SumInflowQuantileMethod) = method.b
get_c(method::SumInflowQuantileMethod) = method.c
get_usedensity(method::SumInflowQuantileMethod) = method.usedensity

get_scenarios(scenmod::NothingScenarioModellingMethod) = nothing
get_scenarios(scenmod::AbstractScenarioModellingMethod) = scenmod.scenarios

get_objects(method::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}) = method.objects
get_scendelta(method::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}) = method.scendelta

get_changes(scenmod::NoScenarioModellingMethod) = scenmod.scenarios
get_changes(scenmod::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}) = (scenmod.scenarios, scenmod.inflowfactors)

function set_changes(scenmod::NoScenarioModellingMethod, changes::Vector{WeatherScenario})
    scenmod.scenarios = changes
end
function set_changes(scenmod::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}, changes::Tuple{Vector{WeatherScenario}, Vector{Float64}})
    scenarios, inflowfactors = changes

    scenmod.scenarios = scenarios
    scenmod.inflowfactors = inflowfactors
end

get_inflowfactors(scenmod::AbstractScenarioModellingMethod) = [1/length(scenmod.scenarios) for s in 1:length(scenmod.scenarios)]
get_inflowfactors(scenmod::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}) = scenmod.inflowfactors

function choose_scenarios!(scenmod::SumInflowQuantileMethod{WeatherScenario}, scenmodmethodoptions::AbstractScenarioModellingMethod, simtime::TuLiPa.ProbTime, input::AbstractJulESInput)
    numscen = length(scenmod.scenarios)
    scenariooptions = get_scenarios(scenmodmethodoptions)
    weightsoptions = [get_probability(scenario) for scenario in scenariooptions]
    factoroptions = get_inflowfactors(scenmodmethodoptions)
    
    # Calculate total energy inflow in the system for the scenariodelta
    totalsumenergyinflow = zeros(length(scenariooptions))
    for obj in scenmod.objects
        if obj isa Balance
            if TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(obj))) == "Hydro"
                enekvglobal = 1.0 # if no energy equivalent, assume inflow is already demoninated in GWh
                if haskey(obj.metadata, TuLiPa.GLOBALENEQKEY)
                    enekvglobal = obj.metadata[TuLiPa.GLOBALENEQKEY]
                end
                for rhsterm in TuLiPa.getrhsterms(obj)
                    for (i, scenario) in enumerate(scenariooptions)
                        time = get_scentnormal(simtime, scenario, input)
                        totalsumenergyinflow[i] += TuLiPa.getparamvalue(rhsterm, time, scenmod.scendelta)*enekvglobal*factoroptions[i]
                    end
                end
            end
        end
    end

    # Fit values to normal distribution
    n = fit(Normal, totalsumenergyinflow, weightsoptions)
    quantiles = [i for i in 1-scenmod.maxquantile:(2*scenmod.maxquantile-1)/(numscen-1):scenmod.maxquantile] # get quantiles from maxquantile
    qvalues = quantile.(n, quantiles) # get quantile values from distribution
    
    # Could also use probability density of the quantiles in the calculation of the weights
    if scenmod.usedensity
        d = pdf.(n, qvalues) # get probability density for each quantile
    else
        d = [1.0 for i in numscen] # ignore probabilty density
    end

    # How much should the inflow in the scenario be adjusted so that it is similar to the quantile?
    for i in 1:numscen
        qvalue = qvalues[i]
        idx = findmin(abs.(totalsumenergyinflow.-qvalue))[2]
        scenmod.scenarios[i] = scenariooptions[idx]
        scenmod.inflowfactors[i] = qvalue/totalsumenergyinflow[idx]
    end

    # How much should each scenario be weighted - combination of weighting function and probability density
    x = collect(-numscen+1:2:numscen-1)
    y = (scenmod.a .* x .^ 2 .+ x .* scenmod.b .+ scenmod.c) .* d
    for i in numscen
        scenmod.scenarios[i].p_weather = y[i]/sum(y)
    end
    @assert sum([scenario.p_weather for scenario in scenmod.scenarios]) == 1

    return
end

function choose_scenarios!(scenmod::InflowClusteringMethod{WeatherScenario}, scenmodmethodoptions::AbstractScenarioModellingMethod, simtime::TuLiPa.ProbTime, input::AbstractJulESInput)
    numscen = length(scenmod.scenarios)
    scenariooptions = get_scenarios(scenmodmethodoptions)
    weightsoptions = [get_probability(scenario) for scenario in scenariooptions]
    inflowfactoroptions = get_inflowfactors(scenmodmethodoptions)

    # Calculate total energy inflow in the system for each part of the scenariodelta
    sumenergyinflow = zeros(length(scenariooptions))
    partsumenergyinflow = zeros(scenmod.parts ,length(scenariooptions))
    scendeltapart = scenmod.scendelta/scenmod.parts

	parts = scenmod.parts

    for obj in scenmod.objects
        if obj isa TuLiPa.Balance
            if TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(obj))) == "Hydro"
                enekvglobal = 1.0 # if no energy equivalent, assume inflow is already demoninated in GWh
                if haskey(obj.metadata, TuLiPa.GLOBALENEQKEY)
                    enekvglobal = obj.metadata[TuLiPa.GLOBALENEQKEY]
                end
                for rhsterm in TuLiPa.getrhsterms(obj)
                    for (i, scenario) in enumerate(scenariooptions)
                        time = get_scentnormal(simtime, scenario, input)
                        sumenergyinflow[i] += TuLiPa.getparamvalue(rhsterm, time, scenmod.scendelta)*enekvglobal*inflowfactoroptions[i]
                        for j in 1:parts
                            partsumenergyinflow[j,i] += TuLiPa.getparamvalue(rhsterm, time + scendeltapart*(j-1), scendeltapart)*enekvglobal*inflowfactoroptions[i]
                        end
                    end
                end
            end
        end
    end

    # Cluster inflow scenarios together
    r = kmeans(partsumenergyinflow, numscen)
    assignments = r.assignments

    # # Take a look at the clustering
    # plots = plot()
    # for i in 1:length(scentimesoptions)
    #     plot!(plots, partsumenergyinflow[:,i], color = palette(:default)[assignments[i]], labels=string(assignments[i]))
    # end
    # display(plots)

    # Find scenario in middle of each cluster
    for i in 1:numscen
        idxs = findall(x -> x == i, assignments)

        # Scenario in middle of cluster
        clustersumenergyinflows = sumenergyinflow[idxs]
        meanclustersumenergyinflow = mean(clustersumenergyinflows)
        clusteridx = findmin(x->abs(x-meanclustersumenergyinflow), clustersumenergyinflows)[2]
        totalidx = findfirst(x -> x == clustersumenergyinflows[clusteridx], sumenergyinflow)
        scenmod.scenarios[i] = deepcopy(scenariooptions[totalidx])

        # Adjust scenario to represent actual middle of cluster
        scenmod.inflowfactors[i] = meanclustersumenergyinflow/sumenergyinflow[totalidx]

        # Weight based on amount of scenarios in cluster and weight of options
        scenmod.scenarios[i].p_weather = sum(weightsoptions[idxs])
    end
    if !isapprox(sum([scenario.p_weather for scenario in scenmod.scenarios]), 1, atol=0.0001)
        println(sum([scenario.p_weather for scenario in scenmod.scenarios]))
        error("Sum of scenarios not 1")
    end


    return
end

# Scale inflow of modelobjects given a scenario modelling method
# Only implemented in subsystem models at the moment. If used in prognosis we loose correlation to rest of market
# In practice this means Prognosis ignores part of scenario generation (SumInflow of cluster), especially if SumInflowQuantileMethod
perform_scenmod!(scenmod::AbstractScenarioModellingMethod, scenix, objects) = nothing # generic fallback
function perform_scenmod!(scenmod::Union{InflowClusteringMethod{WeatherScenario},SumInflowQuantileMethod{WeatherScenario}}, scenix, objects) # inflow methods has field factor
    inflowfactors = get_inflowfactors(scenmod)
    inflowfactors[scenix] == 1.0 && return
    for obj in objects
        if obj isa TuLiPa.BaseBalance
            if TuLiPa.getinstancename(TuLiPa.getid(TuLiPa.getcommodity(obj))) == "Hydro"
                for rhsterm in TuLiPa.getrhsterms(obj)
                    if rhsterm.param isa TuLiPa.TwoProductParam
                        rhsterm.param = TuLiPa.TwoProductParam(rhsterm.param.param1, TuLiPa.ConstantParam(inflowfactors[scenix]))
                    else
                        rhsterm.param = TuLiPa.TwoProductParam(rhsterm.param, TuLiPa.ConstantParam(inflowfactors[scenix]))
                    end
                end
            end
        end
    end
    return
end

# Renumber scenarios
function renumber_scenmodmethod!(scenmod::AbstractScenarioModellingMethod)
    for (i, scenario) in enumerate(get_scenarios(scenmod))
        scenario.parentscenario = i
    end
end