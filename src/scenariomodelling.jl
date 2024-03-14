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
    scendelta::Millisecond

    function InflowClusteringMethod(scenarios, objects, parts, scendelta)
        inflowfactors = Vector{Float64}(undef, length(scenarios))
        return new{eltype(scenarios)}(scenarios, inflowfactors, objects, parts, scendelta)
    end
end
mutable struct SumInflowQuantileMethod{T <: AbstractScenario} <: AbstractScenarioModellingMethod
    scenarios::Vector{T}
    inflowfactors::Vector{Float64}
    objects::Vector
    maxquantile::Float64 # parameter
    a::Float64 # parameter
    b::Float64 # parameter
    c::Float64 # parameter
    usedensity::Bool # parameter
    scendelta::Millisecond # parameter
    function SumInflowQuantileMethod(scenarios, objects, maxquantile, a, b, c, scendelta; usedensity=false)
        inflowfactors = Vector{Float64}(undef, length(scenarios))
        return new{eltype(scenarios)}(scenarios, inflowfactors, objects, maxquantile, a, b, c, scendelta, usedensity)
    end
end

getchanges(scenmodmethod::NoScenarioModellingMethod) = (scenmodmethod.scenarios)
getchanges(scenmodmethod::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}) = (scenmodmethod.scenarios, scenmodmethod.inflowfactors)

function setchanges(scenmodmethod::NoScenarioModellingMethod, changes::Tuple{Vector{WeatherScenario}})
    scenmodmethod.scenarios = first(changes)
end
function setchanges(scenmodmethod::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}, changes::Tuple{Vector{WeatherScenario}, Vector{Float64}})
    scenarios, inflowfactors = changes

    scenmodmethod.scenarios = scenarios
    scenmodmethod.inflowfactors = inflowfactors
end

getinflowfactors(scenmodmethod::AbstractScenarioModellingMethod) = [1/length(scenmodmethod.scenarios) for s in 1:length(scenmodmethod.scenarios)]
getinflowfactors(scenmodmethod::Union{SumInflowQuantileMethod{WeatherScenario},InflowClusteringMethod{WeatherScenario}}) = scenmodmethod.inflowfactors

function choose_scenarios!(scenmodmethod::SumInflowQuantileMethod{WeatherScenario}, scenmodmethodoptions::AbstractScenarioModellingMethod, simtime::ProbTime, input::AbstractJulESInput)
    numscen = length(scenmodmethod.scenarios)
    scenariooptions = scenmodmethodoptions.scenarios
    weightsoptions = [getprobability(scenario) for scenario in scenariooptions]
    factoroptions = getinflowfactors(scenmodmethodoptions)
    
    # Calculate total energy inflow in the system for the scenariodelta
    totalsumenergyinflow = zeros(length(scenariooptions))
    for obj in scenmodmethod.objects
        if obj isa Balance
            if getinstancename(getid(getcommodity(obj))) == "Hydro"
                enekvglobal = 1.0 # if no energy equivalent, assume inflow is already demoninated in GWh
                if haskey(obj.metadata, GLOBALENEQKEY)
                    enekvglobal = obj.metadata[GLOBALENEQKEY]
                end
                for rhsterm in getrhsterms(obj)
                    for (i, scenario) in enumerate(scenariooptions)
                        time = getscenariotnormal(simtime, scenario, input)
                        totalsumenergyinflow[i] += getparamvalue(rhsterm, time, scenmodmethod.scendelta)*enekvglobal*factoroptions[i]
                    end
                end
            end
        end
    end

    # Fit values to normal distribution
    n = fit(Normal, totalsumenergyinflow, weightsoptions)
    quantiles = [i for i in 1-scenmodmethod.maxquantile:(2*scenmodmethod.maxquantile-1)/(numscen-1):scenmodmethod.maxquantile] # get quantiles from maxquantile
    qvalues = quantile.(n, quantiles) # get quantile values from distribution
    
    # Could also use probability density of the quantiles in the calculation of the weights
    if scenmodmethod.usedensity
        d = pdf.(n, qvalues) # get probability density for each quantile
    else
        d = [1.0 for i in numscen] # ignore probabilty density
    end

    # How much should the inflow in the scenario be adjusted so that it is similar to the quantile?
    for i in 1:numscen
        qvalue = qvalues[i]
        idx = findmin(abs.(totalsumenergyinflow.-qvalue))[2]
        scenmodmethod.scenarios[i] = scenariooptions[idx]
        scenmodmethod.inflowfactors[i] = qvalue/totalsumenergyinflow[idx]
    end

    # How much should each scenario be weighted - combination of weighting function and probability density
    x = collect(-numscen+1:2:numscen-1)
    y = (scenmodmethod.a .* x .^ 2 .+ x .* scenmodmethod.b .+ scenmodmethod.c) .* d
    for i in numscen
        scenmodmethod.scenarios[i].p_weather = y[i]/sum(y)
    end
    @assert sum([scenario.p_weather for scenario in scenmodmethod.scenarios]) == 1

    return
end

function choose_scenarios!(scenmodmethod::InflowClusteringMethod{WeatherScenario}, scenmodmethodoptions::AbstractScenarioModellingMethod, simtime::ProbTime, input::AbstractJulESInput)
    numscen = length(scenmodmethod.scenarios)
    scenariooptions = scenmodmethodoptions.scenarios
    weightsoptions = [getprobability(scenario) for scenario in scenariooptions]
    inflowfactoroptions = getinflowfactors(scenmodmethodoptions)

    # Calculate total energy inflow in the system for each part of the scenariodelta
    sumenergyinflow = zeros(length(scenariooptions))
    partsumenergyinflow = zeros(scenmodmethod.parts ,length(scenariooptions))
    scendeltapart = scenmodmethod.scendelta/scenmodmethod.parts

	parts = scenmodmethod.parts

    for obj in scenmodmethod.objects
        if obj isa Balance
            if getinstancename(getid(getcommodity(obj))) == "Hydro"
                enekvglobal = 1.0 # if no energy equivalent, assume inflow is already demoninated in GWh
                if haskey(obj.metadata, GLOBALENEQKEY)
                    enekvglobal = obj.metadata[GLOBALENEQKEY]
                end
                for rhsterm in getrhsterms(obj)
                    for (i, scenario) in enumerate(scenariooptions)
                        time = getscenariotnormal(simtime, scenario, input)
                        sumenergyinflow[i] += getparamvalue(rhsterm, time, scenmodmethod.scendelta)*enekvglobal*inflowfactoroptions[i]
                        for j in 1:parts
                            partsumenergyinflow[j,i] += getparamvalue(rhsterm, time + scendeltapart*(j-1), scendeltapart)*enekvglobal*factoroptions[i]
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
        scenmodmethod.scenarios[i] = scenariooptions[totalidx]

        # Adjust scenario to represent actual middle of cluster
        scenmodmethod.inflowfactors[i] = meanclustersumenergyinflow/sumenergyinflow[totalidx]

        # Weight based on amount of scenarios in cluster and weight of options
        scenmodmethod.scenarios[i].p_weather = sum(weightsoptions[idxs])
    end
    @assert sum([scenario.p_weather for scenario in scenmodmethod.scenarios]) == 1

    return
end

# Scale inflow of modelobjects given a scenario modelling method
# Only implemented in subsystem models at the moment. If used in prognosis we loose correlation to rest of market
# In practice this means Prognosis ignores part of scenario generation (SumInflow of cluster), especially if SumInflowQuantileMethod
perform_scenmod!(scenmodmethod::AbstractScenarioModellingMethod, scenarioix, objects) = nothing # generic fallback
function perform_scenmod!(scenmodmethod::Union{InflowClusteringMethod{WeatherScenario},SumInflowQuantileMethod{WeatherScenario}}, scenarioix, objects) # inflow methods has field factor
    inflowfactors = getinflowfactors
    inflowfactors[scenarioix] == 1.0 && return
    for obj in objects
        if obj isa BaseBalance
            if getinstancename(getid(getcommodity(obj))) == "Hydro"
                for rhsterm in getrhsterms(obj)
                    if rhsterm.param isa TwoProductParam
                        rhsterm.param = TwoProductParam(rhsterm.param.param1, ConstantParam(inflowfactors[scenarioix]))
                    else
                        rhsterm.param = TwoProductParam(rhsterm.param, ConstantParam(inflowfactors[scenarioix]))
                    end
                end
            end
        end
    end
    return
end

# Renumber scenarios
function renumber_scenmodmethod!(scenmodmethod::AbstractScenarioModellingMethod)
    for (i, scenario) in enumerate(scenmodmethod.scenarios)
        scenario.parentscenario = i
    end
end