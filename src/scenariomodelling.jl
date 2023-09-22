abstract type ScenarioModellingMethod end
mutable struct NoScenarioModellingMethod <: ScenarioModellingMethod
    scentimes::Vector{Tuple{Any, Any, Int64}}
    weights::Vector{Float64}
    function NoScenarioModellingMethod(numscen, totalscentimes)
        @assert numscen == length(totalscentimes)
        new(totalscentimes, [1/length(totalscentimes) for i in 1:length(totalscentimes)])
    end
end
# mutable struct ResidualLoadMethod <: ScenarioModellingMethod # choose scenario based on residual load (also energy inflow)
# end
mutable struct InflowClusteringMethod <: ScenarioModellingMethod
    scentimes::Vector{Tuple{Any, Any, Int64}}
    weights::Vector{Float64}
    factors::Vector{Float64}
    parts::Int
    function InflowClusteringMethod(numscen, parts)
        return new(Vector{Tuple{Any, Any, Int64}}(undef, numscen), Vector{Float64}(undef, numscen), Vector{Float64}(undef, numscen), parts)
    end
end
mutable struct SumInflowQuantileMethod <: ScenarioModellingMethod
    scentimes::Vector{Tuple{Any, Any, Int64}}
    weights::Vector{Float64}
    factors::Vector{Float64}
    maxquantile::Float64 # parameter
    a::Float64 # parameter
    b::Float64 # parameter
    c::Float64 # parameter
    usedensity::Bool # parameter
    function SumInflowQuantileMethod(numscen, maxquantile, a, b, c; usedensity=false)
        return new(Vector{Tuple{Any, Any, Int64}}(undef, numscen), Vector{Float64}(undef, numscen), Vector{Float64}(undef, numscen), maxquantile, a, b, c, usedensity)
    end
end

function scenariomodelling!(scenmodmethod::ScenarioModellingMethod, objects, numscen, totalscentimes, scendelta)
    scenmodmethod.scentimes = totalscentimes
end

function scenariomodelling!(scenmodmethod::SumInflowQuantileMethod, objects, numscen, totalscentimes, scendelta)
    # Calculate total energy inflow in the system for the scenariodelta
    totalsumenergyinflow = zeros(length(totalscentimes))
    for obj in objects
        if obj isa Balance
            if getinstancename(getid(getcommodity(obj))) == "Hydro"
                enekvglobal = obj.metadata[GLOBALENEQKEY]
                for rhsterm in getrhsterms(obj)
                    for (i, (scentnormal,scentphasein)) in enumerate(totalscentimes)
                        totalsumenergyinflow[i] += getparamvalue(rhsterm, scentnormal, scendelta)*enekvglobal
                    end
                end
            end
        end
    end

    # Fit values to normal distribution
    n = fit(Normal, totalsumenergyinflow)
    quantiles = [i for i in 1-scenmodmethod.maxquantile:(2*scenmodmethod.maxquantile-1)/(numscen-1):scenmodmethod.maxquantile] # get quantiles from maxquantile
    qvalues = quantile.(n, quantiles) # get quantile values from distribution
    
    # Could also use probability density of the quantiles in the calculation of the weights
    if scenmodmethod.usedensity
        d = pdf.(n, qvalues) # get probability density for each quantile
    else
        d = [1.0 for i in numscen] # ignore probabilty density
    end

    # How much should each scenario be weighted - combination of weighting function and probability density
    x = collect(-numscen+1:2:numscen-1)
    y = (scenmodmethod.a .* x .^ 2 .+ x .* scenmodmethod.b .+ scenmodmethod.c) .* d
    scenmodmethod.weights = y/sum(y)

    # How much should the inflow in the scenario be adjusted so that it is similar to the quantile?
    for i in 1:numscen
        qvalue = qvalues[i]
        idx = findmin(abs.(totalsumenergyinflow.-qvalue))[2]
        scenmodmethod.scentimes[i] = totalscentimes[idx]
        scenmodmethod.factors[i] = qvalue/totalsumenergyinflow[idx]
    end
    return
end

function scenariomodelling!(scenmodmethod::InflowClusteringMethod, objects, numscen, totalscentimes, scendelta)
    # Calculate total energy inflow in the system for each part of the scenariodelta
    sumenergyinflow = zeros(length(totalscentimes))
    partsumenergyinflow = zeros(scenmodmethod.parts ,length(totalscentimes))
    scendeltapart = scendelta/scenmodmethod.parts

    for obj in objects
        if obj isa Balance
            if getinstancename(getid(getcommodity(obj))) == "Hydro"
                enekvglobal = 1.0 # if no energy equivalent, assume inflow is already demoninated in GWh
                if haskey(obj.metadata, GLOBALENEQKEY)
                    enekvglobal = obj.metadata[GLOBALENEQKEY]
                end
                for rhsterm in getrhsterms(obj)
                    for (i, (scentnormal,scentphasein)) in enumerate(totalscentimes)
                        sumenergyinflow[i] += getparamvalue(rhsterm, scentnormal, scendelta)*enekvglobal
                        for j in 1:parts
                            partsumenergyinflow[j,i] += getparamvalue(rhsterm, scentnormal + scendeltapart*(j-1), scendeltapart)*enekvglobal
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
    # for i in 1:length(totalscentimes)
    #     plot!(plots, partsumenergyinflow[:,i], color = palette(:default)[assignments[i]], labels=string(assignments[i]))
    # end
    # display(plots)

    # Find scenario in middle of each cluster
    for i in 1:numscen
        # Weight based on amount of scenarios in cluster
        idxs = findall(x -> x == i, assignments)
        scenmodmethod.weights[i] = length(idxs)/length(totalscentimes)

        # Scenario in middle of cluster
        clustersumenergyinflows = sumenergyinflow[idxs]
        meanclustersumenergyinflow = mean(clustersumenergyinflows)
        clusteridx = findmin(x->abs(x-meanclustersumenergyinflow), clustersumenergyinflows)[2]
        totalidx = findfirst(x -> x == clustersumenergyinflows[clusteridx], sumenergyinflow)
        scenmodmethod.scentimes[i] = totalscentimes[totalidx]

        # Adjust scenario to represent actual middle of cluster
        scenmodmethod.factors[i] = meanclustersumenergyinflow/sumenergyinflow[totalidx]
    end
    return
end

# Scale inflow of modelobjects given a scenario modelling method (only )
scaleinflow!(scenmodmethod::ScenarioModellingMethod, scenario, objects) = nothing # generic fallback
function scaleinflow!(scenmodmethod::Union{InflowClusteringMethod,SumInflowQuantileMethod}, scenario, objects) # inflow methods has field factor
    for obj in objects
        if obj isa BaseBalance
            if getinstancename(getid(getcommodity(obj))) == "Hydro"
                for rhsterm in getrhsterms(obj)
                    if rhsterm.param isa TwoProductParam
                        rhsterm.param = TwoProductParam(rhsterm.param.param1, ConstantParam(scenmodmethod.factors[scenario]))
                    else
                        rhsterm.param = TwoProductParam(rhsterm.param, ConstantParam(scenmodmethod.factors[scenario]))
                    end
                end
            end
        end
    end
end;