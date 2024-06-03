# Initialize market clearing problem and solve for first time step
function clearing_init(probmethod::ProbMethod, elements::Vector{DataElement}, t::ProbTime, clearingduration::Millisecond, cpdp::Millisecond, cpdh::Millisecond, startstates::Dict{String, Float64}, masterslocal::Vector{Prob}, cutslocal::Vector{SimpleSingleCuts}, nonstoragestateslocal::Vector{Dict}, settings::Dict)
    elements1 = copy(elements)

    # Add horizon to dataelements
    battery_horizon = SequentialHorizon(ceil(Int64, clearingduration/cpdp), cpdp)
    hydro_horizon = SequentialHorizon(ceil(Int64, clearingduration/cpdh), cpdh) # TODO: Higher time resolution for small PHS
    power_horizon = SequentialHorizon(ceil(Int64, clearingduration/cpdp), cpdp)

    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Power", 
            (HORIZON_CONCEPT, power_horizon)))
    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Hydro", 
            (HORIZON_CONCEPT, hydro_horizon)))
    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Battery", 
            (HORIZON_CONCEPT, battery_horizon)))

    # Make modelobjects and add upper slack variable for power production
    modelobjects = getmodelobjects(elements1)
    addPowerUpperSlack!(modelobjects)

    # Initialize cuts (has to be added to modelobjects so that all the variables are built at the same time)
    varendperiod = Dict()
    for cuts in cutslocal
        for statevar in cuts.statevars
            (varid, varix) = getvarout(statevar)
            varendperiod[varid] = getnumperiods(gethorizon(modelobjects[varid]))
        end
        cutid = getid(cuts)
        modelobjects[cutid] = cuts
    end

    # Build problem from modelobjects and optimizer
    clearing = buildprob(probmethod, modelobjects)

    # Set start storages
    setstartstates!(clearing, getstorages(getobjects(clearing)), startstates)

    # Update cuts
    for cuts in cutslocal
        updatecuts!(clearing, cuts, varendperiod)
    end

    # Set outgoing states for non-storage variables
    nonstoragestatesmean = Dict{StateVariableInfo, Float64}()
    # for nonstoragestate in keys(nonstoragestateslocal[1])
    #     # (id, ix) = getvarout(nonstoragestate)
    #     # newnonstoragestate = StateVariableInfo(getvarin(nonstoragestate), (id, ix*2)) # Quick fix different resolution short prognosis (2-hourly) and clearing (hourly)
    #     # nonstoragestatesmean[newnonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate]/2 for i in eachindex(nonstoragestateslocal)])
    #     nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    # end
    # setoutgoingstates!(clearing, nonstoragestatesmean)
    setoutgoingstates!(clearing, nonstoragestateslocal[1])

    # State dependent hydropower production and pumping.
    getstatedependentprod(settings["problems"]["clearing"]) && statedependentprod!(clearing, startstates, init=true)
    getstatedependentpump(settings["problems"]["clearing"]) && statedependentpump!(clearing, startstates)

    # Update and solve
    update!(clearing, t)
    getheadlosscost(settings["problems"]["clearing"]) && updateheadlosscosts!(ReservoirCurveSlopeMethod(), clearing, masterslocal, t)
    setminstoragevalue!(clearing, minstoragevaluerule) # add spill cost to water value so that the reservoir is not emptied when watervalue = 0

    solve!(clearing)
    
    return clearing, nonstoragestatesmean, varendperiod
end

# Run market clearing for new time step
function clearing!(clearing::Prob, t::ProbTime, startstates::Dict{String, Float64}, masterslocal::Vector{Prob}, cutslocal::Vector{SimpleSingleCuts}, nonstoragestateslocal::Vector{Dict}, nonstoragestatesmean, detailedrescopl::Dict, enekvglobaldict::Dict, varendperiod::Dict, clearingtimes::Matrix{Float64}, step::Int, settings::Dict)
    clearingtimes[step-1, 3] = @elapsed begin
        # Update startstates for all state variables, equals end state of last market clearing
        setstartstates!(clearing, getobjects(clearing), startstates) # TODO: Also store actual statevariables to update more efficiently?

        # Update cuts in problem
        for cuts in cutslocal
            updatecuts!(clearing, cuts, varendperiod)
        end

        # Update end states for non storage variables
        # for nonstoragestate in keys(nonstoragestateslocal[1])
        #     # (id, ix) = getvarout(nonstoragestate)
        #     # newnonstoragestate = StateVariableInfo(getvarin(nonstoragestate), (id, ix*2)) # Quick fix different resolution short prognosis (2-hourly) and clearing (hourly)
        #     # nonstoragestatesmean[newnonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate]/2 for i in eachindex(nonstoragestateslocal)])
        #     nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
        # end
        # setoutgoingstates!(clearing, nonstoragestatesmean)
        setoutgoingstates!(clearing, nonstoragestateslocal[1])

        # Statedependent hydropower production
        getstatedependentprod(settings["problems"]["clearing"]) && statedependentprod!(clearing, startstates)
        getstatedependentpump(settings["problems"]["clearing"]) && statedependentpump!(clearing, startstates)

        # Update and solve
        clearingtimes[step-1, 1] = @elapsed update!(clearing, t)
        getheadlosscost(settings["problems"]["clearing"]) && updateheadlosscosts!(ReservoirCurveSlopeMethod(), clearing, masterslocal, t)
        setminstoragevalue!(clearing, minstoragevaluerule) # add spill cost to water value so that the reservoir is not emptied when watervalue = 0

        clearingtimes[step-1, 2] = @elapsed solve!(clearing)

        # Get start states for next iteration
        getstartstates!(clearing, detailedrescopl, enekvglobaldict, startstates)
    end

end

function getstartstates!(clearing::Prob, detailedrescopl::Dict, enekvglobaldict::Dict, startstates::Dict{String, Float64})
    startstates_ = getstates(getobjects(clearing))
    getoutgoingstates!(clearing, startstates_)
    
    for var in keys(startstates_)
        value = round(startstates_[var], digits=10) # avoid approx 0 negative values, ignored by solvers so no problem?
        startstates[getinstancename(first(getvarout(var)))] = value
    end

    for area in Set(values(detailedrescopl))
        startstates["Reservoir_" * area * "_hydro_reservoir"] = 0.0
    end
    for res in keys(detailedrescopl)
        resname = "Reservoir_" * res
        areaname = "Reservoir_" * detailedrescopl[res] * "_hydro_reservoir"
        startstates[areaname] += startstates[resname] * enekvglobaldict[res]
    end

    # Avoid reservoirs being filled more than max, gives infeasible solution
    # - If aggregated reservoir capacity is lower than the sum capacities
    # - If reservoir is full in model, numerical tolerance can bring variable value slightly over cap
    # - TODO: Add warning/logging if this happens
    for resname in keys(startstates)
        resmax = resname * "_max"
        if haskey(startstates, resmax)
            if startstates[resname] > startstates[resmax]
                startstates[resname] = startstates[resmax]
            end
        end
    end
end

function minstoragevaluerule(storage::Storage)
    minstoragevalues = Dict{String, Float64}()
    minstoragevalues["Battery"] = 0.0
    minstoragevalues["Hydro"] = 0.001
    commodity = getinstancename(getid(getcommodity(getbalance(storage))))
    return get(minstoragevalues, commodity, 0.0)
end

function setminstoragevalue!(problem::Prob, costrule::Function)
    for modelobject in getobjects(problem)
        if modelobject isa Storage
            id = getid(modelobject)
            balance = getbalance(modelobject)
            horizon = gethorizon(balance)
            T = getnumperiods(horizon)
            coeff = getobjcoeff(problem, id, T)
            cost = costrule(modelobject)
            newcoeff = min(-cost, coeff)
            if !(coeff ≈ newcoeff)
                setobjcoeff!(problem, id, T, newcoeff)
            end
        end
    end
    return
end