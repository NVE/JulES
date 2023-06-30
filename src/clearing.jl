# Initialize market clearing problem and solve for first time step
function clearing_init(elements, t, clearingduration, cpdp, cpdh, masterslocal, cutslocal, nonstoragestateslocal)
    elements1 = copy(elements)

    # Add horizon to dataelements
    battery_horizon = SequentialHorizon(ceil(Int64, clearingduration/cpdp), cpdp) # TODO: Replace with user settings
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
    remove_hydrorampingwithout!(modelobjects)

    # Initialize cuts
    varendperiod = Dict()
    for cuts in cutslocal
        # Replace local stochastic object with version from clearing, also store numperiods of clearingobject
        for (i,obj) in enumerate(cuts.objects)
            objid = getid(obj)
            clearingobj = modelobjects[objid]
            varendperiod[objid] = getnumperiods(gethorizon(clearingobj))
        end

        cutid = getid(cuts)
        modelobjects[cutid] = cuts
    end

    # Build problem from modelobjects and optimizer
    # model = Model(HiGHS.Optimizer)
    # clearing = JuMP_Prob(modelobjects, model)
    clearing = HiGHS_Prob(modelobjects)

    # Set start storages
    shorttermstorages = getshorttermstorages(getobjects(clearing), Hour(10))
    clearingstorages = getstorages(getobjects(clearing))
    longtermstorages = setdiff(clearingstorages, shorttermstorages)
    setstartstoragepercentage!(clearing, shorttermstorages, t, 50)
    setstartstoragepercentage!(clearing, longtermstorages, t, 65)

    # Update cuts
    for cuts in cutslocal
        updatecuts!(clearing, cuts, varendperiod)
    end

    # Set outgoing states for non-storage variables
    nonstoragestatesmean = Dict{StateVariableInfo, Float64}()
    for nonstoragestate in keys(nonstoragestateslocal[1])
        # (id, ix) = getvarout(nonstoragestate)
        # newnonstoragestate = StateVariableInfo(getvarin(nonstoragestate), (id, ix*2)) # Quick fix different resolution short prognosis (2-hourly) and clearing (hourly)
        # nonstoragestatesmean[newnonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate]/2 for i in eachindex(nonstoragestateslocal)])
        nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    end
    setoutgoingstates!(clearing, nonstoragestatesmean)

    # State dependent hydropower production and pumping.
    statedependentprod_init!(clearing, 65, t)
    statedependentpump_init!(clearing, 65, t)

    # Update and solve
    update!(clearing, t)
    updateheadlosscosts!(ReservoirCurveSlopeMethod(), clearing, masterslocal, t)
    setminstoragevalue!(clearing, minstoragevaluerule) # add spill cost to water value so that the reservoir is not emptied when watervalue = 0

    solve!(clearing)
    
    return clearing, nonstoragestatesmean, varendperiod
end

# Run market clearing for new time step
function clearing!(clearing, t, startstates, masterslocal, cutslocal, nonstoragestateslocal, nonstoragestatesmean, detailedrescopl, enekvglobaldict, varendperiod)
        
    # Update startstates for all state variables, equals end state of last market clearing
    setstartstates!(clearing, clearing.objects, startstates) # TODO: Also store actual statevariables to update more efficiently?

    # Update cuts in problem
    for cuts in cutslocal
        updatecuts!(clearing, cuts, varendperiod)
    end

    # Update end states for non storage variables
    for nonstoragestate in keys(nonstoragestateslocal[1])
        # (id, ix) = getvarout(nonstoragestate)
        # newnonstoragestate = StateVariableInfo(getvarin(nonstoragestate), (id, ix*2)) # Quick fix different resolution short prognosis (2-hourly) and clearing (hourly)
        # nonstoragestatesmean[newnonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate]/2 for i in eachindex(nonstoragestateslocal)])
        nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    end
    setoutgoingstates!(clearing, nonstoragestatesmean)

    # Statedependent hydropower production
    statedependentprod!(clearing, startstates)
    statedependentpump!(clearing, startstates)

    # Update and solve
    update!(clearing, t)
    updateheadlosscosts!(ReservoirCurveSlopeMethod(), clearing, masterslocal, t)
    setminstoragevalue!(clearing, minstoragevaluerule) # add spill cost to water value so that the reservoir is not emptied when watervalue = 0

    solve!(clearing)

    # Get start states for next iteration
    getstartstates!(clearing, detailedrescopl, enekvglobaldict, startstates)
end

# Collect startstates for next iteration. Also calculate for aggregated reservoirs
function startstates_init(clearing, detailedrescopl, enekvglobaldict, prob, t)

    startstates_ = getstatevariables(clearing.objects)
    getoutgoingstates!(clearing, startstates_)
    startstates = Dict{String, Float64}()

    for var in keys(startstates_)
        startstates[getinstancename(first(getvarout(var)))] = startstates_[var]
    end

    for area in Set(values(detailedrescopl))
        resname = "Reservoir_" * area * "_hydro_reservoir"
        startstates[resname] = 0.0

        for obj in prob.objects
            if getinstancename(getid(obj)) == resname
                startstates[resname * "_max"] = getparamvalue(obj.ub, t, MsTimeDelta(Millisecond(0)))
            end
        end
    end

    for res in keys(detailedrescopl)
        resname = "Reservoir_" * res
        areaname = "Reservoir_" * detailedrescopl[res] * "_hydro_reservoir"
        startstates[areaname] += startstates[resname] * enekvglobaldict[res]
    end

    for area in Set(values(detailedrescopl)) # aggregated reservoirs cannot be filled more than max
        resname = "Reservoir_" * area * "_hydro_reservoir"
        if startstates[resname] > startstates[resname * "_max"]
            startstates[resname] = startstates[resname * "_max"]
        end
    end
    
    return startstates
end

function getstartstates!(clearing, detailedrescopl, enekvglobaldict, startstates)
    startstates_ = getstatevariables(clearing.objects)
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

    for area in Set(values(detailedrescopl)) # aggregated reservoirs cannot be filled more than max
        resname = "Reservoir_" * area * "_hydro_reservoir"
        if startstates[resname] > startstates[resname * "_max"]
            startstates[resname] = startstates[resname * "_max"]
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