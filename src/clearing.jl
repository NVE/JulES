# Initialize market clearing problem and solve for first time step
function clearing_init(elements, t, clearingdays, cutslocal, nonstoragestateslocal)
    elements1 = copy(elements)

    # Add horizon to dataelements
    battery_horizon = SequentialHorizon(clearingdays*12, Hour(2)) # TODO: Replace with user settings
    hydro_horizon = SequentialHorizon(clearingdays, Hour(24)) # TODO: Higher time resolution for small PHS
    power_horizon = SequentialHorizon(clearingdays*12, Hour(2))

    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Power", 
            (HORIZON_CONCEPT, power_horizon)))
    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Hydro", 
            (HORIZON_CONCEPT, hydro_horizon)))
    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Battery", 
            (HORIZON_CONCEPT, battery_horizon)))

    # Make modelobjects and add upper slack variable for power production
    modelobjects = getmodelobjects(elements1)
    addPowerUpperSlack!(modelobjects)

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

    # Set end values from cuts generated in stochastic subsystem problems
    clearingendvaluesdict = Dict()
    for cuts in cutslocal
        for (state, value) in cuts.slopes[cuts.cutix] # slopes of last cut
            clearingendvaluesdict[first(getvarout(state))] = value
        end
    end
    endvalues = [-clearingendvaluesdict[getid(obj)] for obj in clearingstorages]
    clearingendvaluesid = Id(BOUNDARYCONDITION_CONCEPT,"EndValue")
    clearingendvaluesobj = EndValues(clearingendvaluesid, clearingstorages)
    push!(clearing.objects, clearingendvaluesobj)
    updateendvalues!(clearing, clearingendvaluesobj, endvalues)

    # Set outgoing states for non-storage variables
    nonstoragestatesmean = Dict{StateVariableInfo, Float64}()
    for nonstoragestate in keys(nonstoragestateslocal[1])
        # (id, ix) = getvarout(nonstoragestate)
        # newnonstoragestate = StateVariableInfo(getvarin(nonstoragestate), (id, ix*2)) # Quick fix different resolution short prognosis (2-hourly) and clearing (hourly)
        # nonstoragestatesmean[newnonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate]/2 for i in eachindex(nonstoragestateslocal)])
        nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    end
    setoutgoingstates!(clearing, nonstoragestatesmean)

    # Update and solve problem
    update!(clearing, t)
    solve!(clearing)
    
    return clearing, clearingstorages, nonstoragestatesmean, clearingendvaluesdict
end

# Run market clearing for new time step
function clearing!(clearing, clearingstorages, startstates, cutslocal, clearingendvaluesdict, nonstoragestateslocal, nonstoragestatesmean, detailedrescopl, enekvglobaldict)
        
    # Update startstates for all state variables, equals end state of last market clearing
    setstartstates!(clearing, clearingstorages, startstates)
        
    # Update storage end values from cuts
    for cuts in cutslocal
        for (state, value) in cuts.slopes[cuts.cutix] # slopes of last cut
            clearingendvaluesdict[first(getvarout(state))] = value
        end
    end
    clearingendvaluesobj = clearing.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), clearing.objects)]
    endvalues = [-clearingendvaluesdict[getid(obj)] for obj in clearingstorages]
    updateendvalues!(clearing, clearingendvaluesobj, endvalues)

    # Update end states for non storage variables
    for nonstoragestate in keys(nonstoragestateslocal[1])
        nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    end
    setoutgoingstates!(clearing, nonstoragestatesmean)

    # Update and solve
    update!(clearing, tnormal)
    solve!(clearing)

    # Get start states for next iteration
    startstates = getstartstates(clearing, detailedrescopl, enekvglobaldict)
    
    return startstates
end

# Collect startstates for next iteration. Also calculate for aggregated reservoirs
function startstates_init(clearing, detailedrescopl, prob, t)

    enekvglobaldict = JSON.parsefile("data_fra_dynmodell/magasin_enekvglobal.json") # TODO: only 2021 values, add global enekv as metadata for storages

    startstates_ = getstatevariables(clearing.objects)
    getoutgoingstates!(clearing, startstates_)
    startstates = Dict()
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
    
    return startstates, enekvglobaldict
end

function getstartstates(clearing, detailedrescopl, enekvglobaldict)
    startstates_ = getstatevariables(clearing.objects)
    getoutgoingstates!(clearing, startstates_)
    
    for var in keys(startstates_)
        # value = round(startstates_[var], digits=10) # avoid approx 0 negative values, ignored by solvers so no problem?
        startstates[getinstancename(first(getvarout(var)))] = startstates_[var]
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
    
    return startstates
end