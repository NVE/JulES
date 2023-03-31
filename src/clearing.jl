function clearing_init(elements, t, clearingdays, cutslocal, nonstoragestateslocal)
    elements1 = copy(elements)

    battery_horizon = SequentialHorizon(clearingdays*12, Hour(2))
    hydro_horizon = SequentialHorizon(clearingdays, Hour(24)) # TODO: Higher time resolution for PHS
    power_horizon = SequentialHorizon(clearingdays*12, Hour(2))

    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Power", 
            (HORIZON_CONCEPT, power_horizon)))
    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Hydro", 
            (HORIZON_CONCEPT, hydro_horizon)))
    push!(elements1, getelement(COMMODITY_CONCEPT, "BaseCommodity", "Battery", 
            (HORIZON_CONCEPT, battery_horizon)))

    modelobjects = getmodelobjects(elements1)
    addPowerUpperSlack!(modelobjects)

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    clearing = JuMP_Prob(modelobjects, model)
    # clearing = HiGHS_Prob(modelobjects)

    shorttermstorages = getshorttermstorages(getobjects(clearing), Hour(10))
    clearingstorages = getstorages(getobjects(clearing))
    longtermstorages = setdiff(clearingstorages, shorttermstorages)
    setstartstoragepercentage!(clearing, shorttermstorages, t, 50)
    setstartstoragepercentage!(clearing, longtermstorages, t, 65)

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

    nonstoragestatesmean = Dict{StateVariableInfo, Float64}()
    for nonstoragestate in keys(nonstoragestateslocal[1])
        nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    end
    setoutgoingstates!(clearing, nonstoragestatesmean)

    update!(clearing, t)

    solve!(clearing)
    
    return clearing, clearingstorages, nonstoragestatesmean, clearingendvaluesdict
end

function clearing!(clearing, clearingstorages, startstates, cutslocal, clearingendvaluesdict, nonstoragestateslocal, nonstoragestatesmean, detailedrescopl, enekvglobaldict)
        
    setstartstates!(clearing, clearingstorages, startstates)
        
    for cuts in cutslocal
        for (state, value) in cuts.slopes[cuts.cutix] # slopes of last cut
            clearingendvaluesdict[first(getvarout(state))] = value
        end
    end
    clearingendvaluesobj = clearing.objects[findfirst(x -> getid(x) == Id(BOUNDARYCONDITION_CONCEPT,"EndValue"), clearing.objects)]
    endvalues = [-clearingendvaluesdict[getid(obj)] for obj in clearingstorages]
    updateendvalues!(clearing, clearingendvaluesobj, endvalues)

    for nonstoragestate in keys(nonstoragestateslocal[1])
        nonstoragestatesmean[nonstoragestate] = mean([nonstoragestateslocal[i][nonstoragestate] for i in eachindex(nonstoragestateslocal)])
    end
    setoutgoingstates!(clearing, nonstoragestatesmean)

    update!(clearing, tnormal)

    @time solve!(clearing)

    startstates = getstartstates(clearing, detailedrescopl, enekvglobaldict)
    
    return startstates
end

function startstates_init(clearing, detailedrescopl, prob, t)
    # Collect startreservoirs for next iteration. Also for aggregated reservoirs
    enekvglobaldict = JSON.parsefile("data_fra_dynmodell/magasin_enekvglobal.json") # TODO: 2021 values

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
        # value = round(startstates_[var], digits=10) # avoid approx 0 negative values
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