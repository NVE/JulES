function run_serial(config, datayear, scenarioyear, dataset)
    # Collect main information
    mainconfig = config["main"]
    settings = config[mainconfig["settings"]]
    numcores = mainconfig["numcores"]

    onlysubsystemmodel = false
    if !haskey(settings["problems"], "prognosis") && haskey(settings["problems"], "stochastic") && !haskey(settings["problems"], "clearing")
        onlysubsystemmodel = true
    end

    println("Time parameters")
    @time begin
        weekstart = mainconfig["weekstart"]
        
        scenarioyearstart = settings["time"]["scenarioyearstart"]
        scenarioyearstop = settings["time"]["scenarioyearstop"]
        datanumscen = scenarioyearstop - scenarioyearstart # scenarios to consider uncertainty for
        
        simulationyears = mainconfig["simulationyears"]
        extrasteps = mainconfig["extrasteps"]
        steplength = Millisecond(Hour(settings["time"]["steplength_hours"]))
        
        # Phasein settings
        phaseinoffset = steplength # phase in straight away from second stage scenarios
        phaseindelta = Millisecond(Day(settings["time"]["probtime"]["phaseindelta_days"])) # Phase in the second stage scenario over 5 weeks
        phaseinsteps = settings["time"]["probtime"]["phaseinsteps"] # Phase in second stage scenario in 5 steps

        # Make standard time and scenario uncertainty times
        tnormaltype = settings["time"]["probtime"]["normaltime"]
        tphaseintype = settings["time"]["probtime"]["phaseintime"]
        tnormal, datascenmodmethod = getprobtimes(datayear, weekstart, scenarioyear, datanumscen, tnormaltype, tphaseintype, phaseinoffset, phaseindelta, phaseinsteps)
        
        # How many time steps to run the simulation for
        steps = Int(ceil((getisoyearstart(datayear + simulationyears) - getisoyearstart(datayear)).value/steplength.value) + extrasteps);
    end

    println("Get data")
    @time begin
        elements = dataset["elements"]
        detailedrescopl = dataset["detailedrescopl"]
        addscenariotimeperiod_vector!(elements, scenarioyearstart, scenarioyearstop)

        if haskey(dataset, "progelements")
            progelements = dataset["progelements"]
            addscenariotimeperiod_vector!(progelements, scenarioyearstart, scenarioyearstop)
        else
            progelements = elements
        end
    end

    println("Make dummy objects") # for use in scenario modelling, validate elements and collect storages
    @time begin
        # Horizons are needed to build modelobjects, but not used in scenario modelling
        dummyperiods = 10
        dummyperiodduration = Millisecond(Hour(24))
        power_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)
        hydro_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)

        # Make dummy elements
        dummyobjects, dummydhh, dummydph = make_obj(elements, hydro_horizon, power_horizon, validate=true)

        # Make dummy prog elements
        if haskey(dataset, "progelements")
            dummyprogobjects, dummyphh, dummypph = make_obj(progelements, hydro_horizon, power_horizon, validate=true)
        else
            dummyprogobjects = dummyobjects
        end
    end

    println("Init scenario modelling for simulation, prognosis and stochastic") 
    @time begin
        # Simulation scenario modelling - choose scenarios for the whole simulation
        simnumscen = settings["scenariogeneration"]["simulation"]["numscen"]; @assert simnumscen <= datanumscen
        if simnumscen == datanumscen
            simscenmodmethod = datascenmodmethod
        else
            global simscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["simulation"], simnumscen)
            simscendelta = MsTimeDelta(Day(settings["scenariogeneration"]["simulation"]["scendelta_days"])) # scenario modelling based on the next 3 years, even though the scenario problems can be longer
            @time scenariomodelling!(simscenmodmethod, values(dummyprogobjects), simnumscen, datascenmodmethod, simscendelta) # see JulES/scenariomodelling.jl
            renumber_scenmodmethod!(simscenmodmethod)
        end

        # Prognosis scenario modelling - choose scenarios for the price prognosis models
        prognumscen = settings["scenariogeneration"]["prognosis"]["numscen"]; @assert prognumscen <= simnumscen
        if prognumscen == simnumscen
            progscenmodmethod = simscenmodmethod
        else
            global progscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["prognosis"], prognumscen)
            progscendelta = MsTimeDelta(Day(settings["scenariogeneration"]["prognosis"]["scendelta_days"])) # scenario modelling based on the next 3 years, even though the scenario problems can be longer
            @time scenariomodelling!(progscenmodmethod, values(dummyprogobjects), prognumscen, simscenmodmethod, progscendelta); # see JulES/scenariomodelling.jl
            prognumscen != simnumscen && renumber_scenmodmethod!(progscenmodmethod)
        end

        # Stochastic scenario modelling - choose scenarios for the price stochastic models
        stochnumscen = settings["scenariogeneration"]["stochastic"]["numscen"]; @assert stochnumscen <= prognumscen
        if stochnumscen == prognumscen
            stochscenmodmethod = progscenmodmethod
        else
            stochscendelta = MsTimeDelta(Day(settings["scenariogeneration"]["stochastic"]["scendelta_days"])) # scenario modelling based on the next 3 years, even though the scenario problems can be longer
            global stochscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["stochastic"], stochnumscen)
            @time scenariomodelling!(stochscenmodmethod, values(dummyobjects), stochnumscen, progscenmodmethod, stochscendelta); # see JulES/scenariomodelling.jl
        end
    end

    medendvaluesdicts = Dict[]
    startstates = Dict{String, Float64}()
    enekvglobaldict = Dict()
    medpriceslocal = nothing
    shortpriceslocal = nothing
    if haskey(settings["problems"], "prognosis")
        println("Init prognosis")
        @time begin
            # Set horizons for price prognosis models
            # All
            shorthorizonduration = Millisecond(Hour(settings["horizons"]["short"]["horizonduration_hours"]))

            # Long
            longhorizonduration = Millisecond(Week(settings["horizons"]["long"]["horizonduration_weeks"]))
            longhydroperiodduration = Millisecond(Day(settings["horizons"]["long"]["hydroperiodduration_days"]))
            longrhsdata = getrhsdata(settings["horizons"]["long"]["rhsdata"], datayear, scenarioyearstart, scenarioyearstop)
            longmethod = parse_methods(settings["horizons"]["long"]["rhsmethod"])
            longclusters = settings["horizons"]["long"]["clusters"]
            longunitduration = Millisecond(Hour(settings["horizons"]["long"]["unitduration_hours"]))

            if settings["problems"]["prognosis"]["shrinkable"] == "both"
                longfirstperiod = shorthorizonduration
                longstartafter = longhydroperiodduration + shorthorizonduration
                longshrinkatleast = longhydroperiodduration - steplength
                longminperiod = steplength
                global longhorizon = (longfirstperiod, longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration, longstartafter, longshrinkatleast, longminperiod) # shrinkable
            elseif settings["problems"]["prognosis"]["shrinkable"] == "both_nophasein"
                longfirstperiod = shorthorizonduration
                longstartafter = shorthorizonduration
                longshrinkatleast = longhydroperiodduration - steplength
                longminperiod = steplength
                global longhorizon = (longfirstperiod, longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration, longstartafter, longshrinkatleast, longminperiod) # shrinkable
            elseif settings["problems"]["prognosis"]["shrinkable"] == "no"
                global longhorizon = (longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration)
            end
            lhh, lph = make_horizons(longhorizon...)

            # Simplify modelobjects
            aggzone = getaggzone(settings)
            aggsupplyn = settings["problems"]["prognosis"]["aggsupplyn"]
            removestoragehours = settings["problems"]["shorttermstoragecutoff_hours"]
            residualarealist = settings["problems"]["prognosis"]["residualarealist"]
            simplifyinputs = (aggzone, aggsupplyn, removestoragehours, residualarealist)

            # Medium
            medhorizonduration = Millisecond(Day(settings["horizons"]["med"]["horizonduration_days"]))
            medhydroperiodduration = Millisecond(Day(settings["horizons"]["med"]["hydroperiodduration_days"])); @assert medhorizonduration.value % longhydroperiodduration.value == 0
            medrhsdata = getrhsdata(settings["horizons"]["med"]["rhsdata"], datayear, scenarioyearstart, scenarioyearstop)
            medmethod = parse_methods(settings["horizons"]["med"]["rhsmethod"])
            medclusters = settings["horizons"]["med"]["clusters"]
            medunitduration = Millisecond(Hour(settings["horizons"]["med"]["unitduration_hours"]))

            if (settings["problems"]["prognosis"]["shrinkable"] == "both") || (settings["problems"]["prognosis"]["shrinkable"] == "both_nophasein")
                medfirstperiod = shorthorizonduration
                medstartafter = longstartafter
                medshrinkatleast = longhydroperiodduration - steplength
                medminperiod = steplength
                global medhorizon = (medfirstperiod, medhorizonduration, medhydroperiodduration, medrhsdata, medmethod, medclusters, medunitduration, medstartafter, medshrinkatleast, medminperiod) # shrinkable
            elseif settings["problems"]["prognosis"]["shrinkable"] == "no"
                global medhorizon = (medhorizonduration, medhydroperiodduration, medrhsdata, medmethod, medclusters, medunitduration)
            end
            mhh, mph = make_horizons(medhorizon...)

            # Short
            shorthydroperiodduration = Millisecond(Hour(settings["horizons"]["short"]["hydroperiodduration_hours"])); @assert medhorizonduration.value % shorthorizonduration.value == 0
            shortpowerparts = settings["horizons"]["short"]["powerparts"]
            shorthorizon = (shorthorizonduration, shorthydroperiodduration, shortpowerparts)
            shh, sph = make_horizons(shorthorizon...)

            # Start storages
            dummyprogstorages = getstorages(dummyprogobjects)
            getstartstates!(startstates, settings["problems"], "prognosis", dataset, dummyprogobjects, dummyprogstorages, tnormal)
            startstates_max!(dummyprogstorages, tnormal, startstates)

            # Preallocate storage for problems and results on different cores. Use package DistributedArrays
            # Distribute scenarios
            progscentimes = distribute(progscenmodmethod.scentimes)

            # Problems are built, updated, solved, and stored on a specific core. Moving a problem between cores is expensive, so we want it to only exist on one core. 
            longprobs = distribute([parse_methods(settings["problems"]["prognosis"]["long"]["prob"]) for i in 1:length(progscentimes)], progscentimes)
            medprobs = distribute([parse_methods(settings["problems"]["prognosis"]["med"]["prob"]) for i in 1:length(progscentimes)], progscentimes)
            shortprobs = distribute([parse_methods(settings["problems"]["prognosis"]["short"]["prob"]) for i in 1:length(progscentimes)], progscentimes)

            # Results are moved between cores. These are much smaller than longprobs/medprobs/shortprobs and are inexpensive to move between cores.
            medprices = distribute([Dict() for i in 1:length(progscentimes)], progscentimes)
            shortprices = distribute([Dict() for i in 1:length(progscentimes)], progscentimes)
            medendvaluesobjs = distribute([EndValues() for i in 1:length(progscentimes)], progscentimes)
            nonstoragestates = distribute([Dict{StateVariableInfo, Float64}() for i in 1:length(progscentimes)], progscentimes)

            # Organise inputs and outputs
            probs = (longprobs, medprobs, shortprobs)
            horizons = (lhh, lph, mhh, mph, shh, sph)
            proginput = (numcores, progscentimes, steplength, startstates, simplifyinputs)
            progoutput = (medprices, shortprices, medendvaluesobjs, nonstoragestates)
            
            # Which solver and settings should we use for each problem? Warmstart for long/med and presolve for short
            probmethodsprognosis = [parse_methods(settings["problems"]["prognosis"]["long"]["solver"]), parse_methods(settings["problems"]["prognosis"]["med"]["solver"]), parse_methods(settings["problems"]["prognosis"]["short"]["solver"])]
            # probmethodsprognosis = [CPLEXSimplexMethod(), CPLEXSimplexMethod(), CPLEXSimplexMethod(warmstart=false)]

            # Initialize price prognosis models and run for first time step. Run scenarios in parallell
            @time pl_prognosis_init!(probmethodsprognosis, probs, progelements, horizons, proginput, progoutput)

            # Convert DistributedArray of prices to local process
            medpriceslocal = convert(Vector{Dict}, medprices)
            shortpriceslocal = convert(Vector{Dict}, shortprices)
        end

        println("Mapping between aggregated and detailed storages")
        @time begin
            # Global energy equivalent detailed reservoirs
            for element in elements
                if element.typename == GLOBALENEQKEY
                    enekvglobaldict[split(element.instancename,"GlobalEneq_")[2]] = element.value["Value"]
                end
            end

            medendvaluesdicts = getendvaluesdicts(medendvaluesobjs, detailedrescopl, enekvglobaldict);
        end
    end

    if haskey(settings["problems"], "stochastic")
        println("Init stochastic")
        @time begin
            # Cut parameters
            maxcuts = settings["problems"]["stochastic"]["maxcuts"] # preallocate fixed number of cuts, no cut selection
            lb = settings["problems"]["stochastic"]["lb"] # lower bound of the future value in the first iteration
            reltol = settings["problems"]["stochastic"]["reltol"] # relative tolerance

            # Inputs
            stochasticelements = removeelements!(copy(elements), aggzone=getaggzone(settings), rm_basebalances=!onlysubsystemmodel)
            storageinfo = (startstates, medendvaluesdicts)
            
            ustoragesystemobjects = Tuple{Vector, Vector{Vector}}[]
            ushorts = Bool[]

            # Make modelobjects for short-term subsystemmodels
            if haskey(settings["horizons"]["stochastic"], "short")
                # Parameters for stochastic subsystemmodel problems (could also split totalduration into master- and subduration)
                smpdp = Millisecond(Hour(settings["horizons"]["stochastic"]["short"]["master"]["power"]["periodduration_hours"])) # short/med - master/sub - period duration - power/hydro (commodity)
                smpdh = Millisecond(Hour(settings["horizons"]["stochastic"]["short"]["master"]["hydro"]["periodduration_hours"]))
                sspdp = Millisecond(Hour(settings["horizons"]["stochastic"]["short"]["subs"]["power"]["periodduration_hours"]))
                sspdh = Millisecond(Hour(settings["horizons"]["stochastic"]["short"]["subs"]["hydro"]["periodduration_hours"])) # both master and subproblems for PHS and batteries has 2 hour resolution
                shorttotalduration = Millisecond(Hour(settings["horizons"]["stochastic"]["short"]["horizonduration_hours"])) # total duration of master and subproblem
                shortterminputs = (stochasticelements, shorttotalduration, smpdp, smpdh, sspdp, sspdh, stochscenmodmethod.scentimes, phaseinoffset, shortpriceslocal, true)
                
                # Make sure time resolution of hydro and power are compatible (TODO: Could add function that makes them compatible)
                @assert ceil(Int64, steplength/smpdp) == ceil(Int64, steplength/smpdh)
                @assert ceil(Int64, (shorttotalduration-steplength)/sspdp) == ceil(Int64, (shorttotalduration-steplength)/sspdh)

                @time stochasticmodelobjects = makemastersubobjects!(shortterminputs, ustoragesystemobjects, ushorts, settings)
            end
            
            if haskey(settings["horizons"]["stochastic"], "med") # Make modelobjects for medium-term subsystemmodels
                mmpdp = Millisecond(Hour(settings["horizons"]["stochastic"]["med"]["master"]["power"]["periodduration_hours"]))
                mmpdh = Millisecond(Hour(settings["horizons"]["stochastic"]["med"]["master"]["hydro"]["periodduration_hours"])) # daily resolution in hydro master problems
                mspdp = Millisecond(Hour(settings["horizons"]["stochastic"]["med"]["subs"]["power"]["periodduration_hours"]))
                mspdh = Millisecond(Hour(settings["horizons"]["stochastic"]["med"]["subs"]["hydro"]["periodduration_hours"])) # 7-day resolution in hydro subproblems
                medtotalduration = Millisecond(Day(settings["horizons"]["stochastic"]["med"]["horizonduration_days"])) # we reuse prices for two weeks, so have to be two weeks shorter than price prognosis problem
                medterminputs = (stochasticelements, medtotalduration, mmpdp, mmpdh, mspdp, mspdh, stochscenmodmethod.scentimes, phaseinoffset, medpriceslocal, false)

                @assert ceil(Int64, steplength/mmpdp) == ceil(Int64, steplength/mmpdh)
                @assert ceil(Int64, (medtotalduration-steplength)/mspdp) == ceil(Int64, (medtotalduration-phaseinoffset)/mspdh)

                @time stochasticmodelobjects = makemastersubobjects!(medterminputs, ustoragesystemobjects, ushorts, settings)
            end
            # TODO: Print info about number of short and med term systems

            # Preallocate storagevalues
            if getheadlosscost(settings["problems"]["stochastic"]["master"])
                wwnumscen = stochnumscen*2 + 2 # scenarios + master operative + master operative after headlosscost adjustment
            else
                wwnumscen = stochnumscen*2 + 1 # scenarios + master operative 
            end
            if settings["results"]["storagevalues"]
                ustoragevalues = [zeros(Float64, steps, Int(wwnumscen), length(getcutobjects(msso))) for (i,(msso,ssso)) in enumerate(ustoragesystemobjects)]
            else
                ustoragevalues = [[0] for (i,sso) in enumerate(ustoragesystemobjects)]
            end

            # Add detailed startstates
            stochasticstorages = getstorages(stochasticmodelobjects)
            getstartstates!(startstates, settings["problems"], "stochastic", dataset, stochasticmodelobjects, stochasticstorages, tnormal)
            startstates_max!(stochasticstorages, tnormal, startstates)

            # Distribute subsystemmodels with inputs and outputs on different cores
            if !haskey(dataset, "progelements")
                storagesystemobjects, shorts, storagevalues = distribute_subsystems_flat(ustoragesystemobjects, ushorts, ustoragevalues)
            else
                storagesystemobjects, shorts, storagevalues = distribute_subsystems(ustoragesystemobjects, ushorts, ustoragevalues) # somewhat smart distribution of subsystemmodels to cores based on how many modelobjects in eac subsystemmodel
            end
            masters = distribute([parse_methods(settings["problems"]["stochastic"]["master"]["prob"]) for i in 1:length(storagesystemobjects)], storagesystemobjects)
            subs = distribute([[] for i in 1:length(storagesystemobjects)], storagesystemobjects)
            states = distribute([Dict{StateVariableInfo, Float64}() for i in 1:length(storagesystemobjects)], storagesystemobjects)
            cuts = distribute([SimpleSingleCuts() for i in 1:length(storagesystemobjects)], storagesystemobjects)
            storagesystems = distribute([Dict() for i in 1:length(storagesystemobjects)], storagesystemobjects)

            # Which solver and settings should we use for each problem?
            # probmethodsstochastic = [CPLEXSimplexMethod(), CPLEXSimplexMethod()]
            probmethodsstochastic = [parse_methods(settings["problems"]["stochastic"]["master"]["solver"]), parse_methods(settings["problems"]["stochastic"]["subs"]["solver"])]

            # Initialize subsystemmodel problems and run for first time step. Run subsystemmodels in parallell
            @time pl_stochastic_init!(probmethodsstochastic, numcores, storagesystemobjects, shorts, masters, subs, states, cuts, storagevalues, storageinfo, lb, maxcuts, reltol, tnormal, stochscenmodmethod, settings)

            # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
            @time masterslocal = convert(Vector{Prob}, masters)
            if !haskey(settings["problems"], "clearing")
                @assert length(masterslocal) == 1
                getstartstates!(masterslocal[1], detailedrescopl, enekvglobaldict, startstates)
            end
        end
    end

    if haskey(settings["problems"], "clearing")
        println("Init clearing")
        @time begin
            # Bring data to local core
            @time cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
            @time nonstoragestateslocal = convert(Vector{Dict}, nonstoragestates)
            if settings["results"]["storagevalues"]
                clearingstoragevalues = zeros(Float64, steps, 1, sum([length(cuts.objects) for cuts in cutslocal]))
            end

            # Initialize market clearing problem and run for first time step
            cpdp = Millisecond(Hour(settings["horizons"]["clearing"]["power"]["periodduration_hours"])) # clearing period duration power/battery
            cnpp = ceil(Int64, steplength/cpdp) # clearing numperiods power/battery
            cpdh = Millisecond(Hour(settings["horizons"]["clearing"]["hydro"]["periodduration_hours"])) # clearing period duration hydro
            # cpdh = Millisecond(Hour(2)) # clearing period duration hydro
            cnph = ceil(Int64, steplength/cpdh) # clearing numperiods hydro
            probmethodclearing = parse_methods(settings["problems"]["clearing"]["solver"])
            # probmethodclearing = HighsSimplexSIPMethod(warmstart=false, concurrency=min(8, numcores)) # Which solver and settings should we use for each problem?
            # probmethodclearing = CPLEXIPMMethod(warmstart=false, concurrency=min(8, numcores))
            @time clearing, nonstoragestatesmean, varendperiod = clearing_init(probmethodclearing, elements, tnormal, steplength, cpdp, cpdh, startstates, masterslocal, cutslocal, nonstoragestateslocal, settings)

            # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
            getstartstates!(clearing, detailedrescopl, enekvglobaldict, startstates)

            # Update clearing storage values
            if settings["results"]["storagevalues"]
                j = 0
                for cuts in cutslocal
                    for obj in cuts.objects
                        j += 1
                        balance = getbalance(obj)
                        clearingstoragevalues[1, 1, j] = getcondual(clearing, getid(balance), varendperiod[getid(obj)])
                        if haskey(balance.metadata, GLOBALENEQKEY)
                            clearingstoragevalues[1, 1, j] = clearingstoragevalues[1, 1, j] / balance.metadata[GLOBALENEQKEY]
                        end
                    end
                end
            end
        end
    end

    println("Init results")
    @time begin
        if haskey(settings["results"], "mainresults")
            # Initialize and collect start states
            statenames = collect(keys(startstates))
            statematrix = zeros(length(values(startstates)), Int(steps))
            statematrix[:,1] .= collect(values(startstates))

            if haskey(settings["problems"], "clearing")
                clearingobjects = Dict(zip([getid(obj) for obj in getobjects(clearing)],getobjects(clearing))) # collect results from all areas
                if settings["results"]["mainresults"] == "all"
                    resultobjects = getobjects(clearing) # collect results for all areas
                else
                    resultobjects = getpowerobjects(clearingobjects, settings["results"]["mainresults"]); # only collect results for one area
                end
                prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, rhstermbalances, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages = init_results(steps, clearing, clearingobjects, resultobjects, cnpp, cnph, cpdp, tnormal, true);
            else
                clearingobjects = Dict(zip([getid(obj) for obj in getobjects(masterslocal[1])],getobjects(masterslocal[1]))) # collect results from all areas
                if settings["results"]["mainresults"] == "all"
                    resultobjects = getobjects(masterslocal[1]) # collect results for all areas
                else
                    resultobjects = getpowerobjects(clearingobjects, settings["results"]["mainresults"]); # only collect results for one area
                end
                
                if haskey(settings["horizons"]["stochastic"], "short")
                    cpdp = smpdp
                    cnpp = ceil(Int64, steplength/cpdp)
                    cpdh = smpdh
                    cnph = ceil(Int64, steplength/cpdh)
                else
                    @assert haskey(settings["horizons"]["stochastic"], "med")
                    cpdp = mmpdp
                    cnpp = ceil(Int64, steplength/cpdp)
                    cpdh = mmpdh
                    cnph = ceil(Int64, steplength/cpdh)
                end
                prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, rhstermbalances, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages = init_results(steps, masterslocal[1], clearingobjects, resultobjects, cnpp, cnph, cpdp, tnormal, true);
            end
        end

        # Time problems
        if haskey(settings["problems"], "prognosis")
            prognosistimes = distribute([zeros(steps-1, 3, 3) for i in 1:length(progscentimes)], progscentimes) # update, solve, total per step for long, med, short
        end
        if haskey(settings["problems"], "stochastic")
            stochastictimes = distribute([zeros(steps-1, 9) for i in 1:length(storagesystemobjects)], storagesystemobjects) # update master, update sub, iterate total, solve master, solve sub, iterations, total per system
        end
        if haskey(settings["problems"], "clearing")
            clearingtimes = zeros(steps-1, 3); # update, solve, total
        end
    end

    # Only do scenario modelling and calculate new cuts every 8 days (other reuse scenarios and cuts)
    skipmed = Millisecond(Hour(0))
    skipmax = Millisecond(Hour(steplength*(settings["time"]["skipmax"]-1)))

    stepnr = 2; # already ran first step in initialization

    println("Simulate forward")
    totaltime = @elapsed while stepnr <= steps # while step <= steps and count elapsed time

        # Increment simulation/main scenario and uncertainty scenarios
        tnormal += steplength
        println(tnormal)
    
        increment_scenmodmethod!(simscenmodmethod, phaseinoffset, phaseindelta, phaseinsteps)

        if (prognumscen < simnumscen) && (skipmed.value == 0)
            @time scenariomodelling!(progscenmodmethod, values(dummyprogobjects), prognumscen, simscenmodmethod, progscendelta)
        elseif prognumscen < simnumscen
            increment_scenmodmethod!(progscenmodmethod, phaseinoffset, phaseindelta, phaseinsteps)
        end
        prognumscen != simnumscen && renumber_scenmodmethod!(progscenmodmethod)

        if (stochnumscen < prognumscen) && (skipmed.value == 0)
            # Choose new scenarios
            @time scenariomodelling!(stochscenmodmethod, values(dummyobjects), stochnumscen, progscenmodmethod, stochscendelta)
        elseif stochnumscen < prognumscen
            increment_scenmodmethod!(stochscenmodmethod, phaseinoffset, phaseindelta, phaseinsteps)
        end
    
        # Increment skipmed - should we reuse storagevalues this time step?
        skipmed += Millisecond(steplength)
        if skipmed > skipmax
            skipmed = Millisecond(0)
        end
    
        # Deterministic long/mid/short - calculate scenarioprices for all 30 
        if haskey(settings["problems"], "prognosis")
            progscentimes = distribute(progscenmodmethod.scentimes, progscentimes) # TODO: Find better solution
            @time pl_prognosis!(numcores, longprobs, medprobs, shortprobs, medprices, shortprices, nonstoragestates, startstates, progscentimes, skipmed, prognosistimes, stepnr)
            shortpriceslocal = convert(Vector{Dict}, shortprices)
            if (stochnumscen < prognumscen) && (skipmed.value == 0)
                medpriceslocal = convert(Vector{Dict}, medprices)
                medendvaluesdicts = getendvaluesdicts(medendvaluesobjs, detailedrescopl, enekvglobaldict)
            end
        end

        # Stochastic sub systems - calculate storage value
        if haskey(settings["problems"], "stochastic")
            @time pl_stochastic!(numcores, masters, subs, states, cuts, storagevalues, startstates, medpriceslocal, shortpriceslocal, medendvaluesdicts, shorts, reltol, tnormal, stochscenmodmethod, skipmed, stochastictimes, stepnr, settings)
            masterslocal = convert(Vector{Prob}, masters)
        end
    
        # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
        if haskey(settings["problems"], "clearing")
            # Market clearing
            cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
            nonstoragestateslocal = convert(Vector{Dict}, nonstoragestates)
        
            @time clearing!(clearing, tnormal, startstates, masterslocal, cutslocal, nonstoragestateslocal, nonstoragestatesmean, detailedrescopl, enekvglobaldict, varendperiod, clearingtimes, stepnr, settings)

            if haskey(settings["results"], "mainresults")
                update_results!(stepnr, clearing, prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages, clearingobjects, cnpp, cnph, cpdp, tnormal)   
                statematrix[:,Int(stepnr)] .= collect(values(startstates))
            end

            if settings["results"]["storagevalues"]
                j = 0
                for cuts in cutslocal
                    for obj in cuts.objects
                        j += 1
                        balance = getbalance(obj)
                        clearingstoragevalues[stepnr, 1, j] = getcondual(clearing, getid(balance), varendperiod[getid(obj)])
                        if haskey(balance.metadata, GLOBALENEQKEY)
                            clearingstoragevalues[stepnr, 1, j] = clearingstoragevalues[stepnr, 1, j] / balance.metadata[GLOBALENEQKEY]
                        end
                    end
                end
            end
        else
            getstartstates!(masterslocal[1], detailedrescopl, enekvglobaldict, startstates)

            if haskey(settings["results"], "mainresults")
                update_results!(stepnr, masterslocal[1], prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages, clearingobjects, cnpp, cnph, cpdp, tnormal)   
                statematrix[:,Int(stepnr)] .= collect(values(startstates))
            end
        end
        
        # Increment step
        stepnr += 1
    end
    
    # Total time use and per step
    println(string("The simulation took: ", totaltime/60, " minutes"))
    println(string("Time usage per timestep: ", totaltime/steps, " seconds"))
    
    println("Handle output")
    @time begin
        skipfactor = (skipmax+Millisecond(steplength))/Millisecond(steplength)
        if haskey(settings["problems"], "prognosis") && haskey(settings["problems"], "clearing") # TODO: Split up
            # Prognosis and clearing times
            clearingtimes1 = mean(clearingtimes, dims=1)
            
            factors = [skipfactor,skipfactor,1]
            dims = size(prognosistimes[1])
            dims = (dims..., length(prognosistimes))
            prognosistimes1 = reshape(cat(prognosistimes..., dims=4), dims)
            prognosistimes2 = transpose(dropdims(mean(prognosistimes1,dims=(1,4)),dims=(1,4))).*factors
            progclear = vcat(prognosistimes2, mean(clearingtimes, dims=1))
            df = DataFrame(model=["long","med","short","clearing"], update=progclear[:,1], solve=progclear[:,2], total=progclear[:,3])
            df[!, :other] = df[!, :total] - df[!, :solve] - df[!, :update]
            display(df[!, [1, 2, 3, 5, 4]])
        end
        
        if haskey(settings["problems"], "stochastic")
            # Stochastic times
            core_dists = distribute([0 for i in 1:length(storagesystemobjects)], storagesystemobjects)
            @sync @distributed for core in 1:max(numcores-1,1)
                core_dist = localpart(core_dists)
            
                localix = 0
                for range in localindices(core_dists)
                    for ix in range
                        localix += 1
                        core_dist[localix] = myid()
                    end
                end
            end
            
            dims = size(stochastictimes[1])
            dims = (dims..., length(stochastictimes))
            st1 = reshape(cat(stochastictimes..., dims=4), dims)
            st2 = transpose(dropdims(mean(st1, dims=1), dims=1))
            df = DataFrame(umaster=st2[:,1], usub=st2[:,2], conv=st2[:,3], count=st2[:,4], smaster=st2[:,5], ssub=st2[:,6], hlmaster=st2[:,7], wwres=st2[:,8], total=st2[:,9], short=ushorts, core=core_dists)
            df[df.short .== false, [:umaster, :usub, :conv, :count, :smaster, :ssub, :hlmaster, :wwres, :total]] .= df[df.short .== false, [:umaster, :usub, :conv, :count, :smaster, :ssub, :hlmaster, :wwres, :total]] .* skipfactor
            dfsort = sort(df, :total, rev=true)
            display(dfsort)
            # display(plot(dfsort[!, :total]))
            
            df1 = combine(groupby(df, :core), 
                        :umaster => sum, 
                        :usub => sum, 
                        :conv => sum, 
                        :count => sum, 
                        :smaster => sum, 
                        :ssub => sum,
                        :hlmaster => sum,
                        :wwres => sum,
                        :total => sum)
            dfsort = sort(df1, :total_sum, rev=true)
            display(dfsort)
            # display(plot(dfsort[!, :total_sum]))
            display(mean.(eachcol(select(df1, Not(:core)))))
        end

        # Store results
        data = Dict()

        if haskey(settings["results"], "mainresults")
            # Only keep rhsterms that have at least one value (TODO: Do the same for sypply and demands)
            rhstermtotals = dropdims(sum(rhstermvalues,dims=1),dims=1)
            rhstermsupplyidx = []
            rhstermdemandidx = []

            for k in 1:length(rhsterms)
                if rhstermtotals[k] > 0
                    push!(rhstermsupplyidx, k)
                elseif rhstermtotals[k] < 0
                    push!(rhstermdemandidx, k)
                end
            end

            # Put rhsterms together with supplies and demands
            rhstermsupplyvalues = rhstermvalues[:,rhstermsupplyidx]
            rhstermdemandvalues = rhstermvalues[:,rhstermdemandidx]*-1

            rhstermsupplynames = [getinstancename(rhsterm) for rhsterm in rhsterms[rhstermsupplyidx]]
            rhstermsupplybalancenames = [split(getinstancename(r), "PowerBalance_")[2] for r in rhstermbalances[rhstermsupplyidx]]
            rhstermdemandnames = [getinstancename(rhsterm) for rhsterm in rhsterms[rhstermdemandidx]]
            rhstermdemandbalancenames = [split(getinstancename(r), "PowerBalance_")[2] for r in rhstermbalances[rhstermdemandidx]]

            supplynames = [[getinstancename(plant) for plant in plants];rhstermsupplynames]
            supplybalancenames = [[split(getinstancename(p), "PowerBalance_")[2] for p in plantbalances];rhstermsupplybalancenames]
            supplyvalues = hcat(production,rhstermsupplyvalues)

            demandnames = [[getinstancename(demand) for demand in demands];rhstermdemandnames]
            demandbalancenames = [[split(getinstancename(p), "PowerBalance_")[2] for p in demandbalances];rhstermdemandbalancenames]
            demandvalues = hcat(consumption, rhstermdemandvalues)

            # Prepare for plotting results
            hydronames = [getinstancename(hydro) for hydro in hydrostorages]
            batterynames = [getinstancename(battery) for battery in batterystorages]
            powerbalancenames = [split(getinstancename(getid(powerbalance)), "PowerBalance_")[2] for powerbalance in powerbalances]

            # Convert reservoir filling to TWh
            hydrolevels1 = copy(hydrolevels)
            for (i,hydroname) in enumerate(hydronames)
                if haskey(getbalance(clearingobjects[hydrostorages[i]]).metadata, GLOBALENEQKEY)
                    hydrolevels1[:,i] .= hydrolevels1[:,i]*getbalance(clearingobjects[hydrostorages[i]]).metadata[GLOBALENEQKEY]
                end
            end

            # Indexes
            dim = getoutputindex(config["main"], datayear, scenarioyear)
            x1 = [getisoyearstart(dim) + Week(weekstart-1) + cpdp*(t-1) for t in 1:first(size(supplyvalues))] # power/load resolution
            x2 = [getisoyearstart(dim) + Week(weekstart-1) + cpdh*(t-1) for t in 1:first(size(hydrolevels))]; # reservoir resolution
            x3 = [getisoyearstart(dim) + Week(weekstart-1) + steplength*(t-1) for t in 1:steps]; # state resolution

            outputformat = config["main"]["outputformat"]
            if outputformat != "juliadict"
                datetimeformat = config["main"]["datetimeformat"]
                x1 = Dates.format.(x1, datetimeformat)
                x2 = Dates.format.(x2, datetimeformat)
                x3 = Dates.format.(x3, datetimeformat)
            end

            data["areanames"] = powerbalancenames |> Vector{String}
            data["pricematrix"] = prices
            data["priceindex"] = x1

            data["resnames"] = hydronames
            data["resmatrix"] = hydrolevels1
            data["resindex"] =  x2

            data["batnames"] = batterynames
            data["batmatrix"] = batterylevels
            data["batindex"] =  x2

            data["statenames"] = statenames
            data["statematrix"] = permutedims(statematrix)
            data["stateindex"] =  x3

            data["supplyvalues"] = supplyvalues
            data["supplynames"] = supplynames
            data["supplybalancenames"] = supplybalancenames

            data["demandvalues"] = demandvalues
            data["demandnames"] = demandnames
            data["demandbalancenames"] = demandbalancenames
        end

        if settings["results"]["times"]
            if haskey(settings["problems"], "prognosis") 
                data["prognosistimes"] = prognosistimes1
            end
            if haskey(settings["problems"], "stochastic") 
                data["stochastictimes"] = st1
            end
            if haskey(settings["problems"], "clearing")
                data["clearingtimes"] = clearingtimes
            end
        end

        if settings["results"]["storagevalues"]
            cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
            if haskey(settings["problems"], "clearing")
                data["storagevalues"] = cat(cat(convert(Vector{Array{Float64, 3}}, storagevalues)..., dims=3), clearingstoragevalues, dims=2)
            else
                data["storagevalues"] = cat(convert(Vector{Array{Float64, 3}}, storagevalues)..., dims=3)
            end
            
            data["svindex"] = x3
            data["storagenames"] = [getinstancename(getid(obj)) for cut in cutslocal for obj in cut.objects]
            data["shorts"] = [shorts[i] for (i, cut) in enumerate(cutslocal) for obj in cut.objects]
            data["skipfactor"] = skipfactor
            data["scenarionames"] = String[]
            for i in 1:stochnumscen
                data["scenarionames"] = vcat(data["scenarionames"], [string(i) * " min", string(i) * " max"])
            end
            push!(data["scenarionames"], "Operative master")
            if getheadlosscost(settings["problems"]["stochastic"]["master"])
                push!(data["scenarionames"], "Operative master after")
            end
            if haskey(settings["problems"], "clearing")
                push!(data["scenarionames"], "Operative clearing")
            end
        end
    end

    return data
end