function run_serial(config, datayear, scenarioyear, dataset)
    mainconfig = config["main"]
    settings = config[mainconfig["settings"]]
    numcores = mainconfig["numcores"]

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
        detailedelements = dataset["detailedelements"]
        detailedrescopl = dataset["detailedrescopl"]
        
        addscenariotimeperiod_vector!(elements, scenarioyearstart, scenarioyearstop)
        if !settings["problems"]["stochastic"]["onlyagghydro"]
            addscenariotimeperiod_vector!(detailedelements, scenarioyearstart, scenarioyearstop)
        end
    end

    println("Make dummy objects") # for use in scenario modelling, validate elements and collect storages
    @time begin
        # Horizons are needed to build modelobjects, but not used in scenario modelling
        dummyperiods = 10
        dummyperiodduration = Millisecond(Hour(24))
        power_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)
        hydro_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)

        # Make dummy detailed elements
        dummydetailedobjects, dummydhh, dummydph = make_obj(detailedelements, hydro_horizon, power_horizon, validate=true)

        # Make dummy aggregated elements
        if (elements != []) && !settings["problems"]["stochastic"]["onlyagghydro"]
            dummyprogobjects, dummyphh, dummypph = make_obj(elements, hydro_horizon, power_horizon, validate=true)
        else
            dummyprogobjects = dummydetailedobjects
        end
    end

    println("Init scenario modelling for simulation, prognosis and stochastic") # choose scenarios for the whole simulation
    @time begin
        simnumscen = settings["scenariogeneration"]["simulation"]["numscen"]; @assert simnumscen <= datanumscen
        simscendelta = MsTimeDelta(Day(settings["scenariogeneration"]["simulation"]["scendelta_days"])) # scenario modelling based on the next 3 years, even though the scenario problems can be longer
        if simnumscen == datanumscen
            simscenmodmethod = datascenmodmethod
        else
            global simscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["simulation"], simnumscen)
        end
        @time scenariomodelling!(simscenmodmethod, values(dummyprogobjects), simnumscen, datascenmodmethod, simscendelta) # see JulES/scenariomodelling.jl
        simnumscen != datanumscen && renumber_scenmodmethod!(simscenmodmethod)

        # Prognosis scenario modelling - choose scenarios for the price prognosis models
        prognumscen = settings["scenariogeneration"]["prognosis"]["numscen"]; @assert prognumscen <= simnumscen
        progscendelta = MsTimeDelta(Day(settings["scenariogeneration"]["prognosis"]["scendelta_days"])) # scenario modelling based on the next 3 years, even though the scenario problems can be longer
        if prognumscen == simnumscen
            progscenmodmethod = simscenmodmethod
        else
            global progscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["prognosis"], prognumscen)
        end
        @time scenariomodelling!(progscenmodmethod, values(dummyprogobjects), prognumscen, simscenmodmethod, progscendelta); # see JulES/scenariomodelling.jl
        prognumscen != simnumscen && renumber_scenmodmethod!(progscenmodmethod)

        # Stochastic scenario modelling - choose scenarios for the price stochastic models
        stochnumscen = settings["scenariogeneration"]["stochastic"]["numscen"]; @assert stochnumscen <= prognumscen
        stochscendelta = MsTimeDelta(Day(settings["scenariogeneration"]["stochastic"]["scendelta_days"])) # scenario modelling based on the next 3 years, even though the scenario problems can be longer
        if stochnumscen == prognumscen
            stochscenmodmethod = progscenmodmethod
        else
            global stochscenmodmethod = getscenmodmethod(settings["scenariogeneration"]["stochastic"], stochnumscen)
        end

        @time scenariomodelling!(stochscenmodmethod, values(dummydetailedobjects), stochnumscen, progscenmodmethod, stochscendelta); # see JulES/scenariomodelling.jl
    end

    medendvaluesdicts = Dict[]
    startstates = Dict{String, Float64}()
    if elements != []
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
                longshrinkatleast = longhydroperiodduration - phaseinoffset
                longminperiod = steplength
                global longhorizon = (longfirstperiod, longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration, longstartafter, longshrinkatleast, longminperiod) # shrinkable
            elseif settings["problems"]["prognosis"]["shrinkable"] == "both_nophasein"
                longfirstperiod = shorthorizonduration
                longstartafter = shorthorizonduration
                longshrinkatleast = longhydroperiodduration - phaseinoffset
                longminperiod = steplength
                global longhorizon = (longfirstperiod, longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration, longstartafter, longshrinkatleast, longminperiod) # shrinkable
            elseif settings["problems"]["prognosis"]["shrinkable"] == "no"
                global longhorizon = (longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration)
            end
            lhh, lph = make_horizons(longhorizon...)

            # Simplify modelobjects
            aggzone = settings["problems"]["prognosis"]["aggzone"]
            aggsupplyn = settings["problems"]["prognosis"]["aggsupplyn"]
            removestoragehours = settings["problems"]["prognosis"]["shorttermstoragecutoff_hours"]
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
                medshrinkatleast = longhydroperiodduration - phaseinoffset
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
            dummystorages = getstorages(dummyprogobjects)
            getstartstates!(startstates, settings["problems"]["prognosis"], dataset, dummyprogobjects, dummystorages, tnormal)
            startstates_max!(dummystorages, tnormal, startstates)

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
            proginput = (numcores, progscentimes, phaseinoffset, startstates, simplifyinputs)
            progoutput = (medprices, shortprices, medendvaluesobjs, nonstoragestates)
            
            # Which solver and settings should we use for each problem? Warmstart for long/med and presolve for short
            probmethodsprognosis = [parse_methods(settings["problems"]["prognosis"]["long"]["solver"]), parse_methods(settings["problems"]["prognosis"]["med"]["solver"]), parse_methods(settings["problems"]["prognosis"]["short"]["solver"])]
            # probmethodsprognosis = [CPLEXSimplexMethod(), CPLEXSimplexMethod(), CPLEXSimplexMethod(warmstart=false)]

            # Initialize price prognosis models and run for first time step. Run scenarios in parallell
            @time pl_prognosis_init!(probmethodsprognosis, probs, elements, horizons, proginput, progoutput)
        end

        println("Mapping between aggregated and detailed storages")
        @time begin
            # Global energy equivalent detailed reservoirs
            enekvglobaldict = Dict()
            for element in detailedelements
                if element.typename == GLOBALENEQKEY
                    enekvglobaldict[split(element.instancename,"GlobalEneq_")[2]] = element.value["Value"]
                end
            end

            medendvaluesdicts = getendvaluesdicts(medendvaluesobjs, detailedrescopl, enekvglobaldict);
        end
    end

    println("Init stochastic")
    @time begin
        # Cut parameters
        maxcuts = settings["problems"]["stochastic"]["maxcuts"] # preallocate fixed number of cuts, no cut selection
        lb = settings["problems"]["stochastic"]["lb"] # lower bound of the future value in the first iteration
        reltol = settings["problems"]["stochastic"]["reltol"] # relative tolerance

        # Parameters for stochastic subsystem problems (could also split totalduration into master- and subduration)
        smpdp = Millisecond(Hour(settings["horizons"]["master"]["short"]["power"]["periodduration_hours"])) # short/med - master/sub - period duration - power/hydro (commodity)
        smpdh = Millisecond(Hour(settings["horizons"]["master"]["short"]["hydro"]["periodduration_hours"]))
        sspdp = Millisecond(Hour(settings["horizons"]["subs"]["short"]["power"]["periodduration_hours"]))
        sspdh = Millisecond(Hour(settings["horizons"]["subs"]["short"]["hydro"]["periodduration_hours"])) # both master and subproblems for PHS and batteries has 2 hour resolution
        mmpdp = Millisecond(Hour(settings["horizons"]["master"]["med"]["power"]["periodduration_hours"]))
        mmpdh = Millisecond(Hour(settings["horizons"]["master"]["med"]["hydro"]["periodduration_hours"])) # daily resolution in hydro master problems
        mspdp = Millisecond(Hour(settings["horizons"]["subs"]["med"]["power"]["periodduration_hours"]))
        mspdh = Millisecond(Hour(settings["horizons"]["subs"]["med"]["hydro"]["periodduration_hours"])) # 7-day resolution in hydro subproblems
        shorttotalduration = Millisecond(Hour(settings["horizons"]["short"]["horizonduration_hours"])) # total duration of master and subproblem
        medtotalduration = Millisecond(Day(settings["horizons"]["med"]["horizonduration_days"])) - Millisecond(Day(settings["horizons"]["subs"]["shorterthanprognosismed_days"])) # we reuse prices for two weeks, so have to be two weeks shorter than price prognosis problem

        # Make sure time resolution of hydro and power are compatible (TODO: Could add function that makes them compatible)
        @assert ceil(Int64, phaseinoffset/smpdp) == ceil(Int64, phaseinoffset/smpdh)
        @assert ceil(Int64, (shorttotalduration-phaseinoffset)/sspdp) == ceil(Int64, (shorttotalduration-phaseinoffset)/sspdh)
        @assert ceil(Int64, phaseinoffset/mmpdp) == ceil(Int64, phaseinoffset/mmpdh)
        @assert ceil(Int64, (medtotalduration-phaseinoffset)/mspdp) == ceil(Int64, (medtotalduration-phaseinoffset)/mspdh)

        # Convert DistributedArray of prices to local process
        medpriceslocal = nothing
        shortpriceslocal = nothing
        if elements != []
            medpriceslocal = convert(Vector{Dict}, medprices)
            shortpriceslocal = convert(Vector{Dict}, shortprices)
        end

        # Inputs
        stochasticelements = removeelements!(copy(detailedelements), aggzone=settings["problems"]["prognosis"]["aggzone"])
        storageinfo = (startstates, medendvaluesdicts)
        shortterminputs = (stochasticelements, shorttotalduration, smpdp, smpdh, sspdp, sspdh, stochscenmodmethod.scentimes, phaseinoffset, shortpriceslocal, true)
        medterminputs = (stochasticelements, medtotalduration, mmpdp, mmpdh, mspdp, mspdh, stochscenmodmethod.scentimes, phaseinoffset, medpriceslocal, false)

        ustoragesystemobjects = Tuple{Vector, Vector{Vector}}[]
        ushorts = Bool[]
        # Make modelobjects for short-term subsystems
        @time stochasticmodelobjects = makemastersubobjects!(shortterminputs, ustoragesystemobjects, ushorts)
        # TODO: Print info about number of short and long term systems

        # Make modelobjects for medium-term subsystems
        @time makemastersubobjects!(medterminputs, ustoragesystemobjects, ushorts)

        # Add detailed startstates
        detailedstorages = getstorages(stochasticmodelobjects)
        getstartstates!(startstates, settings["problems"]["stochastic"], dataset, stochasticmodelobjects, detailedstorages, tnormal)
        startstates_max!(detailedstorages, tnormal, startstates)

        # Distribute subsystems with inputs and outputs on different cores
        if settings["problems"]["stochastic"]["onlyagghydro"]
            storagesystemobjects, shorts = distribute_subsystems_flat(ustoragesystemobjects, ushorts)
        else
            storagesystemobjects, shorts = distribute_subsystems(ustoragesystemobjects, ushorts) # somewhat smart distribution of subsystems to cores based on how many modelobjects in eac subsystem
        end
        masters = distribute([parse_methods(settings["problems"]["stochastic"]["master"]["prob"]) for i in 1:length(storagesystemobjects)], storagesystemobjects)
        subs = distribute([[] for i in 1:length(storagesystemobjects)], storagesystemobjects)
        states = distribute([Dict{StateVariableInfo, Float64}() for i in 1:length(storagesystemobjects)], storagesystemobjects)
        cuts = distribute([SimpleSingleCuts() for i in 1:length(storagesystemobjects)], storagesystemobjects)
        storagesystems = distribute([Dict() for i in 1:length(storagesystemobjects)], storagesystemobjects)

        # Which solver and settings should we use for each problem?
        # probmethodsstochastic = [CPLEXSimplexMethod(), CPLEXSimplexMethod()]
        probmethodsstochastic = [parse_methods(settings["problems"]["stochastic"]["master"]["solver"]), parse_methods(settings["problems"]["stochastic"]["subs"]["solver"])]

        # Initialize subsystem problems and run for first time step. Run subsystems in parallell
        @time pl_stochastic_init!(probmethodsstochastic, numcores, storagesystemobjects, shorts, masters, subs, states, cuts, storageinfo, lb, maxcuts, reltol, tnormal, stochscenmodmethod)

        # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
        @time masterslocal = convert(Vector{Prob}, masters)
        if elements == []
            @assert length(masterslocal) == 1
            getstartstates!(masterslocal[1], detailedrescopl, enekvglobaldict, startstates)
        end
    end

    if elements != []
        println("Init clearing")
        @time begin
            # Bring data to local core
            @time cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
            @time nonstoragestateslocal = convert(Vector{Dict}, nonstoragestates)

            # Initialize market clearing problem and run for first time step
            cpdp = Millisecond(Hour(settings["horizons"]["clearing"]["power"]["periodduration_hours"])) # clearing period duration power/battery
            cnpp = ceil(Int64, steplength/cpdp) # clearing numperiods power/battery
            cpdh = Millisecond(Hour(settings["horizons"]["clearing"]["hydro"]["periodduration_hours"])) # clearing period duration hydro
            # cpdh = Millisecond(Hour(2)) # clearing period duration hydro
            cnph = ceil(Int64, phaseinoffset/cpdh) # clearing numperiods hydro
            probmethodclearing = parse_methods(settings["problems"]["clearing"]["solver"])
            # probmethodclearing = HighsSimplexSIPMethod(warmstart=false, concurrency=min(8, numcores)) # Which solver and settings should we use for each problem?
            # probmethodclearing = CPLEXIPMMethod(warmstart=false, concurrency=min(8, numcores))
            @time clearing, nonstoragestatesmean, varendperiod = clearing_init(probmethodclearing, detailedelements, tnormal, phaseinoffset, cpdp, cpdh, startstates, masterslocal, cutslocal, nonstoragestateslocal)

            # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
            getstartstates!(clearing, detailedrescopl, enekvglobaldict, startstates)
        end
    end

    println("Init results")
    @time begin
        # Initialize and collect start states
        statenames = collect(keys(startstates))
        statematrix = zeros(length(values(startstates)), Int(steps))
        statematrix[:,1] .= collect(values(startstates));

        if elements == []
            clearingobjects = Dict(zip([getid(obj) for obj in getobjects(masterslocal[1])],getobjects(masterslocal[1]))) # collect results from all areas
            # resultobjects = getpowerobjects(clearingobjects,["SORLAND"]); # only collect results for one area
            resultobjects = getobjects(masterslocal[1]) # collect results for all areas
            
            prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, rhstermbalances, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages = init_results(steps, masterslocal[1], clearingobjects, resultobjects, cnpp, cnph, cpdp, tnormal, true);
        else
            clearingobjects = Dict(zip([getid(obj) for obj in getobjects(clearing)],getobjects(clearing))) # collect results from all areas
            # resultobjects = getpowerobjects(clearingobjects,["SORLAND"]); # only collect results for one area
            resultobjects = getobjects(clearing) # collect results for all areas
            
            prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, rhstermbalances, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages = init_results(steps, clearing, clearingobjects, resultobjects, cnpp, cnph, cpdp, tnormal, true);
        end

        # Time problems
        prognosistimes = distribute([zeros(steps-1, 3, 3) for i in 1:length(progscentimes)], progscentimes) # update, solve, total per step for long, med, short
        stochastictimes = distribute([zeros(steps-1, 7) for i in 1:length(storagesystemobjects)], storagesystemobjects) # update master, update sub, iterate total, solve master, solve sub, iterations, total per system
        clearingtimes = zeros(steps-1, 3); # update, solve, total
    end

    # Only do scenario modelling and calculate new cuts every 8 days (other reuse scenarios and cuts)
    skipmed = Millisecond(Hour(0))
    skipmax = Millisecond(Hour(steplength*(settings["time"]["skipmax"]-1)))

    stepnr = 2; # already ran first step in initialization

    println("Simulate forward")
    totaltime = @elapsed while stepnr <= steps # while step <= steps and count elapsed time

        # Increment simulation/main scenario and uncertainty scenarios
        tnormal += phaseinoffset
        println(tnormal)
    
        increment_scenmodmethod!(simscenmodmethod, phaseinoffset, phaseindelta, phaseinsteps)
    
        # Increment skipmed - should we reuse watervalues this time step?
        skipmed += Millisecond(phaseinoffset)
        if skipmed > skipmax
            skipmed = Millisecond(0)
        end
    
        # Deterministic long/mid/short - calculate scenarioprices for all 30 
        if (prognumscen < simnumscen) && (skipmed.value == 0)
            @time scenariomodelling!(progscenmodmethod, values(dummyprogobjects), prognumscen, simscenmodmethod, progscendelta)
        elseif prognumscen < simnumscen
            increment_scenmodmethod!(progscenmodmethod, phaseinoffset, phaseindelta, phaseinsteps)
        end
        prognumscen != simnumscen && renumber_scenmodmethod!(progscenmodmethod)

        if elements != []
            progscentimes = distribute(progscenmodmethod.scentimes, progscentimes) # TODO: Find better solution
            @time pl_prognosis!(numcores, longprobs, medprobs, shortprobs, medprices, shortprices, nonstoragestates, startstates, progscentimes, skipmed, prognosistimes, stepnr-1)
        end

        # Stochastic sub systems - calculate storage value    
        if (stochnumscen < prognumscen) && (skipmed.value == 0)
            # Choose new scenarios
            @time scenariomodelling!(stochscenmodmethod, values(dummydetailedobjects), stochnumscen, progscenmodmethod, stochscendelta)
            
            if elements != []
                medpriceslocal = convert(Vector{Dict}, medprices)
                medendvaluesdicts = getendvaluesdicts(medendvaluesobjs, detailedrescopl, enekvglobaldict)
            end
        elseif stochnumscen < prognumscen
            increment_scenmodmethod!(stochscenmodmethod, phaseinoffset, phaseindelta, phaseinsteps)
        end
        if elements != []
            shortpriceslocal = convert(Vector{Dict}, shortprices)
        end
    
        @time pl_stochastic!(numcores, masters, subs, states, cuts, startstates, medpriceslocal, shortpriceslocal, medendvaluesdicts, shorts, reltol, tnormal, stochscenmodmethod, skipmed, stochastictimes, stepnr-1)
    
        # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
        @time masterslocal = convert(Vector{Prob}, masters)
        if elements == []
            getstartstates!(masterslocal[1], detailedrescopl, enekvglobaldict, startstates)

            update_results!(stepnr, masterslocal[1], prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages, clearingobjects, cnpp, cnph, cpdp, tnormal)   
        else
            # Market clearing
            cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
            nonstoragestateslocal = convert(Vector{Dict}, nonstoragestates)
        
            @time clearing!(clearing, tnormal, startstates, masterslocal, cutslocal, nonstoragestateslocal, nonstoragestatesmean, detailedrescopl, enekvglobaldict, varendperiod, clearingtimes, stepnr-1)

            update_results!(stepnr, clearing, prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages, clearingobjects, cnpp, cnph, cpdp, tnormal)   
        end
        statematrix[:,Int(stepnr)] .= collect(values(startstates))
        
        # Increment step
        stepnr += 1
    end
    
    # Total time use and per step
    println(string("The simulation took: ", totaltime/60, " minutes"))
    println(string("Time usage per timestep: ", totaltime/steps, " seconds"))
    
    println("Handle output")
    @time begin
        # Prognosis and clearing times
        clearingtimes1 = mean(clearingtimes, dims=1)
        
        skipfactor = (skipmax+Millisecond(phaseinoffset))/Millisecond(phaseinoffset)
        factors = [skipfactor,skipfactor,1]
        dims = size(prognosistimes[1])
        dims = (dims..., length(prognosistimes))
        prognosistimes1 = reshape(cat(prognosistimes..., dims=4), dims)
        prognosistimes2 = transpose(dropdims(mean(prognosistimes1,dims=(1,4)),dims=(1,4))).*factors
        progclear = vcat(prognosistimes2, mean(clearingtimes, dims=1))
        df = DataFrame(model=["long","med","short","clearing"], update=progclear[:,1], solve=progclear[:,2], total=progclear[:,3])
        df[!, :other] = df[!, :total] - df[!, :solve] - df[!, :update]
        display(df[!, [1, 2, 3, 5, 4]])
        
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
        df = DataFrame(umaster=st2[:,1], usub=st2[:,2], conv=st2[:,3], count=st2[:,4], smaster=st2[:,5], ssub=st2[:,6], total=st2[:,7], short=ushorts, core=core_dists)
        df[df.short .== false, [:umaster, :usub, :conv, :count, :smaster, :ssub, :total]] .= df[df.short .== false, [:umaster, :usub, :conv, :count, :smaster, :ssub, :total]] .* skipfactor
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
                    :total => sum)
        dfsort = sort(df1, :total_sum, rev=true)
        display(dfsort)
        # display(plot(dfsort[!, :total_sum]))
        display(mean.(eachcol(select(df1, Not(:core)))))


        
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

        # Time
        dim = getoutputindex(mainconfig, datayear, scenarioyear)
        x1 = [getisoyearstart(dim) + Week(weekstart-1) + cpdp*(t-1) for t in 1:first(size(supplyvalues))] # power/load resolution
        x2 = [getisoyearstart(dim) + Week(weekstart-1) + cpdh*(t-1) for t in 1:first(size(hydrolevels))]; # reservoir resolution
        x3 = [getisoyearstart(dim) + Week(weekstart-1) + steplength*(t-1) for t in 1:steps]; # state resolution


        # Store results with binary h5 format
        datetimeformat = "yyyy-mm-ddTHH:MM:SS"

        data = Dict()
        data["areanames"] = powerbalancenames |> Vector{String}
        data["pricematrix"] = prices
        data["priceindex"] = Dates.format.(x1, datetimeformat) # not necessary to store as string

        data["resnames"] = hydronames
        data["resmatrix"] = hydrolevels1
        data["resindex"] =  Dates.format.(x2, datetimeformat)

        data["statenames"] = statenames
        data["statematrix"] = permutedims(statematrix)
        data["stateindex"] =  Dates.format.(x3, datetimeformat)

        data["supplyvalues"] = supplyvalues
        data["supplynames"] = supplynames
        data["supplybalancenames"] = supplybalancenames

        data["demandvalues"] = demandvalues
        data["demandnames"] = demandnames
        data["demandbalancenames"] = demandbalancenames

        data["prognosistimes"] = prognosistimes1
        data["stochastictimes"] = st1
        data["clearingtimes"] = clearingtimes
    end

    return data
end