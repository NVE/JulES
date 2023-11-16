module JulesPrognose

using DataFrames, Statistics, JSON, Distributed, Clustering, FileIO, HDF5, CSV

using TuLiPa, Dates
#include(joinpath(dirname(pwd()),raw"TuLiPa/src/TuLiPa.jl"));
include("JulES.jl");    

# Get dictionary with each detailed reservoir and their water value for each scenario
# TODO: Detailed run-of-river reservoirs get water value from aggregated reservoir hydro
function getendvaluesdicts(endvaluesobjs::Any, detailedrescopl::Dict, enekvglobaldict::Dict)
    endvaluesdicts = Dict[]
    for endvaluesobj in endvaluesobjs
        instance = [getinstancename(getid(obj)) for obj in endvaluesobj.objects]
        endvalues = endvaluesobj.values
        endvaluesdict = Dict(instance .=> endvalues)

        for (k,v) in detailedrescopl
            endvaluesdict["Reservoir_" * k] = endvaluesdict["Reservoir_" * v * "_hydro_reservoir"] * enekvglobaldict[k]
        end
        push!(endvaluesdicts, endvaluesdict)
    end
    
    return endvaluesdicts
end

function get_data(prognoser_path, scenarioyear, weekstart)
    
    sti_dataset = joinpath(prognoser_path, "static_input")
    sti_dataset1 = joinpath(prognoser_path, "Uke_$weekstart", "input")

    # Aggregated data
    thermal = getelements(JSON.parsefile(joinpath(sti_dataset, "termisk1.json")), sti_dataset)
    windsol = getelements(JSON.parsefile(joinpath(sti_dataset, "vindsol.json")), sti_dataset)
    cons = getelements(JSON.parsefile(joinpath(sti_dataset, "forbruk5.json")), sti_dataset)
    aggdetd = getelements(JSON.parsefile(joinpath(sti_dataset, "aggdetd2.json")), sti_dataset)

    nuclear = getelements(JSON.parsefile(joinpath(sti_dataset1, "nuclear.json")), sti_dataset1)
    exogen = getelements(JSON.parsefile(joinpath(sti_dataset1, "exogenprices_prognose1.json")), sti_dataset1)
    brenselspriser = getelements(JSON.parsefile(joinpath(sti_dataset1, "brenselspriser.json")), sti_dataset1)
    transm = getelements(JSON.parsefile(joinpath(sti_dataset1, "nett.json")))   
    agginflow = getelements(JSON.parsefile(joinpath(sti_dataset1, "tilsigsprognoseragg$(scenarioyear).json")), sti_dataset1) 

    elements = vcat(exogen, aggdetd, windsol, transm, cons, agginflow, thermal, nuclear, brenselspriser)

    # Detailed data
    detdseries = getelements(JSON.parsefile(joinpath(sti_dataset, "tidsserier_detd.json")), sti_dataset)
    detdstructure = getelements(JSON.parsefile(joinpath(sti_dataset, "dataset_detd.json")))

    inflow = getelements(JSON.parsefile(joinpath(sti_dataset1, "tilsigsprognoser$(scenarioyear).json")), sti_dataset1)

    detailedelements = vcat(exogen,detdseries,detdstructure, windsol,transm,cons,inflow,thermal,nuclear,brenselspriser)

    # Reservoirs
    detailedrescopl = JSON.parsefile(joinpath(sti_dataset, "magasin_elspot.json"))
    startmagdict_json = JSON.parsefile(joinpath(sti_dataset1, "startmagdict.json"))
    startstates = JSON.parsefile(joinpath(sti_dataset1, "aggstartmagdict.json"), dicttype=Dict{String, Float64})

    return elements, startstates, detailedrescopl, startmagdict_json, detailedelements
end

function run(numcores, prognoser_path, datayearstart, weekstart, scenarioyear; simulationyears = 0, steps = 0)

    sti_dataset = joinpath(prognoser_path, "Uke_$(weekstart)")
    sti_output = joinpath(sti_dataset, "output")
    mkpath(sti_output)


    scenarioyearstart = 1991
    scenarioyearstop = 2021

    totalscen = 30 # scenarios to consider uncertainty for

    elements, startstates, detailedrescopl, startmagdict_json, detailedelements = get_data(prognoser_path, scenarioyear, weekstart)

    # Standard time for market clearing - perfect information so simple time type
    datatime = getisoyearstart(datayearstart) + Week(weekstart - 1)
    tnormal = PrognosisTime(datatime, datatime, getisoyearstart(scenarioyear) + Week(weekstart- 1))

    

    # Phasein settings
    phaseinoffsetdays = 2 # also simulation step length
    phaseinoffset = Millisecond(Day(phaseinoffsetdays)) # phase in straight away from second stage scenarios
    phaseindelta = Millisecond(Week(5)) # Phase in the second stage scenario over 5 weeks
    phaseinsteps = 5 # Phase in second stage scenario in 5 steps

    # Make scenario times for all uncertainty scenarios. List of tuples with tnormal, tphasein and scenarionumber
    totalscentimes = []
    for scen in 1:totalscen
        scentnormal = PrognosisTime(datatime, datatime, getisoyearstart(scenarioyear + scen - 1) + Week(weekstart - 1))
        scentphasein = PhaseinPrognosisTime(datatime, datatime, getisoyearstart(scenarioyear) + Week(weekstart - 1), getisoyearstart(scenarioyear + scen - 1) + Week(weekstart - 1), phaseinoffset, phaseindelta, phaseinsteps);
        push!(totalscentimes, (scentnormal, scentphasein, scen))
    end

    # How many time steps to run the simulation for
    if simulationyears != 0
        steps = Int(ceil((getisoyearstart(datayearstart + simulationyears) - getisoyearstart(datayearstart)).value/phaseinoffset.value))
    end

    addscenariotimeperiod_vector!(elements, scenarioyearstart, scenarioyearstop);

    # Set horizons for price prognosis models
    # Long
    longhorizonduration = Millisecond(Week(5*52))
    longhydroperiodduration = Millisecond(Week(6))
    longrhsdata = DynamicExogenPriceAHData(Id("Balance", "PowerBalance_TYSKLAND")) # TODO: If dynamic use tphasein
    longmethod = KMeansAHMethod()
    longclusters = 4
    longunitduration = Millisecond(Hour(6))
    longhorizon = (longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration)

    longobjects, lhh, lph = make_obj(elements, longhorizon...)
    simplify!(longobjects; aggsupplyn=4, removestoragehours=10, residualarealist=[])
    addPowerUpperSlack!(longobjects)

    # Medium
    medhorizonduration = Millisecond(Week(54))
    medhydroperiodduration = Millisecond(Day(7)); @assert medhorizonduration.value % longhydroperiodduration.value == 0
    medrhsdata = DynamicExogenPriceAHData(Id("Balance", "PowerBalance_TYSKLAND"))
    medmethod = KMeansAHMethod()
    medclusters = 4
    medunitduration = Millisecond(Hour(4))
    medhorizon = (medhorizonduration, medhydroperiodduration, medrhsdata, medmethod, medclusters, medunitduration)

    medobjects, mhh, mph = make_obj(elements, medhorizon...)
    simplify!(medobjects; aggsupplyn=4, removestoragehours=10, residualarealist=[])
    addPowerUpperSlack!(medobjects)

    # Short
    shorthorizonduration = Millisecond(Week(1))
    shorthydroperiodduration = Millisecond(Day(1)); @assert medhorizonduration.value % shorthorizonduration.value == 0
    shortpowerparts = 12
    shorthorizon = (shorthorizonduration, shorthydroperiodduration, shortpowerparts)

    shortobjects, shh, sph = make_obj(elements, shorthorizon...)
    simplify!(shortobjects; removestartup=false, residualarealist=[])
    addPowerUpperSlack!(shortobjects)

    # Start storages
    startstates_max!(getstorages(shortobjects), tnormal, startstates)

    # Preallocate storage for problems and results on different cores. Use package DistributedArrays
    # Distribute scenarios
    allscenarios = distribute(totalscentimes)

    # Problems are built, updated, solved, and stored on a specific core. Moving a problem between cores is expensive, so we want it to only exist on one core. 
    longprobs = distribute([HiGHS_Prob() for i in 1:length(allscenarios)], allscenarios)
    medprobs = distribute([HiGHS_Prob() for i in 1:length(allscenarios)], allscenarios)
    shortprobs = distribute([HiGHS_Prob() for i in 1:length(allscenarios)], allscenarios)

    # Results are moved between cores. These are much smaller than longprobs/medprobs/shortprobs and are inexpensive to move between cores.
    medprices = distribute([Dict() for i in 1:length(allscenarios)], allscenarios)
    shortprices = distribute([Dict() for i in 1:length(allscenarios)], allscenarios)
    medendvaluesobjs = distribute([EndValues() for i in 1:length(allscenarios)], allscenarios)
    nonstoragestates = distribute([Dict{StateVariableInfo, Float64}() for i in 1:length(allscenarios)], allscenarios)

    # Organise inputs and outputs
    probs = (longprobs, medprobs, shortprobs)
    objects = (longobjects, medobjects, shortobjects)
    horizons = (lhh, lph, mhh, mph, shh, sph)
    proginput = (numcores, allscenarios, phaseinoffset, startstates)
    progoutput = (medprices, shortprices, medendvaluesobjs, nonstoragestates)
    
    # Which solver and settings should we use for each problem? Warmstart for long/med and presolve for short
    probmethodsprognosis = [HighsSimplexMethod(), HighsSimplexMethod(), HighsSimplexMethod(warmstart=false)]
    # probmethodsprognosis = [CPLEXSimplexMethod(), CPLEXSimplexMethod(), CPLEXSimplexMethod(warmstart=false)]

    # Initialize price prognosis models and run for first time step. Run scenarios in parallell
    @time pl_prognosis_init!(probmethodsprognosis, probs, objects, horizons, proginput, progoutput)

    #=
    scenario = 1
    index = collect(medprices[scenario]["steprange"])
    prices = medprices[scenario]["matrix"]
    labels = [name for name in medprices[scenario]["names"]]
    p = plot(index,prices*100,label=reshape(labels, 1, length(labels)),legend=:outertopright)

    for scenario in 2:5 # length(allscenarios)
        prices = medprices[scenario]["matrix"]
        labels = [name for name in medprices[scenario]["names"]]
        plot!(p,index,prices*100,label=reshape(labels, 1, length(labels)),legend=:outertopright)
    end
    display(p)

    scenario = 1
    index = collect(shortprices[scenario]["steprange"])
    prices = shortprices[scenario]["matrix"]
    labels = [name for name in shortprices[scenario]["names"]]
    p = plot(index,prices,label=reshape(labels, 1, length(labels)), ylabel="€/MWh",legend=:outertopright)

    for scenario in 2:5 # length(allscenarios)
        prices = shortprices[scenario]["matrix"]
        labels = [name for name in shortprices[scenario]["names"]]
        plot!(p, index, prices*100 ,label=reshape(labels, 1, length(labels)),legend=:outertopright)
    end
    display(p)
    =#

   
    addscenariotimeperiod_vector!(detailedelements, scenarioyearstart, scenarioyearstop);

    # Mapping between aggregated and detailed storages
    

    # Global energy equivalent detailed reservoirs
    enekvglobaldict = Dict()
    for element in detailedelements
        if element.typename == GLOBALENEQKEY
            enekvglobaldict[split(element.instancename,"GlobalEneq_")[2]] = element.value["Value"]
        end
    end

    # Detailed dataset has reservoirs for SE4, aggregated does not, TODO: Improve aggregation/mapping
    for k in keys(detailedrescopl)
        if detailedrescopl[k] == "SVER-SE4"
            detailedrescopl[k] = "SVER-SE3"
        end
    end


    medendvaluesdicts = getendvaluesdicts(medendvaluesobjs, detailedrescopl, enekvglobaldict);

    # Scenario reduction to this amount
    numscen = 7

    # Modelobjects that can be used to reduce scenarios
    scenarioelements = copy(detailedelements)

    # Horizons are needed to build modelobjects, but not used in scenario modelling
    dummyperiods = 10
    dummyperiodduration = Millisecond(Hour(24))
    power_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)
    hydro_horizon = SequentialHorizon(dummyperiods, dummyperiodduration)

    set_horizon!(scenarioelements, "Power", power_horizon)
    set_horizon!(scenarioelements, "Battery", power_horizon)
    set_horizon!(scenarioelements, "Hydro", hydro_horizon)

    scenarioobjects = collect(values(getmodelobjects(scenarioelements)))

    # Scenario modelling method
    scendelta = MsTimeDelta(Day(364)) # scenario modelling based on the next year, even though the scenario problems can be longer
    if numscen >= totalscen
        global scenmodmethod = NoScenarioModellingMethod(totalscen, totalscentimes)
    else
        parts = 4 # divide scendelta into this many parts, calculate sum inflow for each part of the inflow series, then use clustering algorithm
        global scenmodmethod = InflowClusteringMethod(numscen, parts)
    end
    @time scenariomodelling!(scenmodmethod, scenarioobjects, numscen, totalscentimes, scendelta); # see JulES/scenariomodelling.jl



    # Cut parameters
    maxcuts = 13 # preallocate fixed number of cuts, no cut selection
    lb = -1e5 # lower bound of the future value in the first iteration
    reltol = 0.0001 # relative tolerance

    # Parameters for stochastic subsystem problems (could also split totalduration into master- and subduration)
    smpdp = Millisecond(Hour(2)) # short/med - master/sub - period duration - power/hydro (commodity)
    smpdh = Millisecond(Hour(2))
    sspdp = Millisecond(Hour(2))
    sspdh = Millisecond(Hour(2)) # both master and subproblems for PHS and batteries has 2 hour resolution
    mmpdp = Millisecond(Hour(24))
    mmpdh = Millisecond(Hour(24)) # daily resolution in hydro master problems
    mspdp = Millisecond(Hour(168))
    mspdh = Millisecond(Hour(168)) # weekly resolution in hydro subproblems
    shorttotalduration = shorthorizonduration # total duration of master and subproblem
    medtotalduration = medhorizonduration - Millisecond(Week(2)) # we reuse prices for two weeks, so have to be two weeks shorter than price prognosis problem

    # Make sure time resolution of hydro and power are compatible (TODO: Could add function that makes them compatible)
    @assert ceil(Int64, phaseinoffset/smpdp) == ceil(Int64, phaseinoffset/smpdh)
    @assert ceil(Int64, (shorttotalduration-phaseinoffset)/sspdp) == ceil(Int64, (shorttotalduration-phaseinoffset)/sspdh)
    @assert ceil(Int64, phaseinoffset/mmpdp) == ceil(Int64, phaseinoffset/mmpdh)
    @assert ceil(Int64, (medtotalduration-phaseinoffset)/mspdp) == ceil(Int64, (medtotalduration-phaseinoffset)/mspdh)

    # Convert DistributedArray of prices to local process
    medpriceslocal = convert(Vector{Dict}, medprices)
    shortpriceslocal = convert(Vector{Dict}, shortprices)

    # Inputs
    stochasticelements = removeelements!(copy(detailedelements))
    storageinfo = (startstates, medendvaluesdicts)
    shortterminputs = (stochasticelements, shorttotalduration, smpdp, smpdh, sspdp, sspdh, scenmodmethod.scentimes, phaseinoffset, shortpriceslocal, true)
    medterminputs = (stochasticelements, medtotalduration, mmpdp, mmpdh, mspdp, mspdh, scenmodmethod.scentimes, phaseinoffset, medpriceslocal, false)

    ustoragesystemobjects = Tuple{Vector, Vector{Vector}}[]
    ushorts = Bool[]
    # Make modelobjects for short-term subsystems
    @time stochasticmodelobjects = makemastersubobjects!(shortterminputs, ustoragesystemobjects, ushorts)
    # Make modelobjects for medium-term subsystems
    @time makemastersubobjects!(medterminputs, ustoragesystemobjects, ushorts)

    # Add detailed startstates
    merge!(startstates, startmagdict_json) # also read detailed startstates used in the other problems
    startstates_max!(getstorages(stochasticmodelobjects), tnormal, startstates)

    # Distribute subsystems with inputs and outputs on different cores
    storagesystemobjects, shorts = distribute_subsystems(ustoragesystemobjects, ushorts) # somewhat smart distribution of subsystems to cores based on how many modelobjects in eac subsystem
    masters = distribute([HiGHS_Prob() for i in 1:length(storagesystemobjects)], storagesystemobjects)
    subs = distribute([[] for i in 1:length(storagesystemobjects)], storagesystemobjects)
    states = distribute([Dict{StateVariableInfo, Float64}() for i in 1:length(storagesystemobjects)], storagesystemobjects)
    cuts = distribute([SimpleSingleCuts() for i in 1:length(storagesystemobjects)], storagesystemobjects)
    storagesystems = distribute([Dict() for i in 1:length(storagesystemobjects)], storagesystemobjects)

    # Which solver and settings should we use for each problem?
    # probmethodsstochastic = [CPLEXSimplexMethod(), CPLEXSimplexMethod()]
    probmethodsstochastic = [HighsSimplexMethod(), HighsSimplexMethod()]

    # Initialize subsystem problems and run for first time step. Run subsystems in parallell
    @time pl_stochastic_init!(probmethodsstochastic, numcores, storagesystemobjects, shorts, masters, subs, states, cuts, storageinfo, lb, maxcuts, reltol, scenmodmethod.scentimes)

    # Bring data to local core
    masterslocal = convert(Vector{Prob}, masters)
    cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
    nonstoragestateslocal = convert(Vector{Dict}, nonstoragestates)

    # Initialize market clearing problem and run for first time step
    cpdp = Millisecond(Hour(2)) # clearing period duration power/battery
    cnpp = ceil(Int64, phaseinoffset/cpdp) # clearing numperiods power/battery
    cpdh = Millisecond(Hour(6)) # clearing period duration hydro
    # cpdh = Millisecond(Hour(2)) # clearing period duration hydro
    cnph = ceil(Int64, phaseinoffset/cpdh) # clearing numperiods hydro
    probmethodclearing = HighsSimplexMethod(warmstart=false) # Which solver and settings should we use for each problem?
    # probmethodclearing = CPLEXIPMMethod(warmstart=false, concurrency=min(8, numcores))
    @time clearing, nonstoragestatesmean, varendperiod = clearing_init(probmethodclearing, detailedelements, tnormal, phaseinoffset, cpdp, cpdh, startstates, masterslocal, cutslocal, nonstoragestateslocal);

    # Update start states for next time step, also mapping to aggregated storages and max capacity in aggregated
    getstartstates!(clearing, detailedrescopl, enekvglobaldict, startstates)

    # Initialize and collect prices and start states
    price = Dict()
    powerhorizonix = argmax(getnumperiods(h) for h in clearing.horizons)
    getareaprices!(price, clearing, clearing.horizons[powerhorizonix], tnormal)
    areanames = price["names"]

    ix = Vector{DateTime}(undef,Int(length(price["steprange"])*steps))
    ix[1:length(price["steprange"])] .= price["steprange"]

    (pricex,pricey) = size(price["matrix"])
    pricematrix = zeros(Int(pricex*steps),pricey)
    pricematrix[1:pricex,:] .= price["matrix"]

    statenames = collect(keys(startstates))
    statematrix = zeros(length(values(startstates)), Int(steps))
    statematrix[:,1] .= collect(values(startstates));

    clearingobjects = Dict(zip([getid(obj) for obj in clearing.objects],clearing.objects)) # collect results from all areas
    # resultobjects = getpowerobjects(clearingobjects,["SORLAND"]); # only collect results for one area
    resultobjects = clearing.objects # collect results for all areas

    prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, rhstermbalances, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages = init_results(steps, clearing, clearingobjects, resultobjects, cnpp, cnph, cpdp, tnormal, true);

    # Only do scenario modelling and calculate new cuts every 8 days (other reuse scenarios and cuts)
    skipmed = Millisecond(Day(6))
    skipmax = Millisecond(Day(6))

    stepnr = 2; # already ran first step in initialization

    totaltime = @elapsed while stepnr <= steps # while step <= steps and count elapsed time

        # Increment simulation/main scenario and uncertainty scenarios
        tnormal += phaseinoffset

        for i in 1:length(totalscentimes)
            (scentnormal, scentphasein, scenario) = totalscentimes[i]
            scentnormal += phaseinoffset
            scentphasein = PhaseinPrognosisTime(getdatatime(scentnormal), getdatatime(scentnormal), getscenariotime(scentnormal), getscenariotime(scentnormal), phaseinoffset, phaseindelta, phaseinsteps)
            totalscentimes[i] = (scentnormal, scentphasein, scenario)
        end

        # Increment skipmed - should we reuse watervalues this time step?
        skipmed += Millisecond(phaseinoffset)
        if skipmed > skipmax
            skipmed = Millisecond(0)
        end

        # Print here to avoid wrong order
        println(tnormal)

        # Deterministic long/mid/short - calculate scenarioprices for all 30 scenarios
        allscenarios = distribute(totalscentimes, allscenarios) # TODO: Find better solution
        @time pl_prognosis!(numcores, longprobs, medprobs, shortprobs, medprices, shortprices, nonstoragestates, startstates, allscenarios, skipmed)

        # Stochastic sub systems - calculate storage value    
        if skipmed.value == 0
            # Choose new scenarios
            @time scenariomodelling!(scenmodmethod, scenarioobjects, numscen, totalscentimes, scendelta)
            
            medpriceslocal = convert(Vector{Dict}, medprices)
            medendvaluesdicts = getendvaluesdicts(medendvaluesobjs, detailedrescopl, enekvglobaldict)
        else
            # Increment existing scenarios
            for i in 1:length(scenmodmethod.scentimes)
                (scentnormal, scentphasein, scenario) = scenmodmethod.scentimes[i]
                scentnormal += phaseinoffset
                scentphasein = PhaseinPrognosisTime(getdatatime(scentnormal), getdatatime(scentnormal), getscenariotime(scentnormal), getscenariotime(scentnormal), phaseinoffset, phaseindelta, phaseinsteps)
                scenmodmethod.scentimes[i] = (scentnormal, scentphasein, scenario)
            end
        end
        shortpriceslocal = convert(Vector{Dict}, shortprices)

        @time pl_stochastic!(numcores, masters, subs, states, cuts, startstates, medpriceslocal, shortpriceslocal, medendvaluesdicts, shorts, reltol, scenmodmethod.scentimes, skipmed)

        # Market clearing
        masterslocal = convert(Vector{Prob}, masters)
        cutslocal = convert(Vector{SimpleSingleCuts}, cuts)
        nonstoragestateslocal = convert(Vector{Dict}, nonstoragestates)

        @time clearing!(clearing, tnormal, startstates, masterslocal, cutslocal, nonstoragestateslocal, nonstoragestatesmean, detailedrescopl, enekvglobaldict, varendperiod)
        
        # Results
        updateareaprices!(price, clearing, clearing.horizons[powerhorizonix], tnormal)
        ix[Int(length(price["steprange"])*(stepnr-1)+1):Int(length(price["steprange"])*stepnr)] .= price["steprange"]
        pricematrix[Int(pricex*(stepnr-1)+1):Int(pricex*(stepnr)),:] .= price["matrix"]
        statematrix[:,Int(stepnr)] .= collect(values(startstates))
        
        update_results!(stepnr, clearing, prices, rhstermvalues, production, consumption, hydrolevels, batterylevels, powerbalances, rhsterms, plants, plantbalances, plantarrows, demands, demandbalances, demandarrows, hydrostorages, batterystorages, clearingobjects, cnpp, cnph, cpdp, tnormal)   
        
        # Increment step
        stepnr += 1
    end

    println(string("The simulation took: ", totaltime/60, " minutes"))
    println(string("Time usage per timestep: ", totaltime/steps, " seconds"))

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

    hydronames = [getinstancename(hydro) for hydro in hydrostorages]
    powerbalancenames = [split(getinstancename(getid(powerbalance)), "PowerBalance_")[2] for powerbalance in powerbalances]

    # Convert reservoir filling to TWh
    hydrolevels1 = copy(hydrolevels)
    for (i,hydroname) in enumerate(hydronames)
        if haskey(getbalance(clearingobjects[hydrostorages[i]]).metadata, GLOBALENEQKEY)
            hydrolevels1[:,i] .= hydrolevels1[:,i]*getbalance(clearingobjects[hydrostorages[i]]).metadata[GLOBALENEQKEY]
        end
    end

    # Time
    x1 = [getisoyearstart(datayearstart) + Week(weekstart-1) + cpdp*(t-1) for t in 1:first(size(supplyvalues))] # power/load resolution
    x2 = [getisoyearstart(datayearstart) + Week(weekstart-1) + cpdh*(t-1) for t in 1:first(size(hydrolevels))]; # reservoir resolution
    x3 = [getisoyearstart(datayearstart) + Week(weekstart-1) + phaseinoffset*(t-1) for t in 1:steps]; # state resolution




    # Store results with binary h5 format
    datetimeformat = "yyyy-mm-ddTHH:MM:SS"
    modelname = scenarioyear

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

    @time h5open(joinpath(sti_output, "$modelname.h5"), "w") do file
        for (k,v) in data
            display(k)
            write(file, k, v)
        end
    end

    # Store as CSV
    areaprices = rename!(DataFrame(prices, :auto),powerbalancenames)
    areaprices[!,:time] = x1
    CSV.write(joinpath(sti_output, "price$scenarioyear.csv"), areaprices)

    demand = rename!(DataFrame(demandvalues, :auto),demandnames)
    demand[!,:time] = x1
    demand = stack(demand,Not(:time))
    demandcopl = DataFrame(variable=demandnames, area=demandbalancenames)
    demand = leftjoin(demand, demandcopl, on=:variable)
    CSV.write(joinpath(sti_output, "demand$scenarioyear.csv"), demand)

    supply = rename!(DataFrame(supplyvalues, :auto),supplynames)
    supply[!,:time] = x1
    supply = stack(supply,Not(:time))
    supplycopl = DataFrame(variable=supplynames, area=supplybalancenames)
    supply = leftjoin(supply, supplycopl, on=:variable)
    CSV.write(joinpath(sti_output, "supply$scenarioyear.csv"), supply)

    hydro = rename!(DataFrame(hydrolevels, :auto),hydronames)
    hydro[!,:time] = x2
    CSV.write(joinpath(sti_output, "hydro$scenarioyear.csv"), hydro)
    return
end

end