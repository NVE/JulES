module JulES

using DistributedArrays
using TuLiPa
using Dates
using Distributed
include("util.jl") # Various useful functions
include("scenariomodelling.jl") # Code for scenario modelling
include("prognosis.jl") # Code for price prognosis problems
include("stochastic.jl") # Code for stochastic subsystem problems
include("clearing.jl") # Code for market clearing problems
include("run_serial.jl")

using PrecompileTools

@setup_workload begin
	using TuLiPa, Dates, JulES
	using DataFrames, Plots, Statistics, JSON, Distributed, Clustering
	numcores = 1
	datayearstart = 2023
	scenarioyear = 1995
	weekstart = 46
	scenarioyearstart = 1991
	scenarioyearstop = 2021
	totalscen = 30 # scenarios to consider uncertainty for
	simulationyears = 3

	# Standard time for market clearing - perfect information so simple time type
	datatime = getisoyearstart(datayearstart) + Week(weekstart-1)
	tnormal = PrognosisTime(datatime, datatime, getisoyearstart(scenarioyear) + Week(weekstart-1))

	# Phasein settings
	phaseinoffsetdays = 2 # also simulation step length
	phaseinoffset = Millisecond(Day(phaseinoffsetdays)) # phase in straight away from second stage scenarios
	phaseindelta = Millisecond(Week(5)) # Phase in the second stage scenario over 5 weeks
	phaseinsteps = 5 # Phase in second stage scenario in 5 steps

	# Make scenario times for all uncertainty scenarios. List of tuples with tnormal, tphasein and scenarionumber
	totalscentimes = []
	for scen in 1:totalscen
		scentnormal = PrognosisTime(datatime, datatime, getisoyearstart(scenarioyear + scen - 1) + Week(weekstart-1))
		scentphasein = PhaseinPrognosisTime(datatime, datatime, getisoyearstart(scenarioyear) + Week(weekstart-1), getisoyearstart(scenarioyear + scen - 1) + Week(weekstart-1), phaseinoffset, phaseindelta, phaseinsteps);
		push!(totalscentimes, (scentnormal, scentphasein, scen))
	end
	
	elements = TuLiPa.gettestdataset()
	JulES.addscenariotimeperiod_vector!(elements, scenarioyearstart, scenarioyearstop)

	@compile_workload begin
        # Set horizons for price prognosis models
        # All
        shorthorizonduration = Millisecond(Day(7))

		# Long
		longfirstperiod = shorthorizonduration
		longhorizonduration = Millisecond(Week(5*52))
		longhydroperiodduration = Millisecond(Day(7*6))
		longrhsdata = DynamicExogenPriceAHData(Id("Balance", "PowerBalance_TYSKLAND")) # TODO: If dynamic use tphasein
		longmethod = KMeansAHMethod()
		longclusters = 4
		longunitduration = Millisecond(Hour(6))
		longstartafter = longhydroperiodduration + shorthorizonduration
		longshrinkatleast = longhydroperiodduration - phaseinoffset
		longminperiod = phaseinoffset

		longhorizon = (longfirstperiod, longhorizonduration, longhydroperiodduration, longrhsdata, longmethod, longclusters, longunitduration, longstartafter, longshrinkatleast, longminperiod) # shrinkable
		lhh, lph = JulES.make_horizons(longhorizon...)

		# Short
		shorthydroperiodduration = Millisecond(Day(1))
		shortpowerparts = 8
		shorthorizon = (shorthorizonduration, shorthydroperiodduration, shortpowerparts)
		shh, sph = JulES.make_horizons(shorthorizon...)

        # Start storages
        dummyshortobjects, dummyshh, dummysph = make_obj(elements, shh, sph)
		simplify!(dummyshortobjects; aggsupplyn=4, removestoragehours=10)
		addPowerUpperSlack!(dummyshortobjects)
		dummyprob = buildprob(HighsSimplexMethod(), dummyshortobjects)

		# Preallocate storage for problems and results on different cores. Use package DistributedArrays
		# Distribute scenarios
		allscenarios = JulES.distribute(totalscentimes)

		# Problems are built, updated, solved, and stored on a specific core. Moving a problem between cores is expensive, so we want it to only exist on one core. 
		longprobs = JulES.distribute([HiGHS_Prob() for i in 1:length(allscenarios)], allscenarios)
		medprobs = JulES.distribute([HiGHS_Prob() for i in 1:length(allscenarios)], allscenarios)
		shortprobs = JulES.distribute([HiGHS_Prob() for i in 1:length(allscenarios)], allscenarios)

		# Results are moved between cores. These are much smaller than longprobs/medprobs/shortprobs and are inexpensive to move between cores.
		medprices = JulES.distribute([Dict() for i in 1:length(allscenarios)], allscenarios)
		shortprices = JulES.distribute([Dict() for i in 1:length(allscenarios)], allscenarios)
		medendvaluesobjs = JulES.distribute([EndValues() for i in 1:length(allscenarios)], allscenarios)
		nonstoragestates = JulES.distribute([Dict{StateVariableInfo, Float64}() for i in 1:length(allscenarios)], allscenarios)

		# Which solver and settings should we use for each problem? Warmstart for long/med and presolve for short
		probmethodsprognosis = [HighsSimplexMethod(), HighsSimplexMethod(), HighsSimplexMethod(warmstart=false)]
	end
end

end