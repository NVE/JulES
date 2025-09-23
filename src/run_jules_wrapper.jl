function getdataset(config, names, filename_clearing, filename_aggregated)
	settings = config[config["main"]["settings"]]

	sti_dataset = joinpath(config["main"]["outputpath"])

	clearing = JulES.JSON.parsefile(joinpath(sti_dataset, filename_clearing))
	clearing = JulES.TuLiPa.getelements(clearing)

	aggregated = JulES.JSON.parsefile(joinpath(sti_dataset, filename_aggregated))
	aggregated = JulES.TuLiPa.getelements(aggregated)

	timevectors = JulES.JSON.parsefile(joinpath(sti_dataset, names["FILENAME_DATAELEMENTS_TIMEVECTORS"] ))
	timevectors = JulES.TuLiPa.getelements(timevectors, sti_dataset)

	elements = vcat(clearing, timevectors)
	elements_ppp = vcat(aggregated, timevectors)
	
	storage_mapping = JulES.JSON.parsefile(
		joinpath(sti_dataset, names["FILENAME_STORAGE_MAPPING"]),
		dicttype=Dict{String, String},
	)

	startmag_aggregated = JulES.JSON.parsefile(
		joinpath(sti_dataset, names["FILENAME_START_STORAGES_AGGREGATED"]),
		dicttype=Dict{String, Float64},
	)
	
	startmag_clearing = JulES.JSON.parsefile(
		joinpath(sti_dataset, names["FILENAME_START_STORAGES_CLEARING"]),
		dicttype=Dict{String, Float64},
	)

	return Dict(
		"elements" => elements,
		"elements_ppp" => elements_ppp,
		"detailedrescopl" => storage_mapping,
		"startmagdict" => startmag_clearing,
		"aggstartmagdict" => startmag_aggregated,
	)
end

function load_ifm_dep()
	mod = @__MODULE__
	@everywhere begin
		@eval $mod using CSV
		@eval $mod using Random
		@eval $mod using OrdinaryDiffEq
		@eval $mod using Lux
		@eval $mod using ComponentArrays
		@eval $mod using Interpolations
		@eval $mod using JLD2
	end
end

function run_jules(
			config_path,
			datayear,
			weatheryear,
			outputpath,
			JulESNames,
			filename_clearing,
			filename_aggregated
			)

        config = YAML.load_file(config_path)

        dataset = getdataset(
			config, 
			JulESNames, 
			filename_clearing, 
			filename_aggregated
		)
		
        input = JulES.DefaultJulESInput(config, dataset, datayear, weatheryear)

		has_ifm_results(input) && load_ifm_dep()
		
		@time data = JulES.run_serial(input)
        println("Total serial time above")

        println("Save output")
        @time h5open(outputpath, "w") do file
            for (k,v) in data
                println(k)
                write(file, k, v)
            end
        end
end
