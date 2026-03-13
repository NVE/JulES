function getdataset(config, names, filename_clearing, filename_aggregated)
    settings = config[config["main"]["settings"]]

    sti_dataset = joinpath(config["main"]["outputpath"])

    clearing = JulES.JSON.parsefile(joinpath(sti_dataset, filename_clearing))
    clearing = JulES.TuLiPa.getelements(clearing)

    aggregated = JulES.JSON.parsefile(joinpath(sti_dataset, filename_aggregated))
    aggregated = JulES.TuLiPa.getelements(aggregated)

    timevectors = JulES.JSON.parsefile(joinpath(sti_dataset, names["FILENAME_DATAELEMENTS_TIMEVECTORS"]))
    timevectors = JulES.TuLiPa.getelements(timevectors, sti_dataset)

    elements = vcat(clearing, timevectors)
    elements_ppp = vcat(aggregated, timevectors)

    storage_mapping = JulES.JSON.parsefile(
        joinpath(sti_dataset, names["FILENAME_STORAGE_MAPPING"]),
        dicttype=Dict{String,String},
    )

    startmag_aggregated = JulES.JSON.parsefile(
        joinpath(sti_dataset, names["FILENAME_START_STORAGES_AGGREGATED"]),
        dicttype=Dict{String,Float64},
    )

    startmag_clearing = JulES.JSON.parsefile(
        joinpath(sti_dataset, names["FILENAME_START_STORAGES_CLEARING"]),
        dicttype=Dict{String,Float64},
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
    if myid() == 1
        function ensure_packages(pkgs::Vector{String})
            deps = values(Pkg.dependencies())
            not_installed = filter(pkg -> !any(d -> d.name == pkg, deps), pkgs)
            if !isempty(not_installed)
                @info "Installing missing packages: ", join(not_installed, ", ")
                Pkg.add(not_installed)
            else
                @info "All packages already installed."
            end
        end
        ensure_packages(["OrdinaryDiffEq", "ComponentArrays", "Interpolations", "JLD2"])
    end

    @everywhere begin
        Pkg.instantiate()
        Base.eval(Main, :(using OrdinaryDiffEq))
        Base.eval(Main, :(using ComponentArrays))
        Base.eval(Main, :(using Interpolations))
        Base.eval(Main, :(using JLD2))
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

    @info "Starting JulES run" datayear = datayear weatheryear = weatheryear outputpath = outputpath workers = nworkers()

    config = YAML.load_file(config_path)

    dataset = getdataset(
        config,
        JulESNames,
        filename_clearing,
        filename_aggregated
    )

    input = JulES.DefaultJulESInput(config, dataset, datayear, weatheryear)

    if has_ifm_results(input)
        @info "Loading IFM dependencies"
        load_ifm_dep()
    else
        @debug "IFM dependency loading skipped" reason = "direct mode"
    end

    @debugtime "Run serial" data = JulES.run_serial(input)
    @debugtime "Save output" h5open(outputpath, "w") do file
        for (k, v) in data
            @debug k
            write(file, k, v)
        end
    end
end
