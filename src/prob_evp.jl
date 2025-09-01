"""
Code for end value problem (evp) used in storage valuation (together with stoch)
"""

struct EndValueProblem
    prob::TuLiPa.Prob
    div::Dict
end

function create_evp(scenix::ScenarioIx, subix::SubsystemIx)
    db = get_local_db()
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = get_duration_evp(subsystem)
    modelobjects = make_modelobjects_evp(db.input, db.horizons, subsystem, scenix, startduration, endduration)

    probmethod = parse_methods(settings["problems"]["endvalue"]["solver"])
    prob = TuLiPa.buildprob(probmethod, modelobjects)

    div = Dict()
    div[MainTiming] = zeros(3)

    db.evp[(scenix, subix)] = EndValueProblem(prob, div)

    return
end

function solve_evp(t::TuLiPa.ProbTime, stepnr::Int, skipmed::Millisecond)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_evp
        if core == db.core
            subsystem = db.subsystems[subix]
            evp = db.evp[(scenix, subix)]
            maintiming = evp.div[MainTiming]
            if skipmed_check(subsystem, skipmed)
                maintiming[3] = @elapsed begin
                    # TODO: set nonstorage startstates
                    set_startstates!(evp.prob, TuLiPa.getstorages(TuLiPa.getobjects(evp.prob)), db.startstates)
                    update_prices_evp(stepnr, db.prices_ppp, db.dist_ppp, scenix, subix, evp, subsystem) # TODO: Do not input db
                    update_endstates_evp(db.input, scenix, subix, evp, subsystem) # TODO: Do not input db

                    scentime = get_scentphasein(t, get_scenarios(db.scenmod_sim)[scenix], db.input)
                    maintiming[1] = @elapsed TuLiPa.update!(evp.prob, scentime)
                    # TODO: perform_scenmod(), also for ppp?
                    maintiming[2] = @elapsed TuLiPa.solve!(evp.prob)
                end
            end
        end
    end
    return
end

function skipmed_check(subsystem, skipmed::Millisecond)
    if get_skipmed_impact(subsystem)
        if skipmed.value != 0
            return false
        end
    end
    return true
end

function update_prices_evp(stepnr::Int, prices_ppp::Dict{Tuple{ScenarioIx, TermName, TuLiPa.Id}, Tuple{Int, Vector{Float64}}}, dist_ppp::Vector{Tuple{ScenarioIx, CoreId}}, scenix::ScenarioIx, subix::SubsystemIx, evp::EndValueProblem, subsystem)
    term_ppp = get_horizonterm_evp(subsystem)
    for obj in TuLiPa.getobjects(evp.prob)
        update_prices_obj(prices_ppp, dist_ppp, scenix, stepnr, obj, term_ppp)
    end

    return
end

function update_endstates_evp(input, scenix::ScenarioIx, subix::SubsystemIx, evp::EndValueProblem, subsystem)
    endvaluemethod_ev = get_endvaluemethod_evp(subsystem)

    storages = TuLiPa.getstorages(TuLiPa.getobjects(evp.prob))
    if endvaluemethod_ev == "startequalstop"
        TuLiPa.setendstates!(evp.prob, storages, startstates)
    elseif endvaluemethod_ev == "ppp"
        detailedrescopl = get_dataset(input)["detailedrescopl"]
        enekvglobaldict = get_dataset(input)["enekvglobaldict"]
        for obj in storages
            balance = TuLiPa.getbalance(obj)
            bid = TuLiPa.getid(balance)
            instancename = split(TuLiPa.getinstancename(bid), "Balance_")
            if haskey(enekvglobaldict, instancename[2])
                if enekvglobaldict[instancename[2]] == 0.0
                    continue
                end
            end
            if haskey(detailedrescopl, instancename[2])
                balancename = detailedrescopl[instancename[2]]
                bid = TuLiPa.Id(bid.conceptname, instancename[1] * "Balance_" * balancename * "_hydro_reservoir") # TODO: This should be in the dataset
            end
            endperiod = TuLiPa.gethorizon(TuLiPa.getbalance(obj)).ix_stop
            term_ppp = get_horizonterm_evp(subsystem)
            core_ppp = get_core_ppp(get_local_db().dist_ppp, scenix)
            future = @spawnat core_ppp get_balancedual_ppp(scenix, bid, endperiod, term_ppp)
            dual_ppp = fetch(future)
            if haskey(enekvglobaldict, instancename[2])
                dual_ppp *= enekvglobaldict[instancename[2]]
            end

            TuLiPa.setobjcoeff!(evp.prob, TuLiPa.getid(obj), endperiod, dual_ppp)
        end
    end

    return
end

function get_balancedual_ppp(scenix::ScenarioIx, bid::TuLiPa.Id, period::Int, term_ppp::TermName)
    db = get_local_db()

    ppp = db.ppp[scenix]
    if term_ppp == LongTermName
        return TuLiPa.getcondual(ppp.longprob, bid, period)
    elseif term_ppp == MedTermName
        return TuLiPa.getcondual(ppp.medprob, bid, period)
    elseif term_ppp == ShortTermName
        return TuLiPa.getcondual(ppp.shortprob, bid, period)
    end
end

# Util functions for create_evp() (see also utils for create_mp/create_sp) ----------------------------------------------------------
function make_modelobjects_evp(input, horizons::Dict{Tuple{ScenarioIx, TermName, CommodityName}}, subsystem, scenix::ScenarioIx, startduration::Millisecond, endduration::Millisecond)
    term_ppp = get_horizonterm_evp(subsystem)
    subelements, numperiods_powerhorizon, probhorizons = get_elements_with_horizons(input, horizons, scenix, subsystem, startduration, endduration, term_ppp, false, false)  # stochastic and master is false when evp
    add_scenix_to_InflowParam(subelements, scenix)

    aggzonecopl = get_aggzonecopl(get_aggzone(get_settings(input)))
    keep_hydroramping = has_keephydroramping_evp(get_settings(input))
    change_elements!(subelements, aggzonecopl=aggzonecopl, keep_hydroramping=keep_hydroramping)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = TuLiPa.getmodelobjects(subelements, validate=false)

    return modelobjects
end

