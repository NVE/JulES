function create_evp(db::LocalDB, scenix::ScenarioIx, subix::SubsystemIx)
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = get_duration_evp(subsystem)
    modelobjects = make_modelobjects_evp(db, scenix, subix, startduration, endduration)

    probmethod = parse_methods(settings["problems"]["endvalue"]["solver"])
    prob = TuLiPa.buildprob(probmethod, modelobjects)

    div = Dict()
    div[MainTiming] = zeros(3)

    db.evp[(scenix, subix)] = EndValueProblem(prob, div)

    return
end

function solve_evp(t, stepnr, skipmed)
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
                    update_prices_evp(stepnr, skipmed, db, scenix, subix, evp, subsystem) # TODO: Do not input db
                    update_endstates_evp(skipmed, db, scenix, subix, evp, subsystem) # TODO: Do not input db

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

function skipmed_check(subsystem, skipmed)
    if get_skipmed_impact(subsystem)
        if skipmed.value != 0
            return false
        end
    end
    return true
end

function update_prices_evp(stepnr, skipmed, db, scenix, subix, evp, subsystem)
    term_ppp = get_horizonterm_evp(subsystem)
    for obj in TuLiPa.getobjects(evp.prob)
        update_prices_obj(db, scenix, subix, stepnr, obj, term_ppp)
    end

    return
end

function update_endstates_evp(skipmed, db, scenix, subix, evp, subsystem)
    endvaluemethod_ev = get_endvaluemethod_evp(subsystem)

    storages = TuLiPa.getstorages(TuLiPa.getobjects(evp.prob))
    if endvaluemethod_ev == "startequalstop"
        TuLiPa.setendstates!(evp.prob, storages, startstates)
    elseif endvaluemethod_ev == "ppp"
        detailedrescopl = get_dataset(db)["detailedrescopl"]
        enekvglobaldict = get_dataset(db)["enekvglobaldict"]
        for obj in storages
            balance = TuLiPa.getbalance(obj)
            bid = TuLiPa.getid(balance)
            instancename = split(TuLiPa.getinstancename(bid), "Balance_")
            if haskey(detailedrescopl, instancename[2])
                balancename = detailedrescopl[instancename[2]]
                bid = TuLiPa.Id(bid.conceptname, instancename[1] * "Balance_" * balancename * "_hydro_reservoir") # TODO: This should be in the dataset
            end
            endperiod = TuLiPa.gethorizon(TuLiPa.getbalance(obj)).ix_stop
            term_ppp = get_horizonterm_evp(subsystem)
            core_ppp = get_core_ppp(db, scenix)
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

function get_balancedual_ppp(scenix, bid, period, term_ppp)
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
function make_modelobjects_evp(db, scenix, subix, startduration, endduration)
    subsystem = get_subsystems(db)[subix]
    term_ppp = get_horizonterm_evp(subsystem)
    subelements, numperiods_powerhorizon, horizons = get_elements_with_horizons(db, scenix, subsystem, startduration, endduration, term_ppp, false, false)  # stochastic and master is false when evp
    add_scenix_to_InflowParam(subelements, scenix)

    aggzonecopl = get_aggzonecopl(get_aggzone(get_settings(db.input)))
    change_elements!(subelements, aggzonecopl=aggzonecopl)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = TuLiPa.getmodelobjects(subelements, validate=false)

    return modelobjects
end

