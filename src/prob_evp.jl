function create_evp(db::LocalDB, scenix::ScenarioIx, subix::SubsystemIx)
    subsystem = get_subsystems(db)[subix]
    settings = get_settings(db)

    startduration = Millisecond(0)
    endduration = get_duration_evp(subsystem)
    modelobjects = make_modelobjects_evp(db, scenix, subix, startduration, endduration)

    probmethod = parse_methods(settings["problems"]["endvalue"]["solver"])
    prob = buildprob(probmethod, modelobjects)

    db.evp[(scenix, subix)] = EndValueProblem(prob, Dict())

    return
end

function solve_evp(t, delta, stepnr, skipmed)
    update_startstates_evp(stepnr, t, skipmed)
    update_endstates_evp(skipmed)
    update_prices_evp(stepnr, skipmed)
    solve_local_evp(t, skipmed)
    # TODO: perform_scenmod()
    return
end

function solve_local_evp(t, skipmed)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_evp
        if core == db.core
            if skipmed_check(subix, skipmed)
                evp = db.evp[(scenix, subix)]
                scentime = get_scentphasein(t, get_scenarios(db.scenmod_evp)[scenix], db.input)
                update!(evp.prob, scentime)
                solve!(evp.prob)
            end
        end
    end
    return
end

function skipmed_check(subix, skipmed)
    subsystem = db.subsystems[subix]
    if skipmed.value == 0
        if get_skipmed_impact(subsystem)        
            return true
        end
    end
    return false
end

function update_prices_evp(stepnr, skipmed)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_evp
        if core == db.core
            if skipmed_check(subix, skipmed)
                subsystem = db.subsystems[subix]
                evp = db.evp[(scenix, subix)]
                duration_evp = get_duration_evp(subsystem)
                for obj in getobjects(evp.prob)
                    update_prices_obj(db, scenix, subix, stepnr, obj, duration_evp)
                end
            end
        end
    end

    return
end

function update_endstates_evp(skipmed)
    db = get_local_db()

    for (scenix, subix, core) in db.dist_evp
        if core == db.core
            if skipmed_check(subix, skipmed)
                evp = db.evp[(scenix, subix)]
                subsystem = get_subsystems(db)[subix]
                endvaluemethod_sevp = get_endvaluemethod_evp(subsystem)

                storages = getstorages(getobjects(evp.prob))
                if endvaluemethod_sp == "startequalstop"
                    setendstates!(sp.prob, storages, startstates)
                elseif endvaluemethod_sp == "ppp"
                    enekvglobaldict = get_dataset(db)["enekvglobaldict"]
                    for obj in storages
                        commodity = getcommodity(getbalance(obj))
                        duration_evp = get_duration_evp(subsystem)
                        term_ppp = get_term_ppp(db, subix, scenix, duration_evp)
                        horizon_ppp = db.horizons[(scenix, term_ppp, commodity)]
                        period = getendperiodfromduration(horizon_ppp, duration_evp)
                        core_ppp = get_core_ppp(db, scenix)
                        bid = getid(getbalance(storage))
                        future = @spawnat core_ppp get_balancedual_ppp(scenix, bid, period, term_ppp)
                        dual_ppp = fetch(future)

                        instancename = getinstancename(getid(obj))
                        if haskey(enekvglobaldict, instancename)
                            dual_ppp *= enekvglobaldict[instancename]
                        end

                        setobjcoeff!(p, getid(obj), period, dual_ppp)
                    end
                end 
            end
        end
    end

    return
end

function get_balancedual_ppp(scenix, bid, period, term_ppp)
    db = get_local_db()

    ppp = db.ppp[scenix]
    if term_ppp == "LongTermName"
        return -getcondual(ppp.longprob, bid, period)
    elseif term_ppp == "MedTermName"
        return -getcondual(ppp.medprob, bid, period)
    elseif term_ppp == "ShortTermName"
        return -getcondual(ppp.shortprob, bid, period)
    end
end

function update_startstates_evp(stepnr, t, skipmed)
    db = get_local_db()
        
    # TODO: set nonstorage startstates
    for (scenix, subix, core) in db.dist_evp
        if core == db.core
            if skipmed_check(subix, skipmed)
                evp = db.evp[(scenix, subix)]
                set_startstates!(evp.prob, get_storages(evp.prob), db.startstates)
            end
        end
    end
end

# Util functions for create_evp() (see also utils for create_mp/create_sp) ----------------------------------------------------------
function make_modelobjects_evp(db, scenix, subix, startduration, endduration)
    subsystem = get_subsystems(db)[subix]
    subelements, numperiods_powerhorizon = get_elements_with_horizons(db, scenix, subix, startduration, endduration)

    aggzonecopl = get_aggzonecopl(get_settings(db.input))
    change_elements!(subelements, aggzonecopl)

    add_prices!(subelements, subsystem, numperiods_powerhorizon, aggzonecopl)

    modelobjects = getmodelobjects(subelements, validate=false)

    return modelobjects
end

