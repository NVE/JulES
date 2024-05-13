

function solve_inflow_models(t, stepnr)
    db = get_local_db()
    for (core, station) in db.dist_inflow_models
        if core == db.core
            inflow_model = db.inflow_models[station]
            S0 = estimate_S0(inflow_model, getscenariotime(t))
            scenarios = get_scenarios(db.scenmod_ppp)
            (stored_stepnr, d) = db.inflow_prognosis[station]
            for (scenix, scen) in enumerate(scenarios)
                scentime = get_scentphasein(t, scen, db.input)
                d[scenix] = predict(inflow_model, S0, getscenariotime(scentime))
            end
            db.inflow_prognosis[station] = (stepnr, d)
        end
    end
end
