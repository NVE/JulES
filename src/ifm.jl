# TODO: Illustrate how data_obs, data_pred and data_forecast are updated over time

# TODO: Interfaces
abstract type AbstractTwoStateIfmDataUpdater end
abstract type AbstractTwoStateIfmPredictor end

struct TwoStateIfmData
    P::Vector{Float32}   
    T::Vector{Float32}  
    Lday::Vector{Float32} 
    timepoints::Vector{Float32}

    function TwoStateIfmData(P::Vector{Float32}, T::Vector{Float32}, Lday::Vector{Float32})
        ndays = length(P)
        @assert ndays == length(T) == length(Lday)
        timepoints = Float32.(1:ndays)
        return new(P, T, Lday, timepoints)
    end

    function TwoStateIfmData(ndays::Int)
        @assert ndays > 0
        P = zeros(Float32, ndays)
        T = zeros(Float32, ndays)
        Lday = zeros(Float32, ndays)
        timepoints = Float32.(1:ndays)
        return new(P, T, Lday, timepoints)
    end
end

length(x::TwoStateIfmData) = length(x.P)

const ONEDAY_MS_TIMEDELTA = MSTimeDelta(Day(1))

mutable struct TwoStateIfmHandler{P <: AbstractTwoStateIfmPredictor, 
                                U <: AbstractTwoStateIfmDataUpdater, 
                                T1 <: TimeVector, T2 <: TimeVector, T3 <: TimeVector}
    predictor::P

    updater::U

    basin_area::Float64
    m3s_per_mm::Float64

    hist_P::T1
    hist_T::T2
    hist_Lday::T3

    ndays_pred::Int
    ndays_obs::Int
    ndays_forecast::Int

    data_pred::TwoStateIfmData
    data_obs::Union{TwoStateIfmData, Nothing}
    data_forecast::Union{TwoStateIfmData, Nothing}

    prev_t::Union{ProbTime, Nothing}
    ndays_forecast_used::Int

    function TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
            ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        @assert ndays_forecast >= 0
        @assert ndays_pred >= ndays_forecast
        @assert ndays_obs > 0
        isnothing(data_obs) || @assert ndays_obs == length(data_obs)
        isnothing(data_forecast) || @assert ndays_forecast == length(data_forecast)
        (ndays_forecast > 0) && @assert !isnothing(data_forecast)
        data_pred = TwoStateIfmData(ndays_pred)
        m3s_per_mm = 1/((1000**3)/(basin_area*10**6)*86400)
        P = typeof(predictor)
        U = typeof(updater)
        T1 = typeof(hist_P)
        T2 = typeof(hist_T)
        T3 = typeof(hist_Lday)
        return new{P, U, T1, T2, T3}(predictor, updater, basin_area, m3s_per_mm, hist_P, hist_T, hist_Lday, 
            ndays_pred, ndays_obs, ndays_forecast, data_pred, 
            data_obs, data_forecast, nothing, 0)
    end
end

# TODO: Use @inbounds in for loops

function _initial_data_obs_update(m::TwoStateIfmHandler, t::ProbTime)
    if isnothing(m.data_obs)
        m.data_obs = TwoStateIfmData(m.ndays_obs)
        for i in 1:m.ndays_obs
            ndays_back = m.ndays_obs + 1 - i
            start = getdatatime(t) - Day(ndays_back)
            m.data_obs.P[i] = getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.T[i] = getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.Lday[i] = getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
        end
    end
end

function _data_obs_update(m::TwoStateIfmHandler, t::ProbTime)
    # calc ndays_update
    diff_ndays_datatime = Day(getdatatime(t) - getdatatime(m.prev_t)).value
    diff_ndays_scentime = Day(getscenariotime(t) - getscenariotime(m.prev_t)).value
    diff_ndays_max = max(diff_ndays_datatime, diff_ndays_scentime)
    ndays_update = min(m.ndays_obs, diff_ndays_max)
    
    # shift backwards to make room for ndays_update new values
    i = m.ndays_obs - ndays_update + 1
    j = m.ndays_obs - ndays_update * 2
    for __ in 1:ndays_update
        m.data_obs.P[j] = m.data_obs.P[i]
        m.data_obs.T[j] = m.data_obs.T[i]
        m.data_obs.Lday[j] = m.data_obs.Lday[i]
        i += 1
        j += 1
    end

    # use forecast for P and T if available
    m.ndays_forecast_used = 0
    if m.ndays_forecast > 0
        startix = length(m.data_forecast) - m.ndays_forecast + 1
        stopix = length(m.data_forecast) - m.ndays_forecast + ndays_update
        stopix = min(stopix, length(m.data_forecast))
        m.ndays_forecast_used = (stopix - startix) + 1
        j = m.ndays_obs - ndays_forecast_used + 1
        for i in startix:stopix
            m.data_obs.P[j] = m.data_forecast.P[i]
            m.data_obs.T[j] = m.data_forecast.T[i]
            j += 1
        end
    end

    # update possible remaining values for P and T
    ndays_remaining = ndays_update - m.ndays_forecast_used
    if ndays_remaining > 0
        start = getscenariotime(t) - Day(ndays_remaining + 1)
        i = m.ndays_obs - ndays_remaining + 1
        for __ in 1:ndays_remaining
            m.data_obs.P[i] = getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.T[i] = getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
            start += Day(1) 
            i += 1
        end
    end

    # always use hist to update Lday
    start = getscenariotime(t) - Day(ndays_update + 1)
    i = m.ndays_obs - ndays_update + 1
    for __ in 1:ndays_update
        m.data_obs.Lday[i] = getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
        start += Day(1) 
        i += 1
    end

    # update struct state
    m.ndays_forecast -= m.ndays_forecast_used
    @assert m.ndays_forecast >= 0   # TODO: remove validation
    
    m.prev_t = t    

    return
end

function estimate_u0(m::TwoStateIfmHandler, t::ProbTime)
    if isnothing(m.prev_t)
        _initial_data_obs_update(m, t)
    else
        _data_obs_update(m, t)
    end

    # do prediction from start of obs up until today
    (S0, G0) = (Float32(0), Float32(0)) 

    # create interpolation input functions
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.data_obs.timepoints, m.data_obs.P, itp_method)
    itp_T = interpolate(m.data_obs.timepoints, m.data_obs.T, itp_method)
    itp_Lday = interpolate(m.data_obs.timepoints, m.data_obs.Lday, itp_method)

    (__, OED_sol) = m.predictor.predict(S0, G0, itp_Lday, itp_P, itp_T, m.data_obs.timepoints)

    # extract states
    est_S0 = Float64(last(OED_sol[1, :]))
    est_G0 = Float64(last(OED_sol[2, :]))

    return (est_S0, est_G0)
end

function predict(m::TwoStateIfmHandler, u0::Vector{Float64}, t::ProbTime)
    update_prediction_data(m, m.updater, t)

    # create interpolation input functions
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.data_pred.timepoints, m.data_pred.P, itp_method)
    itp_T = interpolate(m.data_pred.timepoints, m.data_pred.T, itp_method)
    itp_Lday = interpolate(m.data_pred.timepoints, m.data_pred.Lday, itp_method)

    (S0, G0) = u0
    (Q, OED_sol) = m.predictor.predict(S0, G0, itp_Lday, itp_P, itp_T, m.data_pred.timepoints)

    Q = Float64.(Q)
    u = Float64.(OED_sol.u)

    Q .= Q .* m.m3s_per_mm

    return (Q, u)
end


struct SimpleIfmDataUpdater <: AbstractTwoStateIfmDataUpdater
end

function update_prediction_data(m::TwoStateIfmHandler, updater::SimpleIfmDataUpdater, t::ProbTime)
    if m.ndays_forecast > 0
        # use forecast for P and T if available
        ndays_before_estimate_u0_call = m.ndays_forecast + m.ndays_forecast_used
        startix = length(m.data_forecast) - ndays_before_estimate_u0_call + 1
        stopix = startix + m.ndays_forecast_used
        for (i, j) in enumerate(startix:stopix)
            m.data_pred.P[i] = m.data_forecast.P[j]
            m.data_pred.T[i] = m.data_forecast.T[j]
        end
        # always use hist to update Lday
        for i in 1:m.ndays_forecast_used
                start = getscenariotime(t) + Day(i-1)
                m.data_pred.Lday[i] = getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
        end
    end
    # use hist to update after forecast period
    for i in (m.ndays_forecast_used + 1):m.pred_ndays
        start = getscenariotime(t) + Day(i-1)
        m.data_pred.P[i] = getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA)
        m.data_pred.T[i] = getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
        m.data_pred.Lday[i] = getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
    end
    return
end

struct TwoStateBucketIfmPredictor{P} <: AbstractTwoStateIfmPredictor
    model_params::P
end

function predict(m::TwoStateBucketIfmPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    basic_bucket_incl_states([S0, G0, m.model_params...], itp_Lday, itp_P, itp_T, timepoints)
end

struct TwoStateBucketIfm{H} <: AbstractInflowModel
    id::Id
    handler::H

    function TwoStateBucketIfm(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        predictor = TwoStateBucketIfmPredictor(model_params)
        handler = TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)        
        return new{typeof(handler)}(id, handler)
    end
end

estimate_u0(m::TwoStateBucketIfm, t::ProbTime) = estimate_u0(m.handler, t)
predict(m::TwoStateBucketIfm, u0::Vector{Float64}, t::ProbTime) = predict(m.handler, u0, t)

function includeTwoStateBucketIfm!(toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    common_includeTwoStateIfm!(TwoStateBucketIfm, toplevel, lowlevel, elkey, value)
end

struct TwoStateNeuralODEIfmPredictor{P, NN} <: AbstractTwoStateIfmPredictor
    nn::NN
    nn_params::P
    mean_S::Float32
    mean_G::Float32
    mean_P::Float32
    mean_T::Float32
    std_S::Float32
    std_G::Float32
    std_P::Float32
    std_T::Float32
    function TwoStateNeuralODEIfmPredictor(model_params)
        (nn, __) = initialize_NN_model()
        (nn_params, moments) = model_params
        return new{typeof(nn_params), typeof(nn)}(nn_params, nn, moments...)
    end
end

function predict(m::TwoStateNeuralODEIfmPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    norm_S = (S) -> (S.-m.mean_S)./m.std_S
    norm_G = (G) -> (G.-m.mean_G)./m.std_G
    norm_P = (P) -> (P.-m.mean_P)./m.std_P
    norm_T = (T) -> (T.-m.mean_T)./m.std_T
    NeuralODE_M100(m.nn_params, norm_S, norm_G, norm_P, norm_T, 
        S0, G0, itp_Lday, itp_P, itp_T, timepoints, m.nn)
end

struct TwoStateNeuralODEIfm{H} <: AbstractInflowModel
    id::Id
    handler::H

    function TwoStateNeuralODEIfm(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                    ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        predictor = TwoStateNeuralODEIfmPredictor(model_params)
        handler = TwoStateIfmHandler(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)        
        return new{typeof(handler)}(id, handler)
    end
end

estimate_u0(m::TwoStateNeuralODEIfm, t::ProbTime) = estimate_u0(m.handler, t)
predict(m::TwoStateNeuralODEIfm, u0::Vector{Float64}, t::ProbTime) = predict(m.handler, u0, t)

function common_includeTwoStateIfm!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    checkkey(toplevel, elkey)

    OBSERVED_PERCIPITATION = "ObservedPercipitation"
    OBSERVED_TEMPERATURE = "ObservedTemperature"
    FORECASTED_PERCIPITATION = "ForecastedPercipitation"
    FORECASTED_TEMPERATURE = "ForecastedTemperature"
    
    model_params = getdictvalue(value, "ModelParams", (Any, ), elkey)
    if model_params isa String
        model_params = JLD2.load_object(model_params)
    end

    hist_P = getdictvalue(value, "HistoricalPercipitation",   TIMEVECTORPARSETYPES, elkey)
    hist_T = getdictvalue(value, "HistoricalTemperature", TIMEVECTORPARSETYPES, elkey)
    hist_Lday = getdictvalue(value, "HistoricalDaylight", TIMEVECTORPARSETYPES, elkey)

    ndays_pred = getdictvalue(value, "NDaysPred", Int, elkey)
    basin_area = getdictvalue(value, "BasinArea", Float64, elkey)

    if haskey(value, OBSERVED_PERCIPITATION)
        obs_P = getdictvalue(value, OBSERVED_PERCIPITATION,   Vector{Real}, elkey)
    else
        obs_P = nothing
    end
    if haskey(value, OBSERVED_TEMPERATURE)
        obs_T = getdictvalue(value, OBSERVED_TEMPERATURE,   Vector{Real}, elkey)
    else
        obs_T = nothing
    end
    
    if haskey(value, FORECASTED_PERCIPITATION)
        forecast_P = getdictvalue(value, FORECASTED_PERCIPITATION,   Vector{Real}, elkey)
    else
        forecast_P = nothing
    end
    if haskey(value, FORECASTED_TEMPERATURE)
        forecast_T = getdictvalue(value, FORECASTED_TEMPERATURE,   Vector{Real}, elkey)
    else
        forecast_T = nothing
    end

    isnothing(obs_P) && (!isnothing(obs_T)) && error("Missing $OBSERVED_PERCIPITATION for $elkey")
    isnothing(obs_T) && (!isnothing(obs_P)) && error("Missing $OBSERVED_TEMPERATURE for $elkey")

    isnothing(obs_P) && (!isnothing(obs_T)) && error("Missing $FORECASTED_PERCIPITATION for $elkey")
    isnothing(obs_T) && (!isnothing(obs_P)) && error("Missing $FORECASTED_TEMPERATURE for $elkey")

    deps = Id[]
    all_ok = true

    (id, percipitation, ok) = getdicttimevectorvalue(lowlevel, percipitation)    
    all_ok = all_ok && ok
    _update_deps(deps, id, ok)
    
    (id, temperature, ok) = getdicttimevectorvalue(lowlevel, temperature)  
    all_ok = all_ok && ok
    _update_deps(deps, id, ok)

    (id, daylight, ok) = getdicttimevectorvalue(lowlevel, daylight)  
    all_ok = all_ok && ok
    _update_deps(deps, id, ok)

    if all_ok == false
        return (false, deps)
    end

    if !(isnothing(obs_P) || isnothing(obs_T))
        dummy_Lday = zeros(eltype(obs_P), length(obs_P))
        data_obs = TwoStateIfmData(obs_P, obs_T, dummy_Lday)
        ndays_obs == length(data_obs)
    else
        data_obs = nothing
        ndays_obs = 365
    end


    if !(isnothing(forecast_P) || isnothing(forecast_T))
        dummy_Lday = zeros(eltype(forecast_P), length(forecast_P))
        data_forecast = TwoStateIfmData(forecast_P, forecast_T, dummy_Lday)
        ndays_forecast = length(data_forecast)
    else
        data_forecast = nothing
        ndays_forecast = 0
    end

    id = getobjkey(elkey)
    toplevel[id] = Constructor(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)

    return (true, deps)
end


# --- Functions used in run_serial in connection with inflow models ---

"""
Create inflow models and store some of them locally according to db.dist_ifm
"""
function create_ifm()
    db = get_local_db()
    elements = get_ifm_elements(db)
    modelobjects = getmodelobjects(elements)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            id = Id(ABSTRACT_INFLOW_MODEL, inflow_name)
            db.ifm[inflow_name] = modelobjects[id]
        end
    end
end

"""
Sequentially solve inflow models stored locally. 
Each inflow model is solved for each scenario.
"""
function solve_ifm(t)
    db = get_local_db()
    normfactors = get_ifm_normfactors(db)
    scenarios = get_scenarios(db.scenmod_ppp)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            if !haskey(db.ifm_output, inflow_name)
                db.ifm_output[inflow_name] = Dict()
            end
            inflow_model = db.ifm[inflow_name]
            normalize_factor = normfactors[inflow_name]
            u0 = estimate_u0(inflow_model, t)
            for (scenix, scen) in enumerate(scenarios)
                scentime = get_scentphasein(t, scen, db.input)
                (Q, u) = predict(inflow_model, u0, scentime)
                Q .= Q .* normalize_factor
                start = getscenariotime(scentime)
                if !haskey(db.ifm_output[inflow_name], scenix)
                    ix = [start + Day(i-1) for i in 1:length(Q)]
                    db.ifm_output[inflow_name][scenix] = (ix, Q, u)
                else
                    (stored_ix, stored_Q, stored_u) = db.ifm_output[inflow_name][scenix]
                    for in in 1:length(Q)
                        stored_ix[i] = start + Day(i-1)
                        stored_Q[i] = Q[i]
                        stored_u[i] = u[i]
                    end
                end
            end
        end
    end
end

"""
Ensure that all cores have the latest output for all inflow models for all scenarios.
A core holding output from an inflow model, copies the output to the local db on all other cores.
"""
function synchronize_ifm_output()
    db = get_local_db()
    cores = get_cores(db)
    @sync for (inflow_name, core) in db.dist_ifm
        if core == db.core
            data = db.ifm_output[inflow_name]
            for other_core in cores
                if other_core != db.core
                    @spawnat other_core set_ifm_output(data, inflow_name)
                end
            end
        end
    end
end

function set_ifm_output(data, inflow_name)
    db = get_local_db()
    db.ifm_output[inflow_name] = data
    return
end

"""
Some inflow profiles are derived from (they are weighted averages of) 
other inflow profiles. This function utdates the derived profiles data 
so as to reflect the latest underlying output.
"""
function update_ifm_derived()
    db = get_local_db()
    ifm_weights = get_ifm_weights(db)
    for derived_name in keys(db.ifm_derived)
        weights = ifm_weights[derived_name]
        for (scenix, (derived_ix, derived_vals)) in db.ifm_derived[derived_name]
            do_ix = true
            fill!(derived_vals, 0.0)
            for (inflow_name, weight) in weights
                (ix, vals) = db.ifm_output[inflow_name][scenix]
                if do_ix
                    derived_ix .= ix
                    do_ix = false
                end
                derived_vals .= derived_vals .+ weight .* vals
            end
        end
    end
end
