# TODO: Illustrate how data_obs, data_pred and data_forecast are updated over time

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

getndays(x::TwoStateIfmData) = length(x.P)

const ONEDAY_MS_TIMEDELTA = TuLiPa.MsTimeDelta(Day(1))

mutable struct TwoStateIfmHandler{P <: AbstractTwoStateIfmPredictor, 
                                U <: AbstractTwoStateIfmDataUpdater, 
                                T1 <: TuLiPa.TimeVector, T2 <: TuLiPa.TimeVector, T3 <: TuLiPa.TimeVector}
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

    prev_t::Union{TuLiPa.ProbTime, Nothing}
    ndays_forecast_used::Int

    function TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
            ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        @assert ndays_forecast >= 0
        @assert ndays_pred >= ndays_forecast
        @assert ndays_obs > 0
        isnothing(data_obs) || @assert ndays_obs == getndays(data_obs)
        isnothing(data_forecast) || @assert ndays_forecast == getndays(data_forecast)
        (ndays_forecast > 0) && @assert !isnothing(data_forecast)
        data_pred = TwoStateIfmData(ndays_pred)
        m3s_per_mm = 1/((1000^3)/(basin_area*(10^6))*86400)
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

function _initial_data_obs_update(m::TwoStateIfmHandler, t::TuLiPa.ProbTime)
    if isnothing(m.data_obs)
        m.data_obs = TwoStateIfmData(m.ndays_obs)
        for i in 1:m.ndays_obs
            ndays_back = m.ndays_obs + 1 - i
            start = TuLiPa.getscenariotime(t) - Day(ndays_back)
            m.data_obs.P[i] = TuLiPa.getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.T[i] = TuLiPa.getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.Lday[i] = TuLiPa.getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
        end
    end
    return
end

function _data_obs_update(m::TwoStateIfmHandler, t::TuLiPa.ProbTime)
    # calc ndays_update
    diff_ndays_scentime = Day(TuLiPa.getscenariotime(t) - TuLiPa.getscenariotime(m.prev_t)).value
    ndays_update = min(m.ndays_obs, diff_ndays_scentime)
    
    # shift backwards to make room for ndays_update new values
    i = 1 + ndays_update
    j = 1
    for __ in 1:(m.ndays_obs - ndays_update)
        m.data_obs.P[j] = m.data_obs.P[i]
        m.data_obs.T[j] = m.data_obs.T[i]
        m.data_obs.Lday[j] = m.data_obs.Lday[i]
        i += 1
        j += 1
    end

    # use forecast for P and T if available
    # TODO: Test this case
    m.ndays_forecast_used = 0
    if m.ndays_forecast > 0
        startix = getndays(m.data_forecast) - m.ndays_forecast + 1
        stopix = getndays(m.data_forecast) - m.ndays_forecast + ndays_update
        stopix = min(stopix, getndays(m.data_forecast))
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
        start = TuLiPa.getscenariotime(t) - Day(ndays_remaining)
        i = m.ndays_obs - ndays_remaining + 1
        for __ in 1:ndays_remaining
            m.data_obs.P[i] = TuLiPa.getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.T[i] = TuLiPa.getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
            start += Day(1) 
            i += 1
        end
    end

    # always use hist to update Lday
    start = TuLiPa.getscenariotime(t) - Day(ndays_update + 1)
    i = m.ndays_obs - ndays_update + 1
    for __ in 1:ndays_update
        m.data_obs.Lday[i] = TuLiPa.getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
        start += Day(1) 
        i += 1
    end

    # update struct state
    m.ndays_forecast -= m.ndays_forecast_used
    @assert m.ndays_forecast >= 0   # TODO: remove validation
    return
end

function estimate_u0(m::TwoStateIfmHandler, t::TuLiPa.ProbTime)
    if isnothing(m.prev_t)
        _initial_data_obs_update(m, t)
    else
        _data_obs_update(m, t)
    end
    m.prev_t = t

    # do prediction from start of obs up until today
    (S0, G0) = (Float32(0), Float32(0))

    # create interpolation input functions
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.data_obs.timepoints, m.data_obs.P, itp_method)
    itp_T = interpolate(m.data_obs.timepoints, m.data_obs.T, itp_method)
    itp_Lday = interpolate(m.data_obs.timepoints, m.data_obs.Lday, itp_method)

    (__, OED_sol) = predict(m.predictor, S0, G0, itp_Lday, itp_P, itp_T, m.data_obs.timepoints)

    # extract states
    est_S0 = Float64(last(OED_sol[1, :]))
    est_G0 = Float64(last(OED_sol[2, :]))

    return [est_S0, est_G0]
end

function predict(m::TwoStateIfmHandler, u0::Vector{Float64}, t::TuLiPa.ProbTime)
    update_prediction_data(m, m.updater, t)

    # create interpolation input functions
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.data_pred.timepoints, m.data_pred.P, itp_method)
    itp_T = interpolate(m.data_pred.timepoints, m.data_pred.T, itp_method)
    itp_Lday = interpolate(m.data_pred.timepoints, m.data_pred.Lday, itp_method)

    (S0, G0) = u0
    (Q, OED_sol) = predict(m.predictor, S0, G0, itp_Lday, itp_P, itp_T, m.data_pred.timepoints)

    Q = Float64.(Q)

    Q .= Q .* m.m3s_per_mm

    return Q
end


struct SimpleIfmDataUpdater <: AbstractTwoStateIfmDataUpdater
end

function update_prediction_data(m::TwoStateIfmHandler, updater::SimpleIfmDataUpdater, t::TuLiPa.ProbTime)
    if m.ndays_forecast > 0
        # use forecast for P and T if available
        ndays_before_estimate_u0_call = m.ndays_forecast + m.ndays_forecast_used
        startix = getndays(m.data_forecast) - ndays_before_estimate_u0_call + 1
        stopix = startix + m.ndays_forecast_used
        for (i, j) in enumerate(startix:stopix)
            m.data_pred.P[i] = m.data_forecast.P[j]
            m.data_pred.T[i] = m.data_forecast.T[j]
        end
        # always use hist to update Lday
        for i in 1:m.ndays_forecast_used
                start = TuLiPa.getscenariotime(t) + Day(i-1)
                m.data_pred.Lday[i] = TuLiPa.getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
        end
    end
    # use hist to update after forecast period
    for i in (m.ndays_forecast_used + 1):m.ndays_pred
        start = TuLiPa.getscenariotime(t) + Day(i-1)
        m.data_pred.P[i] = TuLiPa.getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA)
        m.data_pred.T[i] = TuLiPa.getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
        m.data_pred.Lday[i] = TuLiPa.getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
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
    id::TuLiPa.Id
    handler::H

    function TwoStateBucketIfm(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        predictor = TwoStateBucketIfmPredictor(model_params)
        handler = TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)        
        return new{typeof(handler)}(id, handler)
    end
end

TuLiPa.getid(m::TwoStateBucketIfm) = m.id
TuLiPa.assemble!(m::TwoStateBucketIfm) = true
numstates(::TwoStateBucketIfm) = 2
estimate_u0(m::TwoStateBucketIfm, t::TuLiPa.ProbTime) = estimate_u0(m.handler, t)
predict(m::TwoStateBucketIfm, u0::Vector{Float64}, t::TuLiPa.ProbTime) = predict(m.handler, u0, t)

function includeTwoStateBucketIfm!(toplevel::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
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
        (nn_params, d) = model_params
        return new{typeof(nn_params), typeof(nn)}(nn_params, nn, 
            d["mean_S"], d["mean_G"], d["mean_P"], d["mean_T"], 
            d["std_S"], d["std_G"], d["std_P"], d["std_T"])
    end
end

function predict(m::TwoStateNeuralODEIfmPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    norm_S = (S) -> (S.-m.mean_S)./m.std_S
    norm_G = (G) -> (G.-m.mean_G)./m.std_G
    norm_P = (P) -> (P.-m.mean_P)./m.std_P
    norm_T = (T) -> (T.-m.mean_T)./m.std_T
    NeuralODE_M100(m.nn_params, norm_S, norm_G, norm_P, norm_T, itp_Lday, itp_P, itp_T, timepoints, m.nn; S_init = [S0, G0])
end

struct TwoStateNeuralODEIfm{H} <: AbstractInflowModel
    id::TuLiPa.Id
    handler::H

    function TwoStateNeuralODEIfm(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                    ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        predictor = TwoStateNeuralODEIfmPredictor(model_params)
        handler = TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                    ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)  
        return new{typeof(handler)}(id, handler)
    end
end

TuLiPa.getid(m::TwoStateNeuralODEIfm) = m.id
TuLiPa.assemble!(m::TwoStateNeuralODEIfm) = true
numstates(::TwoStateNeuralODEIfm) = 2
estimate_u0(m::TwoStateNeuralODEIfm, t::TuLiPa.ProbTime) = estimate_u0(m.handler, t)
predict(m::TwoStateNeuralODEIfm, u0::Vector{Float64}, t::TuLiPa.ProbTime) = predict(m.handler, u0, t)

function includeTwoStateNeuralODEIfm!(toplevel::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
    common_includeTwoStateIfm!(TwoStateNeuralODEIfm, toplevel, lowlevel, elkey, value)
end

function common_includeTwoStateIfm!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
    TuLiPa.checkkey(toplevel, elkey)

    OBSERVED_PERCIPITATION = "ObservedPercipitation"
    OBSERVED_TEMPERATURE = "ObservedTemperature"
    FORECASTED_PERCIPITATION = "ForecastedPercipitation"
    FORECASTED_TEMPERATURE = "ForecastedTemperature"

    model_params = TuLiPa.getdictvalue(value, "ModelParams", String, elkey)

    moments = nothing
    if haskey(value, "Moments")
        moments = TuLiPa.getdictvalue(value, "Moments", String, elkey)
    end

    hist_P = TuLiPa.getdictvalue(value, "HistoricalPercipitation",   TuLiPa.TIMEVECTORPARSETYPES, elkey)
    hist_T = TuLiPa.getdictvalue(value, "HistoricalTemperature", TuLiPa.TIMEVECTORPARSETYPES, elkey)
    hist_Lday = TuLiPa.getdictvalue(value, "HistoricalDaylight", TuLiPa.TIMEVECTORPARSETYPES, elkey)

    ndays_pred = TuLiPa.getdictvalue(value, "NDaysPred", Real, elkey)
    try 
        ndays_pred = Int(ndays_pred)
        @assert ndays_pred >= 0
    catch e
        error("Value for key NDaysPred must be positive integer for $elkey")
    end

    basin_area = TuLiPa.getdictvalue(value, "BasinArea", Float64, elkey)

    if haskey(value, OBSERVED_PERCIPITATION)
        obs_P = TuLiPa.getdictvalue(value, OBSERVED_PERCIPITATION,   Vector{Real}, elkey)
    else
        obs_P = nothing
    end
    if haskey(value, OBSERVED_TEMPERATURE)
        obs_T = TuLiPa.getdictvalue(value, OBSERVED_TEMPERATURE,   Vector{Real}, elkey)
    else
        obs_T = nothing
    end
    
    if haskey(value, FORECASTED_PERCIPITATION)
        forecast_P = TuLiPa.getdictvalue(value, FORECASTED_PERCIPITATION,   Vector{Real}, elkey)
    else
        forecast_P = nothing
    end
    if haskey(value, FORECASTED_TEMPERATURE)
        forecast_T = TuLiPa.getdictvalue(value, FORECASTED_TEMPERATURE,   Vector{Real}, elkey)
    else
        forecast_T = nothing
    end

    isnothing(obs_P) && (!isnothing(obs_T)) && error("Missing $OBSERVED_PERCIPITATION for $elkey")
    isnothing(obs_T) && (!isnothing(obs_P)) && error("Missing $OBSERVED_TEMPERATURE for $elkey")

    isnothing(obs_P) && (!isnothing(obs_T)) && error("Missing $FORECASTED_PERCIPITATION for $elkey")
    isnothing(obs_T) && (!isnothing(obs_P)) && error("Missing $FORECASTED_TEMPERATURE for $elkey")

    deps = TuLiPa.Id[]
    # errs = String[]

    all_ok = true

    (id, hist_P, ok) = TuLiPa.getdicttimevectorvalue(lowlevel, hist_P)    
    all_ok = all_ok && ok
    TuLiPa._update_deps(deps, id, ok)
    
    (id, hist_T, ok) = TuLiPa.getdicttimevectorvalue(lowlevel, hist_T)  
    all_ok = all_ok && ok
    TuLiPa._update_deps(deps, id, ok)

    (id, hist_Lday, ok) = TuLiPa.getdicttimevectorvalue(lowlevel, hist_Lday)  
    all_ok = all_ok && ok
    TuLiPa._update_deps(deps, id, ok)

    if all_ok == false
        return (false, deps)
        # return (false, deps, errs)
    end

    if !(isnothing(obs_P) || isnothing(obs_T))
        dummy_Lday = zeros(eltype(obs_P), length(obs_P))
        data_obs = TwoStateIfmData(obs_P, obs_T, dummy_Lday)
        ndays_obs == getndays(data_obs)
    else
        data_obs = nothing
        ndays_obs = 365
    end

    if !(isnothing(forecast_P) || isnothing(forecast_T))
        dummy_Lday = zeros(eltype(forecast_P), length(forecast_P))
        data_forecast = TwoStateIfmData(forecast_P, forecast_T, dummy_Lday)
        ndays_forecast = getndays(data_forecast)
    else
        data_forecast = nothing
        ndays_forecast = 0
    end

    # TODO: Maybe make this user input in future?
    updater = SimpleIfmDataUpdater()

    model_params = JLD2.load_object(model_params)

    if !isnothing(moments)
        moments = JLD2.load_object(moments)
        model_params = (model_params, moments)
    end

    id = TuLiPa.getobjkey(elkey)
    toplevel[id] = Constructor(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)

    return (true, deps)
    # return (true, deps, errs)
end


# --- Functions used in run_serial in connection with inflow models ---

const IFM_DB_STATE_KEY = "ifm_step_u0"

"""
Create inflow models and store some of them locally according to db.dist_ifm
"""
function create_ifm()
    db = get_local_db()
    @assert !haskey(db.div, IFM_DB_STATE_KEY)
    db.div[IFM_DB_STATE_KEY] = Dict{String, Tuple{Int, Vector{Float64}}}()
    elements = get_ifm_elements(db)
    t0 = time()
    modelobjects = TuLiPa.getmodelobjects(elements)
    t1 = time()
    sec = t1 - t0
    println("getmodelobjects(ifm_elements) took $sec")
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            id = TuLiPa.Id(ABSTRACT_INFLOW_MODEL, inflow_name)
            db.ifm[inflow_name] = modelobjects[id]
        end
    end
end

function save_ifm_u0(db, inflow_name, stepnr, u0)
    if !haskey(db.div[IFM_DB_STATE_KEY], inflow_name)
        db.div[IFM_DB_STATE_KEY][inflow_name] = (stepnr, u0)
    else
        (stored_stepnr, __) = db.div[IFM_DB_STATE_KEY][inflow_name]
        if stored_stepnr != stepnr
            db.div[IFM_DB_STATE_KEY][inflow_name] = (stepnr, u0)
        end
    end
end

"""
Sequentially solve inflow models stored locally. 
Each inflow model is solved for each scenario.
"""
function solve_ifm(t, stepnr)
    db = get_local_db()
    normfactors = get_ifm_normfactors(db)
    scenarios = get_scenarios(db.scenmod_sim)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            if !haskey(db.ifm_output, inflow_name)
                db.ifm_output[inflow_name] = Dict()
            end
            inflow_model = db.ifm[inflow_name]
            normalize_factor = normfactors[inflow_name]
            u0 = estimate_u0(inflow_model, t)
            save_ifm_u0(db, inflow_name, stepnr, u0)
            for (scenix, scen) in enumerate(scenarios)
                scentime = get_scentphasein(t, scen, db.input)
                Q = predict(inflow_model, u0, scentime)
                Q .= Q .* normalize_factor
                start = TuLiPa.getdatatime(scentime)
                if !haskey(db.ifm_output[inflow_name], scenix)
                    ix = [start + Day(i-1) for i in 1:length(Q)]
                    db.ifm_output[inflow_name][scenix] = (ix, Q)
                else
                    (stored_ix, stored_Q) = db.ifm_output[inflow_name][scenix]
                    if length(stored_Q) != length(Q)
                        # ModeledInflowParam added empty vectors which it points to
                        # Therefore, it is important not to replace the stored vectors, 
                        # we therefore push values to it
                        @assert length(stored_ix) == 0
                        @assert length(stored_Q) == 0
                        for i in eachindex(Q)
                            push!(stored_ix, start + Day(i-1))
                            push!(stored_Q, Q[i])
                        end
                    else
                        for i in eachindex(Q)
                            stored_ix[i] = start + Day(i-1)
                            stored_Q[i] = Q[i]
                        end
                    end
                end
            end

            # Ifm for stoch scenarios TODO: Quickfix, improve?
            stochscenarios = get_scenarios(db.scenmod_stoch)
            if length(stochscenarios) != length(scenarios)
                for (j, stochscen) in enumerate(stochscenarios)
                    stochscenix = length(scenarios) + j
                    parentscenix = stochscen.parentscenario
                    (ix, Q) = db.ifm_output[inflow_name][parentscenix]
                    if !haskey(db.ifm_output[inflow_name], stochscenix)
                        db.ifm_output[inflow_name][stochscenix] = (ix, Q)
                    else
                        (stored_ix, stored_Q) = db.ifm_output[inflow_name][stochscenix]
                        if length(stored_Q) != length(Q)
                            @assert length(stored_ix) == 0
                            @assert length(stored_Q) == 0
                            for i in eachindex(Q)
                                push!(stored_ix, ix[i])
                                push!(stored_Q, Q[i])
                            end
                        else
                            for i in eachindex(Q)
                                stored_ix[i] = ix[i]
                                stored_Q[i] = Q[i]
                            end
                        end
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
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            data = db.ifm_output[inflow_name]
            @sync for other_core in cores
                if other_core != db.core
                    @spawnat other_core set_ifm_output(data, inflow_name)
                end
            end
        end
    end
end

"""
Write content of data to stored vectors in local ifm_output slot
Make sure to reuse existing vectors, since other objects hold
pointers to them, and expect them to be hold in sync by this code
"""
function set_ifm_output(data, inflow_name)
    db = get_local_db()
    for scenix in keys(data)
        if !haskey(db.ifm_output, inflow_name)
            db.ifm_output[inflow_name] = Dict()
        end
        if !haskey(db.ifm_output[inflow_name], scenix)
            db.ifm_output[inflow_name][scenix] = (DateTime[], Float64[])
        end
        (incoming_ix, incoming_vals) = data[scenix]
        (stored_ix, stored_vals) = db.ifm_output[inflow_name][scenix]
        @assert length(stored_ix) == length(stored_vals)
        @assert length(incoming_ix) == length(incoming_vals)
        if length(incoming_ix) != length(stored_ix)
            @assert length(stored_ix) == 0
            for i in eachindex(incoming_ix)
                push!(stored_ix, incoming_ix[i])
                push!(stored_vals, incoming_vals[i])
            end
        else
            for i in eachindex(incoming_ix)
                stored_ix[i] = incoming_ix[i]
                stored_vals[i] = incoming_vals[i]
            end
        end
    end
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
                if length(vals) != length(derived_vals)
                    @assert length(derived_vals) == 0
                    for (i, v) in enumerate(vals)
                        push!(derived_vals, v * weight)
                        push!(derived_ix, ix[i])
                        do_ix = false
                    end
                else
                    if do_ix
                        derived_ix .= ix
                        do_ix = false
                    end
                    derived_vals .= derived_vals .+ weight .* vals
                end
            end
        end
    end
end


"""
ModeledInflowParam is actually a stateful PrognosisSeriesParam under-the-hood,
but where the prognosis part is created and managed by JulES and not 
given as user input
"""
function includeModeledInflowParam!(::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
    TuLiPa.checkkey(lowlevel, elkey)
    
    deps = TuLiPa.Id[]
    all_ok = true

    # Not part of user input 
    # This info is added by JulES with add_scenix_to_ModeledInflowParam
    # See e.g. prob_stoch.get_elements_with_horizons
    scenix = TuLiPa.getdictvalue(value, "ScenarioIndex", Int, elkey)

    level = TuLiPa.getdictvalue(value, "Level", TuLiPa.TIMEVECTORPARSETYPES, elkey)
    (id, level, ok) = TuLiPa.getdicttimevectorvalue(lowlevel, level)
    all_ok = all_ok && ok
    TuLiPa._update_deps(deps, id, ok)

    hist_profile = TuLiPa.getdictvalue(value, "HistoricalProfile", TuLiPa.TIMEVECTORPARSETYPES, elkey)
    (id, hist_profile, ok) = TuLiPa.getdicttimevectorvalue(lowlevel, hist_profile)
    all_ok = all_ok && ok
    TuLiPa._update_deps(deps, id, ok)

    if all_ok == false
        return (false, deps)
    end

    # Creates an InfiniteTimeVector that refers to vectors stored in local db, 
    # which will be updated by JulES each step after running inflow models.
    # This way, model objects holding reference to such Param, will use updated 
    # prognosis from db when called upon by update!(prob, t)
    db = get_local_db()
    replacemap = get_ifm_replacemap(db.input)
    haskey(replacemap, elkey.instancename) || error("Instance name not found in replacemap for $elkey")
    inflow_name = replacemap[elkey.instancename]
    ifm_weights = get_ifm_weights(db)
    ifm_names = get_ifm_names(db.input)
    if haskey(ifm_weights, inflow_name)
        d = db.ifm_derived
    else
        @assert (inflow_name in ifm_names)
        d = db.ifm_output
    end
    if !haskey(d, inflow_name)
        d[inflow_name] = Dict()
    end
    if !haskey(d[inflow_name], scenix)
        d[inflow_name][scenix] = (DateTime[], Float64[])
    end
    (ix, vals) = d[inflow_name][scenix]
    prognosis_profile = TuLiPa.InfiniteTimeVector(ix, vals)    

    param = TuLiPa.PrognosisSeriesParam(level, hist_profile, prognosis_profile)
    param = TuLiPa.StatefulParam(param)

    lowlevel[TuLiPa.getobjkey(elkey)] = param

    return (true, deps)
end

"""
MixedInflowParam...
"""
# TODO: Complete this
# struct MixedInflowParam{P1, P2} <: Param
#     directparam::P1
#     ifmparam::P2
#     ndays::Int
# end
# function includeMixedInflowParam!(::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
# end
# TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(TuLiPa.PARAM_CONCEPT, "MixedInflowParam")] = includeMixedInflowParam!

"""
Used by JulES in appropriate places to embed scenix info 
into data elements of type ModeledInflowParam or MixedInflowParam
"""
function add_scenix_to_InflowParam(elements, scenix)
    for e in elements
        if e.typename == "ModeledInflowParam" || e.typename == "MixedInflowParam"
            e.value["ScenarioIndex"] = scenix
        end
    end
end

"""
Return copy of elements with replacement of PrognosisSeriesParam 
in ifm_replacemap in accordance with value of iprogtype
"""
function copy_elements_iprogtype(elements, iprogtype, ifm_replacemap)
    if iprogtype == "ifm"
        elements1 = TuLiPa.DataElement[]
        for e in elements
            if e.typename == "PrognosisSeriesParam" && haskey(ifm_replacemap, e.instancename)
                new_e = TuLiPa.DataElement(e.conceptname, "ModeledInflowParam", e.instancename,
                    Dict("Level" => e.value["Level"], "HistoricalProfile" => e.value["Profile"]))
                push!(elements1, new_e)
            else
                push!(elements1, e)
            end
        end

    elseif startswith(iprogtype, "mix")
        # TODO: Validate ndays > 0 in constructor of DefaultJulESInput
        ndays = parse(Int, iprogtype[4:end])
        elements1 = TuLiPa.DataElement[]
        for e in elements
            if e.typename == "PrognosisSeriesParam" && haskey(ifm_replacemap, e.instancename)
                new_value = copy(e.value::Dict)
                new_value["ndays"] = ndays
                new_e = TuLiPa.DataElement(e.conceptname, "MixedInflowParam", e.instancename, new_value)
                push!(elements1, new_e)
            else
                push!(elements1, e)
            end
        end

    else
        @assert iprogtype == "direct"
        elements1 = copy(elements)
    end

    return elements1
end
