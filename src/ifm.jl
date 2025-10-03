"""
Implementation of inflow models (ifm) TwoStateBucketIfm and TwoStateNeuralODEIfm
and integration with JulES
"""

# TODO: Add return nothing to many functions
# TODO: Illustrate how data_obs, data_pred and data_forecast are updated over time
# TODO: add data element for observations for an ifm (and extend interface with set_observations)
# TODO: add data element for forecast for an ifm (and extend interface with set_forecast)
# TODO: Use Float32 instead of Float64 in historical time vectors?
# TODO: Replace all calls to Day(n) with Millisecond(86400000*n) for better performance 

const ONEDAY_MS_TIMEDELTA = TuLiPa.MsTimeDelta(Day(1))

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

    scen_start::Any
    scen_stop::Any

    function TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
            ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)
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
            data_obs, data_forecast, nothing, 0, scen_start, scen_stop)
    end
end

get_numstates(m::TwoStateIfmHandler) = 2
get_basin_area_m2(m::TwoStateIfmHandler) = m.basin_area

function get_statename(::TwoStateIfmHandler, i::Int)
    if i == 1
        return "snow"
    elseif i == 2
        return "ground"
    else
        error("Only state ix 1 and 2 supported")
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
    (Q, __) = predict(m.predictor, S0, G0, itp_Lday, itp_P, itp_T, m.data_pred.timepoints)

    Q = Float64.(Q)

    Q .= Q .* m.m3s_per_mm

    return Q
end

# TODO: Add AutoCorrIfmDataUpdater that use value w(t)*x(t0) + (1-w(t-t0))*x(t), where w(0) = 1 and w -> 0 for larger inputs
struct SimpleIfmDataUpdater <: AbstractTwoStateIfmDataUpdater
end

function update_prediction_data(m::TwoStateIfmHandler, ::SimpleIfmDataUpdater, t::TuLiPa.ProbTime)
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
                                ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)
        predictor = TwoStateBucketIfmPredictor(model_params)
        handler = TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)        
        return new{typeof(handler)}(id, handler)
    end
end

TuLiPa.getid(m::TwoStateBucketIfm) = m.id
TuLiPa.assemble!(m::TwoStateBucketIfm) = true
estimate_u0(m::TwoStateBucketIfm, t::TuLiPa.ProbTime) = estimate_u0(m.handler, t)
predict(m::TwoStateBucketIfm, u0::Vector{Float64}, t::TuLiPa.ProbTime) = predict(m.handler, u0, t)
get_numstates(m::TwoStateBucketIfm) =  get_numstates(m.handler)
get_statename(m::TwoStateBucketIfm, i::Int) = get_statename(m.handler, i)
get_basin_area_m2(m::TwoStateBucketIfm) = get_basin_area_m2(m.handler)

function includeTwoStateBucketIfm!(toplevel::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
    common_includeTwoStateIfm!(TwoStateBucketIfm, toplevel, lowlevel, elkey, value)
end

struct TwoStateNeuralODEIfmPredictor{NN, P} <: AbstractTwoStateIfmPredictor
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
        return new{typeof(nn), typeof(nn_params)}(nn, nn_params, 
            d["mean_S0"], d["mean_S1"], d["mean_P"], d["mean_T"], 
            d["std_S0"], d["std_S1"], d["std_P"], d["std_T"])
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
                                    ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)
        predictor = TwoStateNeuralODEIfmPredictor(model_params)
        handler = TwoStateIfmHandler(predictor, updater, basin_area, hist_P, hist_T, hist_Lday, 
                                    ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)  
        return new{typeof(handler)}(id, handler)
    end
end

TuLiPa.getid(m::TwoStateNeuralODEIfm) = m.id
TuLiPa.assemble!(m::TwoStateNeuralODEIfm) = true
estimate_u0(m::TwoStateNeuralODEIfm, t::TuLiPa.ProbTime) = estimate_u0(m.handler, t)
predict(m::TwoStateNeuralODEIfm, u0::Vector{Float64}, t::TuLiPa.ProbTime) = predict(m.handler, u0, t)
get_numstates(m::TwoStateNeuralODEIfm) =  get_numstates(m.handler)
get_statename(m::TwoStateNeuralODEIfm, i::Int) = get_statename(m.handler, i)
get_basin_area_m2(m::TwoStateNeuralODEIfm) = get_basin_area_m2(m.handler)

function includeTwoStateNeuralODEIfm!(toplevel::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
    common_includeTwoStateIfm!(TwoStateNeuralODEIfm, toplevel, lowlevel, elkey, value)
end

function common_includeTwoStateIfm!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::TuLiPa.ElementKey, value::Dict)
    TuLiPa.checkkey(toplevel, elkey)

    model_params = TuLiPa.getdictvalue(value, "ModelParams", String, elkey)

    moments = nothing
    if haskey(value, "Moments")
        moments = TuLiPa.getdictvalue(value, "Moments", String, elkey)
    end

    hist_P = TuLiPa.getdictvalue(value, "HistoricalPercipitation", TuLiPa.TIMEVECTORPARSETYPES, elkey)
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

    deps = TuLiPa.Id[]

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
    end

    # TODO: Maybe make this user input in future?
    updater = SimpleIfmDataUpdater()

    model_params = JLD2.load_object(model_params)

    is_nn = !isnothing(moments)
    if is_nn
        # convert model_params, stored with simpler data structure for stability between versions,
        # into ComponentArray, which the NN-model needs
        # (the simpler data structure is Vector{Tuple{Vector{Float32}, Vector{Float32}}})
        _subarray(i) = ComponentArray(weight = model_params[i][1], bias = model_params[i][2])
        _tuple(i) = (Symbol("layer_", i), _subarray(i))
        model_params = ComponentArray(NamedTuple(_tuple(i) for i in eachindex(model_params)))
        # add moments, which is needed to normalize state inputs to the NN-model
        moments = JLD2.load_object(moments)
        model_params = (model_params, moments)
    end

    data_forecast = nothing
    data_obs = nothing
    ndays_obs = 365
    data_forecast = nothing
    ndays_forecast = 0


    
    periodkey = TuLiPa.Id(TuLiPa.TIMEPERIOD_CONCEPT, "ScenarioTimePeriod")
    period = lowlevel[periodkey]
    scen_start = period["Start"]
    scen_stop  = period["Stop"]


    id = TuLiPa.getobjkey(elkey)
    toplevel[id] = Constructor(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)

    return (true, deps)
end

# --- Functions used in run_serial in connection with inflow models ---

const IFM_DB_STATE_KEY = "ifm_step_u0"
const IFM_DB_FLOW_KEY = "ifm_step_Q"
const IFM_DB_ALLFLOW_KEY = "ifm_step_allQ"
const IFM_STATION_ID_KEY = "ifm_station_id"

"""
Create inflow models and store some of them locally according to db.dist_ifm
"""
function create_ifm()
    db = get_local_db()

    @assert !haskey(db.div, IFM_DB_STATE_KEY)
    db.div[IFM_DB_STATE_KEY] = Dict{String, Tuple{Int, Vector{Float64}}}()
    @assert !haskey(db.div, IFM_DB_FLOW_KEY)
    db.div[IFM_DB_FLOW_KEY] = Dict{String, Tuple{Int, Float64}}()
    @assert !haskey(db.div, IFM_DB_ALLFLOW_KEY)
    db.div[IFM_DB_ALLFLOW_KEY] = Dict{String, Vector{Float64}}()

    elements = get_ifm_elements(db)
    modelobjects = TuLiPa.getmodelobjects(elements)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            id = TuLiPa.Id(ABSTRACT_INFLOW_MODEL, inflow_name)
            db.ifm[inflow_name] = modelobjects[id]
        end
    end
end

function save_ifm_u0(div_db, inflow_name, stepnr, u0)
    if !haskey(div_db[IFM_DB_STATE_KEY], inflow_name)
        div_db[IFM_DB_STATE_KEY][inflow_name] = (stepnr, u0)
    else
        (stored_stepnr, __) = div_db[IFM_DB_STATE_KEY][inflow_name]
        if stored_stepnr != stepnr
            div_db[IFM_DB_STATE_KEY][inflow_name] = (stepnr, u0)
        end
    end
end

function save_ifm_Q(div_db, inflow_name, stepnr, Q)
    if !haskey(div_db[IFM_DB_FLOW_KEY], inflow_name)
        div_db[IFM_DB_FLOW_KEY][inflow_name] = (stepnr, Q)
    else
        (stored_stepnr, __) = div_db[IFM_DB_FLOW_KEY][inflow_name]
        if stored_stepnr != stepnr
            div_db[IFM_DB_FLOW_KEY][inflow_name] = (stepnr, Q)
        end
    end
end

function calculate_normalize_factor(ifm_model)
    start_with_buffer = ifm_model.handler.scen_start - Day(ifm_model.handler.ndays_obs) # Add ndays_obs days buffer
    days_with_buffer = Day(ifm_model.handler.scen_stop - start_with_buffer) |> Dates.value
    timepoints_with_buffer = (1:days_with_buffer)
    days = Dates.value(Day(ifm_model.handler.scen_stop - ifm_model.handler.scen_start))    
    timepoints_start = days_with_buffer - days 

    P = zeros(length(timepoints_with_buffer))
    T = zeros(length(timepoints_with_buffer))
    Lday = zeros(length(timepoints_with_buffer))
    for i in timepoints_with_buffer
        start = ifm_model.handler.scen_start + Day(i - 1)
        P[i] = TuLiPa.getweightedaverage(ifm_model.handler.hist_P, start, JulES.ONEDAY_MS_TIMEDELTA)
        T[i] = TuLiPa.getweightedaverage(ifm_model.handler.hist_T, start, JulES.ONEDAY_MS_TIMEDELTA)
        Lday[i] = TuLiPa.getweightedaverage(ifm_model.handler.hist_Lday, start, JulES.ONEDAY_MS_TIMEDELTA)
    end

    itp_method = JulES.SteffenMonotonicInterpolation()
    itp_P = JulES.interpolate(timepoints_with_buffer, P, itp_method)
    itp_T = JulES.interpolate(timepoints_with_buffer, T, itp_method)
    itp_Lday = JulES.interpolate(timepoints_with_buffer, Lday, itp_method)
    Q, _ = JulES.predict(ifm_model.handler.predictor, 0, 0, itp_Lday, itp_P, itp_T, timepoints_with_buffer)
    Q = Float64.(Q)[timepoints_start:end]
    Q .= Q .* ifm_model.handler.m3s_per_mm
    return 1 / mean(Q)
end

"""
Sequentially solve inflow models stored locally. 
Each inflow model is solved for each scenario.
"""
function solve_ifm(t, stepnr)
    db = get_local_db()

    # setup utility variables to save results below
    steplen_ms = Millisecond(get_steplength(db.input))
    ifmstep_ms = TuLiPa.getduration(ONEDAY_MS_TIMEDELTA)
    @assert steplen_ms >= ifmstep_ms
    nifmsteps = div(steplen_ms.value, ifmstep_ms.value)
    remainder_ms = steplen_ms - ifmstep_ms * nifmsteps
    steplen_f = float(steplen_ms.value)
    ifmstep_f = float(ifmstep_ms.value)
    remainder_f = float(remainder_ms.value)
    normfactors = get_ifm_normfactors(db)
    scenarios = get_scenarios(db.scenmod_sim)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            if !haskey(db.ifm_output, inflow_name)
                db.ifm_output[inflow_name] = Dict()
            end
            inflow_model = db.ifm[inflow_name]

            if haskey(normfactors, inflow_name) == false
                normalize_factor = calculate_normalize_factor(inflow_model) 
                normfactors[inflow_name] = normalize_factor
            else
                normalize_factor = normfactors[inflow_name]
            end

            u0 = estimate_u0(inflow_model, t)

            # save in familiar unit mm3 (u0 in mm converted to m in calculation)
            u0_mm3 = (u0 ./ 1000.0) .* get_basin_area_m2(inflow_model) ./ 1e6 
            save_ifm_u0(db.div, inflow_name, stepnr, u0_mm3)

            # predict mean Q for clearing period and store result
            # can be used to measure goodness of ifm model
            Q = predict(inflow_model, u0, t)
            @assert nifmsteps <= length(Q)
            mean_Q = 0.0
            for i in 1:nifmsteps
                mean_Q += Q[i] * ifmstep_f
            end
            if nifmsteps > 0
                mean_Q += Q[nifmsteps] * remainder_f
            end
            mean_Q /= steplen_f
            save_ifm_Q(db.div, inflow_name, stepnr, mean_Q)

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
Ensure that all cores have the latest output for all inflow models 
for all scenarios. A core holding output from an inflow model, 
copies the output to the local db on all other cores.
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
ModeledInflowParam is actually a stateful TuLiPa.PrognosisSeriesParam under-the-hood,
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
    inflow_name = value[IFM_STATION_ID_KEY]
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
    
    phaseindays = TuLiPa.getdictvalue(value, "phaseindays", Int, elkey)
    preddays = TuLiPa.getdictvalue(value, "preddays", Int, elkey)
    phaseinvalues = ones(preddays)
    if phaseindays == 0
        confidence = TuLiPa.ConstantTimeVector(1.0)
    else
        it = min(phaseindays, preddays)
        for i in 1:it
            phaseinvalues[i] = round((i-1)/phaseindays,digits=3)
        end
        confidence = TuLiPa.InfiniteTimeVector(ix, phaseinvalues)
    end

    param = TuLiPa.DynamicPrognosisSeriesParam(level, hist_profile, prognosis_profile, confidence)
    param = TuLiPa.StatefulParam(param)

    lowlevel[TuLiPa.getobjkey(elkey)] = param

    return (true, deps)
end

# TODO: Maybe support MixedInflowParam to use external forecast first year and ifm for rest?

"""
Used by JulES in appropriate places to embed scenix info 
into data elements of type ModeledInflowParam
"""
function add_scenix_to_InflowParam(elements, scenix::Int)
    for e in elements
        if e.typename == "ModeledInflowParam"
            e.value["ScenarioIndex"] = scenix
        end
    end
end

"""
Return copy of elements with replacement of Param with ModeledInflowParam 
if Param is embedded with a known station_id behind IFM_STATION_ID_KEY
The original Param must have Level and Profile keys, but this is how we 
normally set up inflow parameters anyway, so this should be ok.
"""
function copy_elements_iprogtype(elements, iprogtype, ifm_names, ifm_derivednames)
    if iprogtype == "ifm"
        elements1 = TuLiPa.DataElement[]
        for e in elements
            maybe_new_element = e
            if e.conceptname == TuLiPa.PARAM_CONCEPT
                if e.value isa Dict
                    if haskey(e.value, IFM_STATION_ID_KEY)
                        station_id = e.value[IFM_STATION_ID_KEY]
                        if (station_id in ifm_names) || (station_id in ifm_derivednames) 
                            maybe_new_element = TuLiPa.DataElement(e.conceptname, "ModeledInflowParam", e.instancename,
                                Dict("Level" => e.value["Level"], "HistoricalProfile" => e.value["Profile"], 
                                    IFM_STATION_ID_KEY => station_id,
                                    "phaseindays" => e.value["phaseindays"],
                                    "preddays" => e.value["preddays"]))
                        end
                    end
                end
            end
            push!(elements1, maybe_new_element)
        end
    else
        @assert iprogtype == "direct"
        elements1 = copy(elements)
    end
    @assert length(elements) == length(elements1)
    return elements1
end

function initialize_NN_model()
    error("Extension not loaded")
end

function basic_bucket_incl_states(p_, itp_Lday, itp_P, itp_T, t_out)
    error("Extension not loaded")
end

function NeuralODE_M100(p, norm_S0, norm_S1, norm_P, norm_T, 
        itp_Lday, itp_P, itp_T, t_out, ann; S_init = [0.0, 0.0])
    error("Extension not loaded")
end