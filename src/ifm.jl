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
        isnothing(data_obs) || @assert ndays_obs == length(data_obs.P)
        isnothing(data_forecast) || @assert ndays_forecast == length(data_forecast.P)
        (ndays_forecast > 0) && @assert !isnothing(data_forecast)
        data_pred = TwoStateIfmData(ndays_pred)
        m3s_per_mm = 1/((1000**3)/(basin_area*10**6)*86400)
        P = typeof(predictor)
        U = typeof(updater)
        T1 = typeof(hist_P)
        T2 = typeof(hist_T)
        T3 = typeof(hist_Lday)
        return new{P, T1, T2, T3}(predictor, basin_area, m3s_per_mm, hist_P, hist_T, hist_Lday, 
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
            m.data_obs.P[i] = getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA) * m.m3s_per_mm
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

    # use forecast if available
    m.ndays_forecast_used = 0
    if m.ndays_forecast > 0
        startix = length(m.data_forecast.P) - m.ndays_forecast + 1
        stopix = length(m.data_forecast.P) - m.ndays_forecast + ndays_update
        stopix = min(stopix, length(m.data_forecast.P))
        m.ndays_forecast_used = (stopix - startix) + 1
        j = m.ndays_obs - ndays_forecast_used + 1
        for i in startix:stopix
            m.data_obs.P[j] = m.data_forecast.P[i]
            m.data_obs.T[j] = m.data_forecast.T[i]
            m.data_obs.Lday[j] = m.data_forecast.Lday[i]
            j += 1
        end
    end

    # update possible remaining values
    ndays_remaining = ndays_update - m.ndays_forecast_used
    if ndays_remaining > 0
        start = getscenariotime(t) - Day(ndays_remaining + 1)
        i = m.ndays_obs - ndays_remaining + 1
        for __ in 1:ndays_remaining
            m.data_obs.P[i] = getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA) * m.m3s_per_mm
            m.data_obs.T[i] = getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
            m.data_obs.Lday[i] = getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
            start += Day(1) 
            i += 1
        end
    end
    @assert m.ndays_forecast >= 0   # TODO: remove validation

    # update struct state
    m.ndays_forecast -= m.ndays_forecast_used
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
    if m.ndays_forecast_used > 0
        ndays_before_estimate_u0_call = m.ndays_forecast + m.ndays_forecast_used
        startix = length(m.data_forecast.P) - ndays_before_estimate_u0_call + 1
        stopix = startix + m.ndays_forecast_used
        for (i, j) in enumerate(startix:stopix)
            m.data_pred.P[i] = m.data_forecast.P[j]
            m.data_pred.T[i] = m.data_forecast.T[j]
            m.data_pred.Lday[i] = m.data_forecast.Lday[j]
        end
    end
    for i in (ndays_forecast_used + 1):m.pred_ndays
        start = getscenariotime(t) + Day(i-1)
        m.data_pred.P[i] = getweightedaverage(m.hist_P, start, ONEDAY_MS_TIMEDELTA) * m.m3s_per_mm
        m.data_pred.T[i] = getweightedaverage(m.hist_T, start, ONEDAY_MS_TIMEDELTA)
        m.data_pred.Lday[i] = getweightedaverage(m.hist_Lday, start, ONEDAY_MS_TIMEDELTA)
    end
    return
end

