"""
This code is copied from https://github.com/marv-in/HydroNODE
and has BSD 3-Clause License
"""
# TODO: Update License file
# TODO: Write license info in each source file

module IfmExt

using OrdinaryDiffEq
using ComponentArrays
using Interpolations
using JLD2
using Dates
using Statistics
using JulES


step_fct(x) = (tanh(5.0*x) + 1.0)*0.5
Ps(P, T, Tmin) = step_fct(Tmin-T)*P
Pr(P, T, Tmin) = step_fct(T-Tmin)*P
M(S0, T, Df, Tmax) = step_fct(T-Tmax)*step_fct(S0)*minimum([S0, Df*(T-Tmax)])
PET(T, Lday) = 29.8 * Lday * 0.611 * exp((17.3*T)/(T+237.3)) / (T + 273.2)
ET(S1, T, Lday, Smax) = step_fct(S1)*step_fct(S1-Smax)*PET(T,Lday) + step_fct(S1)*step_fct(Smax-S1)*PET(T,Lday)*(S1/Smax)
Qb(S1,f,Smax,Qmax) = step_fct(S1)*step_fct(S1-Smax)*Qmax + step_fct(S1)*step_fct(Smax-S1)*Qmax*exp(-f*(Smax-S1))
Qs(S1, Smax) = step_fct(S1)*step_fct(S1-Smax)*(S1-Smax)

function JulES.predict(m::JulES.TwoStateIfmHandler, u0::Vector{Float64}, t::JulES.TuLiPa.ProbTime)
    JulES.update_prediction_data(m, m.updater, t)

    # create interpolation input functions
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.data_pred.timepoints, m.data_pred.P, itp_method)
    itp_T = interpolate(m.data_pred.timepoints, m.data_pred.T, itp_method)
    itp_Lday = interpolate(m.data_pred.timepoints, m.data_pred.Lday, itp_method)

    (S0, G0) = u0
    (Q, __) = JulES.predict(m.predictor, S0, G0, itp_Lday, itp_P, itp_T, m.data_pred.timepoints)

    Q = Float64.(Q)

    Q .= Q .* m.m3s_per_mm

    return Q
end

function JulES.estimate_u0(m::JulES.TwoStateIfmHandler, t::JulES.TuLiPa.ProbTime)
    if isnothing(m.prev_t)
        JulES._initial_data_obs_update(m, t)
    else
        JulES._data_obs_update(m, t)
    end
    m.prev_t = t

    # do prediction from start of obs up until today
    (S0, G0) = (Float32(0), Float32(0))

    # create interpolation input functions
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.data_obs.timepoints, m.data_obs.P, itp_method)
    itp_T = interpolate(m.data_obs.timepoints, m.data_obs.T, itp_method)
    itp_Lday = interpolate(m.data_obs.timepoints, m.data_obs.Lday, itp_method)

    (__, OED_sol) = JulES.predict(m.predictor, S0, G0, itp_Lday, itp_P, itp_T, m.data_obs.timepoints)

    # extract states
    est_S0 = Float64(last(OED_sol[1, :]))
    est_G0 = Float64(last(OED_sol[2, :]))

    return [est_S0, est_G0]
end

function JulES.calculate_normalize_factor(ifm_model)
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
        P[i] = JulES.TuLiPa.getweightedaverage(ifm_model.handler.hist_P, start, JulES.ONEDAY_MS_TIMEDELTA)
        T[i] = JulES.TuLiPa.getweightedaverage(ifm_model.handler.hist_T, start, JulES.ONEDAY_MS_TIMEDELTA)
        Lday[i] = JulES.TuLiPa.getweightedaverage(ifm_model.handler.hist_Lday, start, JulES.ONEDAY_MS_TIMEDELTA)
    end

    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(timepoints_with_buffer, P, itp_method)
    itp_T = interpolate(timepoints_with_buffer, T, itp_method)
    itp_Lday = interpolate(timepoints_with_buffer, Lday, itp_method)
    Q, _ = JulES.predict(ifm_model.handler.predictor, 0, 0, itp_Lday, itp_P, itp_T, timepoints_with_buffer)
    Q = Float64.(Q)[timepoints_start:end]
    Q .= Q .* ifm_model.handler.m3s_per_mm
    return 1 / mean(Q)
end

function JulES.common_includeTwoStateIfm!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::JulES.TuLiPa.ElementKey, value::Dict)
    JulES.TuLiPa.checkkey(toplevel, elkey)

    model_params = JulES.TuLiPa.getdictvalue(value, "ModelParams", String, elkey)

    moments = nothing
    if haskey(value, "Moments")
        moments = JulES.TuLiPa.getdictvalue(value, "Moments", String, elkey)
    end

    hist_P = JulES.TuLiPa.getdictvalue(value, "HistoricalPercipitation", JulES.TuLiPa.TIMEVECTORPARSETYPES, elkey)
    hist_T = JulES.TuLiPa.getdictvalue(value, "HistoricalTemperature", JulES.TuLiPa.TIMEVECTORPARSETYPES, elkey)
    hist_Lday = JulES.TuLiPa.getdictvalue(value, "HistoricalDaylight", JulES.TuLiPa.TIMEVECTORPARSETYPES, elkey)

    ndays_pred = JulES.TuLiPa.getdictvalue(value, "NDaysPred", Real, elkey)
    try 
        ndays_pred = Int(ndays_pred)
        @assert ndays_pred >= 0
    catch e
        error("Value for key NDaysPred must be positive integer for $elkey")
    end

    basin_area = JulES.TuLiPa.getdictvalue(value, "BasinArea", Float64, elkey)

    deps = JulES.TuLiPa.Id[]

    all_ok = true

    (id, hist_P, ok) = JulES.TuLiPa.getdicttimevectorvalue(lowlevel, hist_P)    
    all_ok = all_ok && ok
    JulES.TuLiPa._update_deps(deps, id, ok)
    
    (id, hist_T, ok) = JulES.TuLiPa.getdicttimevectorvalue(lowlevel, hist_T)  
    all_ok = all_ok && ok
    JulES.TuLiPa._update_deps(deps, id, ok)

    (id, hist_Lday, ok) = JulES.TuLiPa.getdicttimevectorvalue(lowlevel, hist_Lday)  
    all_ok = all_ok && ok
    JulES.TuLiPa._update_deps(deps, id, ok)

    if all_ok == false
        return (false, deps)
    end

    # TODO: Maybe make this user input in future?
    updater = JulES.SimpleIfmDataUpdater()

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

    periodkey = JulES.TuLiPa.Id(JulES.TuLiPa.TIMEPERIOD_CONCEPT, "ScenarioTimePeriod")
    period = lowlevel[periodkey]
    scen_start = period["Start"]
    scen_stop  = period["Stop"]

    id = JulES.TuLiPa.getobjkey(elkey)
    toplevel[id] = Constructor(id, model_params, updater, basin_area, hist_P, hist_T, hist_Lday, 
        ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast, scen_start, scen_stop)

    return (true, deps)
end

function JulES.basic_bucket_incl_states(p_, itp_Lday, itp_P, itp_T, t_out)
    function exp_hydro_optim_states!(dS,S,ps,t)
        f, Smax, Qmax, Df, Tmax, Tmin = ps
        Lday = itp_Lday(t)
        P    = itp_P(t)
        T    = itp_T(t)
        Q_out = Qb(S[2],f,Smax,Qmax) + Qs(S[2], Smax)
        dS[1] = Ps(P, T, Tmin) - M(S[1], T, Df, Tmax)
        dS[2] = Pr(P, T, Tmin) + M(S[1], T, Df, Tmax) - ET(S[2], T, Lday, Smax) - Q_out
    end

    prob = ODEProblem(exp_hydro_optim_states!, p_[1:2], Float64.((t_out[1], maximum(t_out))))
    # sol = solve(prob, BS3(), u0=p_[1:2], p=p_[3:end], saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=ForwardDiffSensitivity())
    sol = solve(prob, BS3(), u0=p_[1:2], p=p_[3:end], saveat=t_out, reltol=1e-3, abstol=1e-3)
    Qb_ = Qb.(sol[2,:], p_[3], p_[4], p_[5])
    Qs_ = Qs.(sol[2,:], p_[4])
    Qout_ = Qb_.+Qs_
    return Qout_, sol
end

end