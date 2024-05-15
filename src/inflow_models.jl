

"""
Interface:
    estimate_S0(::InflowModel, ::ProbTime) -> S0
    predict(::InflowModel, ::ProbTime) -> (Q, S)
"""
abstract type InflowModel end

using CSV
using DataFrames, Dates, Statistics
using OrdinaryDiffEq, DiffEqFlux, Lux
using ComponentArrays
using SciMLSensitivity
using Optimization, OptimizationOptimisers, OptimizationBBO
using Zygote
using Interpolations
using Random
using JLD2


mutable struct _InflowModelHandler{P, T1 <: TimeVector, T2 <: TimeVector, T3 <: TimeVector, F1, F2, F3, F4, F5, F6}
    predictor::P

    basin_area::Float64

    hist_P::T1
    hist_T::T2
    hist_Lday::T3

    pred_ndays::Int
    pred_P::Vector{Float32}   
    pred_T::Vector{Float32}  
    pred_Lday::Vector{Float32} 
    pred_timepoints::Vector{Float32}
    pred_itp_P::F1
    pred_itp_T::F2
    pred_itp_Lday::F3

    obs_ndays::Int
    obs_P::Vector{Float32}   
    obs_T::Vector{Float32}  
    obs_Lday::Vector{Float32} 
    obs_timepoints::Vector{Float32}
    obs_itp_P::F4
    obs_itp_T::F5
    obs_itp_Lday::F6

    prev_t::Union{ProbTime, Nothing}

    function _InflowModelHandler(predictor, basin_area, hist_P, hist_T, hist_Lday, 
                               pred_ndays, obs_ndays;
                               obs_P::Union{Nothing, Vector{Float32}}=nothing,
                               obs_T::Union{Nothing, Vector{Float32}}=nothing,
                               obs_Lday::Union{Nothing, Vector{Float32}}=nothing)
        @assert obs_ndays >= 200
        @assert pred_ndays >= 1
        isnothing(obs_P) || @assert len(obs_P) == obs_ndays
        isnothing(obs_T) || @assert len(obs_T) == obs_ndays
        isnothing(obs_Lday) || @assert len(obs_Lday) == obs_ndays

        itp_method = SteffenMonotonicInterpolation()

        pred_timepoints = Vector{Float32}(1:pred_ndays)
        pred_P = zeros(Float32, pred_ndays)
        pred_T = zeros(Float32, pred_ndays)
        pred_Lday = zeros(Float32, pred_ndays)
        pred_itp_P = interpolate(pred_timepoints, pred_P, itp_method)
        pred_itp_T = interpolate(pred_timepoints, pred_T, itp_method)
        pred_itp_Lday = interpolate(pred_timepoints, pred_Lday, itp_method)

        obs_usehist = false
        obs_timepoints = Vector{Float32}(1:obs_ndays)
        if isnothing(obs_P) || isnothing(obs_T) || isnothing(obs_Lday)
            obs_usehist = true
            obs_P = zeros(Float32, obs_ndays)
            obs_T = zeros(Float32, obs_ndays)
            obs_Lday = zeros(Float32, obs_ndays)
        end
        obs_itp_P = interpolate(obs_timepoints, obs_P, itp_method)
        obs_itp_T = interpolate(obs_timepoints, obs_T, itp_method)
        obs_itp_Lday = interpolate(obs_timepoints, obs_Lday, itp_method)

        prev_t = nothing

        P = typeof(predictor)
        T1 = typeof(hist_P)
        T2 = typeof(hist_T)
        T3 = typeof(hist_Lday)
        F1 = typeof(pred_itp_P)
        F2 = typeof(pred_itp_T)
        F3 = typeof(pred_itp_Lday)        
        F4 = typeof(obs_itp_P)
        F5 = typeof(obs_itp_T)
        F6 = typeof(obs_itp_Lday)

        return new{P, T1, T2, T3, F1, F2, F3, F4, F5, F6}(
            predictor, basin_area, hist_P, hist_T, hist_Lday,
            pred_ndays, pred_P, pred_T, pred_Lday,
            pred_timepoints, pred_itp_P, pred_itp_T, pred_itp_Lday,
            obs_usehist, obs_ndays, obs_P, obs_T, obs_Lday, 
            obs_timepoints, obs_itp_P, obs_itp_T, obs_itp_Lday, prev_t)
    end
end

function estimate_initial_state(m::_InflowModelHandler, t::ProbTime)
    mm_per_m3s = ((1000**3)/(m.basin_area*10**6)*86400)

    is_first_step = isnothing(m.prev_t)

    if is_first_step
        if m.obs_usehist
            m.obs_usehist = false
            delta = MSTimeDelta(Day(1))
            @inbounds for i in 1:m.obs_ndays
                ndays_back = m.obs_ndays + 1 - i
                start = getdatatime(t) - Day(ndays_back)
                m.obs_P[i] = getweightedaverage(m.hist_P, start, delta) / mm_per_m3s
                m.obs_T[i] = getweightedaverage(m.hist_T, start, delta)
                m.obs_Lday[i] = getweightedaverage(m.hist_Lday, start, delta)
            end
        else
            # just use observations from user
        end
    else
        # update observations since last time
        N = min(m.obs_ndays,
                max(Day(getdatatime(t) - getdatatime(m.prev_t)).value, 
                    Day(getscenariotime(t) - getscenariotime(m.prev_t)).value))
        
        # TODO: Use @inbounds 

        # shift backwards to make room for N new values
        for i in m.obs_ndays:-1:(N+1)
            m.obs_P[i-N] = m.obs_P[i]
        end
        # update the N last values
        for ndays_back in 1:N
            i = m.obs_ndays + 1 - ndays_back
            start = getscenariotime(t) - Day(ndays_back)
            m.obs_P[i] = getweightedaverage(m.hist_P, start, delta) / mm_per_m3s
            m.obs_T[i] = getweightedaverage(m.hist_T, start, delta)
            m.obs_Lday[i] = getweightedaverage(m.hist_Lday, start, delta)
        end
    end
    m.prev_t = t

    # do hindcast up until today
    (S0, G0) = (Float32(0), Float32(0)) 
    (__, OED_sol) = m.predictor.predict(S0, G0, m.obs_itp_Lday, m.obs_itp_P, m.obs_itp_T, m.obs_timepoints)

    # extract states
    est_S0 = last(OED_sol[1, :])
    est_G0 = last(OED_sol[2, :])

    return (est_S0, est_G0)
end

function predict(m::_InflowModelHandler, initial_State, t::ProbTime)
    # update prediction data
    # TODO: Different update methods here? E.g. to enable use of historical autocorrelation to phase in uncertainty?
    delta = MSTimeDelta(Day(1))
    mm_per_m3s = ((1000**3)/(m.basin_area*10**6)*86400)
    @inbounds for i in 1:m.pred_ndays
        start = getscenariotime(t) + Day(i-1)
        m.pred_P[i] = getweightedaverage(m.hist_P, start, delta) / mm_per_m3s
        m.pred_T[i] = getweightedaverage(m.hist_T, start, delta)
        m.pred_Lday[i] = getweightedaverage(m.hist_Lday, start, delta)
    end

    (S0, G0) = initial_State
    (__, OED_sol) = m.predictor.predict(S0, G0, m.pred_itp_Lday, m.pred_itp_P, m.pred_itp_T, m.pred_timepoints)

    Q = Float64.(Q)
    Q .= Q ./ mm_per_m3s

    return Q
end

# ---- Common functions for NeuralOED model and Bucket -----
step_fct(x) = (tanh(5.0*x) + 1.0)*0.5
Ps(P, T, Tmin) = step_fct(Tmin-T)*P
Pr(P, T, Tmin) = step_fct(T-Tmin)*P
M(S0, T, Df, Tmax) = step_fct(T-Tmax)*step_fct(S0)*minimum([S0, Df*(T-Tmax)])
PET(T, Lday) = 29.8 * Lday * 0.611 * exp((17.3*T)/(T+237.3)) / (T + 273.2)
ET(S1, T, Lday, Smax) = step_fct(S1)*step_fct(S1-Smax)*PET(T,Lday) + step_fct(S1)*step_fct(Smax-S1)*PET(T,Lday)*(S1/Smax)
Qb(S1,f,Smax,Qmax) = step_fct(S1)*step_fct(S1-Smax)*Qmax + step_fct(S1)*step_fct(Smax-S1)*Qmax*exp(-f*(Smax-S1))
Qs(S1, Smax) = step_fct(S1)*step_fct(S1-Smax)*(S1-Smax)


# ---- BucketInflowModel ---

function bucket_predict(p, S0, G0, itp_Lday, itp_P, itp_T, t_out)

    function exp_hydro_optim_states!(dS,S,ps,t)
        f, Smax, Qmax, Df, Tmax, Tmin = ps
        Lday = itp_Lday(t)
        P    = itp_P(t)
        T    = itp_T(t)
        Q_out = Qb(S[2],f,Smax,Qmax) + Qs(S[2], Smax)
        dS[1] = Ps(P, T, Tmin) - M(S[1], T, Df, Tmax)
        dS[2] = Pr(P, T, Tmin) + M(S[1], T, Df, Tmax) - ET(S[2], T, Lday, Smax) - Q_out
    end

    prob = ODEProblem(exp_hydro_optim_states!, [S0, G0], Float64.((t_out[1], maximum(t_out))))
    sol = solve(prob, BS3(), u0=[S0, G0], p=p, saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=ForwardDiffSensitivity())
    Qb_ = Qb.(sol[2,:], p[1], p[2], p[3])
    Qs_ = Qs.(sol[2,:], p[2])
    Qout_ = Qb_.+Qs_
    return Qout_, sol
end

struct _BucketPredictor{P}
    model_params::P
end
function predict(m::_BucketPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    bucket_predict(m.model_params, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
end

struct BucketInflowModel{H} <: InflowModel
    id::Id
    handler::H

    function BucketInflowModel(id, model_params, basin_area, hist_P, hist_T, hist_Lday, 
                                pred_ndays, obs_ndays;
                                obs_P::Union{Nothing, Vector{Float32}}=nothing,
                                obs_T::Union{Nothing, Vector{Float32}}=nothing,
                                obs_Lday::Union{Nothing, Vector{Float32}}=nothing)
        predictor = _BucketPredictor(model_params)
        handler = _InflowModelHandler(predictor, basin_area, hist_P, hist_T, hist_Lday, 
                                        pred_ndays, obs_ndays;
                                        obs_P::Union{Nothing, Vector{Float32}}=nothing,
                                        obs_T::Union{Nothing, Vector{Float32}}=nothing,
                                        obs_Lday::Union{Nothing, Vector{Float32}}=nothing)        
        return new{typeof(handler)}(id, handler)
    end
end

function includeBucketInflowModel!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    _common_includeInflowModel!(BucketInflowModel, toplevel, lowlevel, elkey, value)
end

# ---- NeuralOEDInflowModel ---

function NeuralODE_M100(p, norm_S0, norm_S1, norm_P, norm_T, itp_Lday, itp_P, itp_T, t_out, ann; S_init = [0.0, 0.0])
    function NeuralODE_M100_core!(dS,S,p,t)
        Lday = itp_Lday(t)
        P    = itp_P(t)
        T    = itp_T(t)
        g = ann([norm_S0(S[1]), norm_S1(S[2]), norm_P(P), norm_T(T)],p)
        melting = relu(step_fct(S[1])*sinh(g[3]))
        dS[1] = relu(sinh(g[4])*step_fct(-T)) - melting
        dS[2] = relu(sinh(g[5])) + melting - step_fct(S[2])*Lday*exp(g[1])- step_fct(S[2])*exp(g[2])
    end
    prob = ODEProblem(NeuralODE_M100_core!, S_init, Float64.((t_out[1], maximum(t_out))), p)
    sol = solve(prob, BS3(), dt=1.0, saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
    P_interp = norm_P.(itp_P.(t_out))
    T_interp = norm_T.(itp_T.(t_out))
    S0_ = norm_S0.(sol[1,:])
    S1_ = norm_S1.(sol[2,:])
    Qout_ =  exp.(ann(permutedims([S0_ S1_ P_interp T_interp]),p)[2,:])
    return Qout_, sol
end

struct _NeuralODEPredictor{P, F1, F2, F3, F4}
    nn_params::P
    norm_S::F1
    norm_G::F2
    norm_P::F3
    norm_T::F4
    function _NeuralODEPredictor(model_params)
        (nn_params, moments) = model_params
        (mean_S, std_S, mean_G, std_G, mean_P, std_P, mean_T, std_T) = moments    
        norm_S = (S) -> (S.-mean_S)./std_S
        norm_G = (G) -> (G.-mean_G)./std_G
        norm_P = (P) -> (P.-mean_P)./std_P
        norm_T = (T) -> (T.-mean_T)./std_T
        P = typeof(nn_params)
        F1 = typeof(norm_S)
        F2 = typeof(norm_G)
        F3 = typeof(norm_P)
        F4 = typeof(norm_T)
        return new{P, F1, F2, F3, F4}(nn_params, norm_S, norm_G, norm_P, norm_T)

end
function predict(m::_NeuralODEPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    NeuralODE_M100(m.nn_params, m.norm_S, m.norm_G, m.norm_P, m.norm_T, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
end

struct NeuralOEDInflowModel{H} <: InflowModel
    id::Id
    handler::H

    function NeuralOEDInflowModel(id, model_params, basin_area, hist_P, hist_T, hist_Lday, 
                                pred_ndays, obs_ndays;
                                obs_P::Union{Nothing, Vector{Float32}}=nothing,
                                obs_T::Union{Nothing, Vector{Float32}}=nothing,
                                obs_Lday::Union{Nothing, Vector{Float32}}=nothing)
        predictor = _NeuralODEPredictor(model_params)
        handler = _InflowModelHandler(predictor, basin_area, hist_P, hist_T, hist_Lday, 
                                        pred_ndays, obs_ndays;
                                        obs_P::Union{Nothing, Vector{Float32}}=nothing,
                                        obs_T::Union{Nothing, Vector{Float32}}=nothing,
                                        obs_Lday::Union{Nothing, Vector{Float32}}=nothing)        
        return new{typeof(handler)}(id, handler)
    end
end

function includeNeuralOEDInflowModel!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    _common_includeInflowModel!(NeuralOEDInflowModel, toplevel, lowlevel, elkey, value)
end

function _common_includeInflowModel!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    checkkey(toplevel, elkey)
    
    model_params = getdictvalue(value, "ModelParams", (Any, ), elkey)
    if model_params isa String
        model_params = JLD2.load_object(model_params)
    end

    hist_P = getdictvalue(value, "Historical_Percipitation",   TIMEVECTORPARSETYPES, elkey)
    hist_T = getdictvalue(value, "Historical_Temperature", TIMEVECTORPARSETYPES, elkey)
    hist_Lday = getdictvalue(value, "Historical_Daylight", TIMEVECTORPARSETYPES, elkey)

    pred_ndays = getdictvalue(value, "NDays_Pred", Int, elkey)
    obs_ndays = getdictvalue(value, "NDays_Obs", Int, elkey)
    basin_area = getdictvalue(value, "Basin_Area", Float64, elkey)

    if haskey(value, "Observed_Percipitation")
        obs_P = getdictvalue(value, "Observed_Percipitation",   Vector{Real}, elkey)
    else
        obs_P = nothing
    end
    if haskey(value, "Observed_Temperature")
        obs_T = getdictvalue(value, "Observed_Temperature",   Vector{Real}, elkey)
    else
        obs_T = nothing
    end
    if haskey(value, "Observed_Daylight")
        obs_Lday = getdictvalue(value, "Observed_Daylight",   Vector{Real}, elkey)
    else
        obs_Lday = nothing
    end

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

    id = getobjkey(elkey)
    toplevel[id] = Constructor(id, model_params, basin_area, hist_P, hist_T, hist_Lday, 
                               pred_ndays, obs_ndays; 
                               obs_P=obs_P, obs_T=obs_T, obs_Lday=Obs_Lday)

    return (true, deps)
end




"""
How inflow_models are distributed on cores initially
"""
function get_dist_inflow_models(input::AbstractJulESInput)
    cores = get_cores(input)
    N = length(cores)
    names = get_inflow_model_names(input)
    dist = Vector{Tuple{String, CoreId}}(undef, length(names))
    for (i, name) in enumerate(names)
        j = (i - 1) % N + 1
        dist[i] = (name, cores[j])
    end
    return dist
end

function solve_inflow_models(t, stepnr)
    db = get_local_db()
    for (station, core) in db.dist_inflow_models
        if core == db.core
            inflow_model = db.inflow_models[station]
            (stored_stepnr, d) = db.inflow_prognosis[station]
            initial_state = estimate_initial_state(inflow_model, t)
            scenarios = get_scenarios(db.scenmod_ppp)
            for (scenix, scen) in enumerate(scenarios)
                scentime = get_scentphasein(t, scen, db.input)
                d[scenix] = predict(inflow_model, initial_state, scentime)
            end
            db.inflow_prognosis[station] = (stepnr, d)
        end
    end
end
