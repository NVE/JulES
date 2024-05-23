

"""
Interface:
    estimate_S0(::InflowModel, ::ProbTime) -> S0
    predict(::InflowModel, ::ProbTime) -> (Q, S)
"""
abstract type InflowModel end

# TODO: Remove unused packages
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


mutable struct _InflowModelHandler{P, T1 <: TimeVector, T2 <: TimeVector, T3 <: TimeVector}
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

    obs_ndays::Int
    obs_P::Vector{Float32}   
    obs_T::Vector{Float32}  
    obs_Lday::Vector{Float32} 
    obs_timepoints::Vector{Float32}

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

        pred_timepoints = Vector{Float32}(1:pred_ndays)
        pred_P = zeros(Float32, pred_ndays)
        pred_T = zeros(Float32, pred_ndays)
        pred_Lday = zeros(Float32, pred_ndays)

        obs_usehist = false
        obs_timepoints = Vector{Float32}(1:obs_ndays)
        if isnothing(obs_P) || isnothing(obs_T) || isnothing(obs_Lday)
            obs_usehist = true
            obs_P = zeros(Float32, obs_ndays)
            obs_T = zeros(Float32, obs_ndays)
            obs_Lday = zeros(Float32, obs_ndays)
        end

        prev_t = nothing

        P = typeof(predictor)
        T1 = typeof(hist_P)
        T2 = typeof(hist_T)
        T3 = typeof(hist_Lday)

        return new{P, T1, T2, T3}(
            predictor, basin_area, hist_P, hist_T, hist_Lday,
            pred_ndays, pred_P, pred_T, pred_Lday, pred_timepoints,
            obs_usehist, obs_ndays, obs_P, obs_T, obs_Lday, obs_timepoints, prev_t)
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
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.obs_timepoints, m.obs_P, itp_method)
    itp_T = interpolate(m.obs_timepoints, m.obs_T, itp_method)
    itp_Lday = interpolate(m.obs_timepoints, m.obs_Lday, itp_method)
    (__, OED_sol) = m.predictor.predict(S0, G0, itp_Lday, itp_P, itp_T, m.obs_timepoints)

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
    itp_method = SteffenMonotonicInterpolation()
    itp_P = interpolate(m.pred_timepoints, m.pred_P, itp_method)
    itp_T = interpolate(m.pred_timepoints, m.pred_T, itp_method)
    itp_Lday = interpolate(m.pred_timepoints, m.pred_Lday, itp_method)
    (Q, OED_sol) = m.predictor.predict(S0, G0, itp_Lday, itp_P, itp_T, m.pred_timepoints)

    Q = Float64.(Q)
    Q .= Q ./ mm_per_m3s

    return Q
end

# ---- Common functions for NeuralOED model and Bucket -----
# TODO: link to source to give credit
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

struct _NeuralODEPredictor{P}
    nn_params::P
    mean_S::Float32
    mean_G::Float32
    mean_P::Float32
    mean_T::Float32
    std_S::Float32
    std_G::Float32
    std_P::Float32
    std_T::Float32
    function _NeuralODEPredictor(model_params)
        (nn_params, moments) = model_params
        return new{typeof(nn_params)}(nn_params, moments...)
    end
end

function predict(m::_NeuralODEPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    norm_S = (S) -> (S.-m.mean_S)./m.std_S
    norm_G = (G) -> (G.-m.mean_G)./m.std_G
    norm_P = (P) -> (P.-m.mean_P)./m.std_P
    norm_T = (T) -> (T.-m.mean_T)./m.std_T
    NeuralODE_M100(m.nn_params, norm_S, norm_G, norm_P, norm_T, 
        S0, G0, itp_Lday, itp_P, itp_T, timepoints)
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
function get_dist_ifm(input::AbstractJulESInput)
    names = get_inflow_names(input)
    cores = get_cores(input)
    N = length(cores)
    dist = Vector{Tuple{String, CoreId}}(undef, length(names))
    for (i, name) in enumerate(names)
        j = (i - 1) % N + 1
        dist[i] = (name, cores[j])
    end
    return dist
end


function solve_ifm(t)
    db = get_local_db()
    normfactors = get_ifm_normfactors(db)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            inflow_model = db.ifm[inflow_name]
            normalize_factor = normfactors[inflow_name]
            initial_state = estimate_initial_state(inflow_model, t)
            scenarios = get_scenarios(db.scenmod_ppp)
            for (scenix, scen) in enumerate(scenarios)
                scentime = get_scentphasein(t, scen, db.input)
                Q = predict(inflow_model, initial_state, scentime)
                Q .= Q .* normalize_factor
                start = getscenariotime(scentime)
                ix = [start + Day(i-1) for i in 1:length(Q)]    # TODO: Allocate this only once, then reuse
                db.ifm_output[inflow_name][scenix] = (ix, Q)
            end
        end
    end
end

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

function update_ifm_derived()
    db = get_local_db()
    ifm_weights = get_ifm_weights(db)
    for derived_name in keys(db.ifm_derived)
        weights = ifm_weights[derived_name]
        do_ix = true
        for (scenix, (derived_ix, derived_vals)) in db.ifm_derived[derived_name]
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

function create_ifm()
    db = get_local_db()
    elements = get_ifm_elements(db)
    modelobjects = getmodelobjects(elements)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            id = Id(TuLiPa.INFLOW_MODEL_CONCEPT, inflow_name)
            db.ifm[inflow_name] = modelobjects[id]
        end
    end
end

function includeModeledInflow!(::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    checkkey(lowlevel, elkey)

    deps = Id[]

    scenix = getdictvalue(value, "ScenarioIndex", Int, elkey)
    inflow_name = getdictvalue(value, "InflowName", String, elkey)

    hist_profile_name = getdictvalue(value, "HistoricalProfile", String, elkey)
    hist_profile_key = Id(TIMEVECTOR_CONCEPT,  hist_profile_name)
    push!(deps, hist_profile_key)

    level_name = getdictvalue(value, "Level", String, elkey)
    level_key = Id(TIMEVECTOR_CONCEPT,  level_name)
    push!(deps, level_key)

    haskey(lowlevel, hist_profile_key)   || return (false, deps)
    haskey(lowlevel, level_key)          || return (false, deps)

    hist_profile = lowlevel[hist_profile_key]
    level = lowlevel[level_key]

    # Creates an infinite timevector that refers to vectors stored in local db, 
    # which will be updated by JulES each step after running inflow models.
    # This way, model objects holding reference to such RHS term, will use updated 
    # prognosis from db when called upon by update!(prob, t)
    prognosis_profile = get_prognosis_from_local_db(inflow_name, scenix)

    steps = 1   # TODO: Is this correct?
    param = PrognosisSeriesParam(level, hist_profile, prognosis_profile, steps)

    isingoing = true
    id = Id(RHSTERM_CONCEPT, inflow_name)
    lowlevel[id] = BaseRHSTerm(id, param, isingoing)

    return (true, deps)
end

function get_prognosis_from_local_db(inflow_model, scenix)
    db = get_local_db()
    ifm_weights = get_ifm_weights(db)
    if haskey(ifm_weights, inflow_name)
        d = db.ifm_derived
    else
        d = db.ifm_output
    end
    if !haskey(d[inflow_name], scenix)
        d[inflow_name][scenix] = (DateTime[], Float64[])
    end
    (ix, vals) = d[inflow_name][scenix]
    return InfiniteTimeVector(ix, vals)
end
