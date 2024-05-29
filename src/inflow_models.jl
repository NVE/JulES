"""
In this file we define:
- BucketInflowModel and NeuralOEDInflowModel

- Functions used in run_serial in connection with inflow models

- The ModeledInflowParam DataElement, which connects output from
  inflow models to model object inflow in optimization problems in JulES,
  through the local db

- The MixedInflowParam DataElement, which first uses PrognosisSeriesParam 
  for N days and thereafter uses ModeledInflowParam
"""

# --- BucketInflowModel and NeuralOEDInflowModel ---

"""
Common component of BucketInflowModel and NeuralOEDInflowModel
"""
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

function estimate_S0(m::_InflowModelHandler, t::ProbTime)
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

# ---- BucketInflowModel ---

function bucket_predict(p, S0, G0, itp_Lday, itp_P, itp_T, t_out)
    basic_bucket_incl_states([S0, G0, p...], itp_Lday, itp_P, itp_T, t_out)
end

struct _BucketPredictor{P}
    model_params::P
end
function predict(m::_BucketPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    bucket_predict(m.model_params, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
end

struct BucketInflowModel{H} <: AbstractInflowModel
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
                                        obs_P=obs_P,
                                        obs_T=obs_T,
                                        obs_Lday=obs_Lday)        
return new{typeof(handler)}(id, handler)
    end
end

estimate_S0(m::BucketInflowModel, t::ProbTime) = estimate_S0(m.handler, t)
predict(m::BucketInflowModel, S0, t::ProbTime) = predict(m.handler, S0, t)

function includeBucketInflowModel!(toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    _common_includeInflowModel!(BucketInflowModel, toplevel, lowlevel, elkey, value)
end

# ---- NeuralOEDInflowModel ---

struct _NeuralODEPredictor{P, NN}
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
    function _NeuralODEPredictor(model_params)
        (nn, __) = initialize_NN_model()
        (nn_params, moments) = model_params
        return new{typeof(nn_params)}(nn_params, nn, moments...)
    end
end

function predict(m::_NeuralODEPredictor, S0, G0, itp_Lday, itp_P, itp_T, timepoints)
    norm_S = (S) -> (S.-m.mean_S)./m.std_S
    norm_G = (G) -> (G.-m.mean_G)./m.std_G
    norm_P = (P) -> (P.-m.mean_P)./m.std_P
    norm_T = (T) -> (T.-m.mean_T)./m.std_T
    NeuralODE_M100(m.nn_params, norm_S, norm_G, norm_P, norm_T, 
        S0, G0, itp_Lday, itp_P, itp_T, timepoints, m.nn)
end

struct NeuralOEDInflowModel{H} <: AbstractInflowModel
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
                                        obs_P=obs_P,
                                        obs_T=obs_T,
                                        obs_Lday=obs_Lday)        
        return new{typeof(handler)}(id, handler)
    end
end

estimate_S0(m::NeuralOEDInflowModel, t::ProbTime) = estimate_S0(m.handler, t)
predict(m::NeuralOEDInflowModel, S0, t::ProbTime) = predict(m.handler, S0, t)

function includeNeuralOEDInflowModel!(toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    _common_includeInflowModel!(NeuralOEDInflowModel, toplevel, lowlevel, elkey, value)
end

function _common_includeInflowModel!(Constructor, toplevel::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    checkkey(toplevel, elkey)
    
    model_params = getdictvalue(value, "ModelParams", (Any, ), elkey)
    if model_params isa String
        model_params = JLD2.load_object(model_params)
    end

    hist_P = getdictvalue(value, "HistoricalPercipitation",   TIMEVECTORPARSETYPES, elkey)
    hist_T = getdictvalue(value, "HistoricalTemperature", TIMEVECTORPARSETYPES, elkey)
    hist_Lday = getdictvalue(value, "HistoricalDaylight", TIMEVECTORPARSETYPES, elkey)

    pred_ndays = getdictvalue(value, "NDaysPred", Int, elkey)
    obs_ndays = getdictvalue(value, "NDaysObs", Int, elkey)
    basin_area = getdictvalue(value, "BasinArea", Float64, elkey)

    if haskey(value, "ObservedPercipitation")
        obs_P = getdictvalue(value, "ObservedPercipitation",   Vector{Real}, elkey)
    else
        obs_P = nothing
    end
    if haskey(value, "ObservedTemperature")
        obs_T = getdictvalue(value, "ObservedTemperature",   Vector{Real}, elkey)
    else
        obs_T = nothing
    end
    if haskey(value, "ObservedDaylight")
        obs_Lday = getdictvalue(value, "ObservedDaylight",   Vector{Real}, elkey)
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


# --- Functions used in run_serial in connection with inflow models ---

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

function solve_ifm(t)
    db = get_local_db()
    normfactors = get_ifm_normfactors(db)
    scenarios = get_scenarios(db.scenmod_ppp)
    for (inflow_name, core) in db.dist_ifm
        if core == db.core
            inflow_model = db.ifm[inflow_name]
            normalize_factor = normfactors[inflow_name]
            S0 = estimate_S0(inflow_model, t)
            for (scenix, scen) in enumerate(scenarios)
                scentime = get_scentphasein(t, scen, db.input)
                Q = predict(inflow_model, S0, scentime) # TODO: Return (Q, S) instead
                Q .= Q .* normalize_factor
                start = getscenariotime(scentime)
                ix = [start + Day(i-1) for i in 1:length(Q)]    # TODO: Allocate this only once, then reuse
                db.ifm_output[inflow_name][scenix] = (ix, Q)
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

# ---- The ModeledInflowParam DataElement -----

"""
ModeledInflowParam is actually a PrognosisSeriesParam under-the-hood,
but where the prognosis part is created and managed by JulES and not 
given as user input
"""
function includeModeledInflowParam!(::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
    checkkey(lowlevel, elkey)

    deps = Id[]

    # Not part of user input 
    # This info is added by JulES with add_scenix_to_ModeledInflowParam
    # See e.g. prob_stoch.get_elements_with_horizons
    scenix = getdictvalue(value, "ScenarioIndex", Int, elkey)

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

    # Creates an InfiniteTimeVector that refers to vectors stored in local db, 
    # which will be updated by JulES each step after running inflow models.
    # This way, model objects holding reference to such Param, will use updated 
    # prognosis from db when called upon by update!(prob, t)
    db = get_local_db()
    replacemap = get_ifm_replacemap(db.input)
    haskey(replacemap, elkey.instancename) || error("Instance name not found in replacemap for $elkey")
    inflow_name = replacemap[elkey.instancename]
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
    prognosis_profile = InfiniteTimeVector(ix, vals)    

    steps = 1   # use prognosis 100% when it applies

    lowlevel[getobjkey(elkey)] = PrognosisSeriesParam(level, hist_profile, prognosis_profile, steps)

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
# function includeMixedInflowParam!(::Dict, lowlevel::Dict, elkey::ElementKey, value::Dict)
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
        elements1 = DataElement[]
        for e in elements
            if e.typename == "PrognosisSeriesParam" && haskey(ifm_replacemap, e.instancename)
                new_e = DataElement(e.conceptname, "ModeledInflowParam", e.instancename,
                    Dict("Level" => e.value["Level"], "HistoricalProfile" => e.value["Profile"]))
                push!(elements1, new_e)
            else
                push!(elements1, e)
            end
        end
    elseif startswith(iprogtype, "mix")
        # TODO: Validate ndays > 0 in constructor of DefaultJulESInput
        ndays = parse(Int, iprogtype[4:end])
        elements1 = DataElement[]
        for e in elements
            if e.typename == "PrognosisSeriesParam" && haskey(ifm_replacemap, e.instancename)
                new_value = copy(e.value::Dict)
                new_value["ndays"] = ndays
                new_e = DataElement(e.conceptname, "MixedInflowParam", e.instancename, new_value)
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

# Register extentions to TuLiPa input system
TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(ABSTRACT_INFLOW_MODEL, "BucketInflowModel")] = includeBucketInflowModel!
TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(ABSTRACT_INFLOW_MODEL, "NeuralOEDInflowModel")] = includeNeuralOEDInflowModel!
TuLiPa.INCLUDEELEMENT[TuLiPa.TypeKey(TuLiPa.PARAM_CONCEPT, "ModeledInflowParam")] = includeModeledInflowParam!
