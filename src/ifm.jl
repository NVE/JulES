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

mutable struct TwoStateIfmHandler{PTyp, T1 <: TimeVector, T2 <: TimeVector, T3 <: TimeVector}
    predictor::PTyp

    basin_area::Float64

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

    function TwoStateIfmHandler(predictor, basin_area, hist_P, hist_T, hist_Lday, 
            ndays_pred, ndays_obs, ndays_forecast, data_obs, data_forecast)
        @assert ndays_pred > 0
        @assert ndays_obs > 0
        @assert ndays_forecast >= 0
        isnothing(data_obs) || @assert ndays_obs == length(data_obs.P)
        isnothing(data_forecast) || @assert ndays_forecast == length(data_forecast.P)
        (ndays_forecast > 0) && @assert !isnothing(data_forecast)
        data_pred = TwoStateIfmData(ndays_pred)
        P = typeof(predictor)
        T1 = typeof(hist_P)
        T2 = typeof(hist_T)
        T3 = typeof(hist_Lday)
        return new{P, T1, T2, T3}(predictor, basin_area, hist_P, hist_T, hist_Lday, 
            ndays_pred, ndays_obs, ndays_forecast, data_pred, data_obs, data_forecast, nothing)
    end
end
