struct ClearingProblem
    prob::TuLiPa.Prob
    endstates::Dict{String, Float64} # startstates in next iteration
    div::Dict
end
# TODO: also preallocate clearingstates?

struct EndValueProblem
    prob::TuLiPa.Prob
    div::Dict
end

struct PricePrognosisProblem
    longprob::TuLiPa.Prob
    medprob::TuLiPa.Prob
    shortprob::TuLiPa.Prob
    nonstoragestates_short::Dict{TuLiPa.StateVariableInfo, Float64}
    div::Dict
end

mutable struct MasterProblem
    prob::TuLiPa.Prob
    cuts::TuLiPa.SimpleSingleCuts
    states::Dict{TuLiPa.StateVariableInfo, Float64}
    div::Dict
end

mutable struct ScenarioProblem # TODO: Should the others be mutable as well?
    prob::TuLiPa.Prob
    scenslopes::Vector{Float64}
    scenconstant::Float64
    div::Dict
end