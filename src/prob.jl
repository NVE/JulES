struct ClearingProblem
    prob::Prob
    endstates::Dict{String, Float64} # startstates in next iteration
    div::Dict
end
# TODO: also preallocate clearingstates?

struct EndValueProblem
    prob::Prob
    div::Dict
end

struct PricePrognosisProblem
    longprob::Prob
    medprob::Prob
    shortprob::Prob
    nonstoragestates_short::Dict{StateVariableInfo, Float64}
    div::Dict
end

mutable struct MasterProblem
    prob::Prob
    cuts::SimpleSingleCuts
    states::Dict{StateVariableInfo, Float64}
    div::Dict
end

mutable struct ScenarioProblem # TODO: Should the others be mutable as well?
    prob::Prob
    scenslopes::Vector{Float64}
    scenconstant::Float64
    div::Dict
end