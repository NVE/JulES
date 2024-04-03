struct ClearingProblem
    prob::Prob
    endstates::Dict{String, Float64} # startstates in next iteration
end
# TODO: also preallocate clearingstates?

struct EndValueProblem
    prob::Prob
end

struct PricePrognosisProblem
    longprob::Prob
    medprob::Prob
    shortprob::Prob
    nonstoragestates_short::Dict{StateVariableInfo, Float64}
end

struct MasterProblem
    prob::Prob
    cuts::SimpleSingleCuts
    states::Dict{StateVariableInfo, Float64}
end

struct ScenarioProblem
    prob::Prob
    scenslopes::Vector{Float64}
    scenconstant::Float64
end