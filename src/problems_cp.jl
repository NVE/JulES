struct ClearingProblem
    prob
    states
    endstates # startstates in next iteration
end
get_prob(cp::ClearingProblem) = cp.prob
get_states(cp::ClearingProblem) = cp.states
get_endstates(cp::ClearingProblem) = cp.endstates

function create_cp(db::LocalDB)
end

get_startstates_from_cp() = get_endstates(get_cp(get_local_db()))

# TODO: Rename to update_startstates
function get_startstates!(clearing::Prob, detailedrescopl::Dict, enekvglobaldict::Dict, startstates::Dict{String, Float64})
    startstates_ = get_states(getobjects(clearing))
    getoutgoingstates!(clearing, startstates_)
    
    for var in keys(startstates_)
        value = round(startstates_[var], digits=10) # avoid approx 0 negative values, ignored by solvers so no problem?
        startstates[getinstancename(first(getvarout(var)))] = value
    end

    for area in Set(values(detailedrescopl))
        startstates["Reservoir_" * area * "_hydro_reservoir"] = 0.0
    end
    for res in keys(detailedrescopl)
        resname = "Reservoir_" * res
        areaname = "Reservoir_" * detailedrescopl[res] * "_hydro_reservoir"
        startstates[areaname] += startstates[resname] * enekvglobaldict[res]
    end

    # Avoid reservoirs being filled more than max, gives infeasible solution
    # - If aggregated reservoir capacity is lower than the sum capacities
    # - If reservoir is full in model, numerical tolerance can bring variable value slightly over cap
    # - TODO: Add warning/logging if this happens
    for resname in keys(startstates)
        resmax = resname * "_max"
        if haskey(startstates, resmax)
            if startstates[resname] > startstates[resmax]
                startstates[resname] = startstates[resmax]
            end
        end
    end
end