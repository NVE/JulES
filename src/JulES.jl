module JulES

macro debugtime(msg, expr)
    quote
        local stats = Base.@timed $(esc(expr))
        @debug $msg elapsed_s = round(stats.time, digits=3) bytes = stats.bytes gctime_s = round(stats.gctime, digits=3)
        stats.value
    end
end

import TuLiPa

using Distributed
using Dates
using Statistics
using Clustering
using Distributions
using DataFrames
using JSON
using YAML
using HDF5
using Logging
# Used by ifm
#using ComponentArrays
#using Interpolations
#using JLD2

# Used by Nerual inflow model but not HBV
# using DiffEqFlux
# using SciMLSensitivity
# using Optimization
# using OptimizationOptimisers
# using OptimizationBBO
# using Zygote

include("python_logger.jl")
include("abstract_types.jl")
include("dimension_types.jl")
include("ifm.jl")
include("generic_io.jl")
include("io.jl")
include("prob_cp.jl")
include("prob_evp.jl")
include("prob_ppp.jl")
include("prob_stoch.jl")
include("prob_util.jl")
include("local_db.jl")
include("run_serial.jl")
include("scenariomodelling.jl")
include("run_jules_wrapper.jl")

end