# Add TuLiPa when it is published as a package: /src/TuLiPa.jl"

using DistributedArrays

include("util.jl") # Various useful functions
include("scenarioreduction.jl") # Code for scenario reduction
include("prognosis.jl") # Code for price prognosis problems
include("stochastic.jl") # Code for stochastic subsystem problems
include("clearing.jl") # Code for market clearing problems