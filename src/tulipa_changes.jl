

struct ExternalHorizon <: Horizon
end

# For horizon-sync 
getchange(h::Horizon)::Dict
setchange(h::Horizon, change::Dict)