"""
We need to transfer master-horizons from one core to other cores, 
and these should not do anything in update!(horizon, t),
because they are updated as part of the data transfer, as the true
horizon-update have already taken place in the master-horizon.
"""
struct ExternalHorizon{H <: Horizon} <: Horizon
    subhorizon::H
    function ExternalHorizon(h::Horizon)
        @assert !(h isa ExternalHorizon)
        new{typeof(h)}(h)
    end
end

# Forwarded methods
isadaptive(h::ExternalHorizon) = isadaptive(h.subhorizon)
getnumperiods(h::ExternalHorizon) = getnumperiods(h.subhorizon)
getstartduration(h::ExternalHorizon, t::Int) = getstartduration(h.subhorizon, t)
getendperiodfromduration(h::ExternalHorizon, d::Millisecond) = getendperiodfromduration(h.subhorizon, d)
getduration(h::ExternalHorizon) = getduration(h.subhorizon)
gettimedelta(h::ExternalHorizon, t::Int) = gettimedelta(h.subhorizon, t)
hasoffset(h::ExternalHorizon) = hasoffset(h.subhorizon)
getoffset(h::ExternalHorizon) = getoffset(h.subhorizon)
getstarttime(h::ExternalHorizon, t::Int, start::ProbTime) = getstarttime(h.subhorizon, t, start)
getsubperiods(coarse::ExternalHorizon, fine::Horizon, coarse_t::Int) = getsubperiods(coarse.subhorizon, fine, coarse_t)
getsubperiods(coarse::Horizon, fine::ExternalHorizon, coarse_t::Int) = getsubperiods(coarse, fine.subhorizon, coarse_t)
getsubperiods(coarse::ExternalHorizon, fine::ExternalHorizon, coarse_t::Int) = getsubperiods(coarse.subhorizon,fine.subhorizon,coarse_t)
hasconstantdurations(h::ExternalHorizon) = hasconstantdurations(h.subhorizon)
mayshiftfrom(h::ExternalHorizon, t::Int) = mayshiftfrom(h.subhorizon, t)
mustupdate(h::ExternalHorizon, t::Int) = mustupdate(h.subhorizon, t)

build!(h::ExternalHorizon, p::Prob) = build!(h.subhorizon, p)  # TODO: Should do nothing?

# Specialized methods
update!(::ExternalHorizon, ::ProbTime) = nothing

# TODO: Change this so it does not have to start in period 1
"""
In JulES, we would like subsystem models to use same horizon 
as price prognosis models, but not neccesary the whole horizon. For many systems, 
the first 2-3 years would be sufficiently long horizon. ShortendHorizon meets
this need, as it wraps another horison, and only use some of the first periods.
"""
struct ShortendHorizon{H <: Horizon} <: Horizon
    subhorizon::H
    n:Int
    function ShortendHorizon(h::Horizon, n::Int)
        @assert 0 < n <= getnumperiods(h) 
        new{typeof(h)}(h, n)
    end
end

# Forwarded methods
isadaptive(h::ShortendHorizon) = isadaptive(h.subhorizon)
getstartduration(h::ShortendHorizon, t::Int) = getstartduration(h.subhorizon, t)
getendperiodfromduration(h::ShortendHorizon, d::Millisecond) = getendperiodfromduration(h.subhorizon, d)
getduration(h::ShortendHorizon) = getduration(h.subhorizon)
gettimedelta(h::ShortendHorizon, t::Int) = gettimedelta(h.subhorizon, t)
hasoffset(h::ShortendHorizon) = hasoffset(h.subhorizon)
getoffset(h::ShortendHorizon) = getoffset(h.subhorizon)
getstarttime(h::ShortendHorizon, t::Int, start::ProbTime) = getstarttime(h.subhorizon, t, start)
getsubperiods(coarse::ShortendHorizon, fine::Horizon, coarse_t::Int) = getsubperiods(coarse.subhorizon, fine, coarse_t)
getsubperiods(coarse::Horizon, fine::ShortendHorizon, coarse_t::Int) = getsubperiods(coarse, fine.subhorizon, coarse_t)
getsubperiods(coarse::ShortendHorizon, fine::ShortendHorizon, coarse_t::Int) = getsubperiods(coarse.subhorizon,fine.subhorizon,coarse_t)
hasconstantdurations(h::ShortendHorizon) = hasconstantdurations(h.subhorizon)
mustupdate(h::ShortendHorizon, t::Int) = mustupdate(h.subhorizon, t)
build!(h::ShortendHorizon, p::Prob) = build!(h.subhorizon, p)

# Specialized methods
getnumperiods(h::ShortendHorizon) = h.n

function mayshiftfrom(h::ShortendHorizon, t::Int)
    (future_t, ok) = mayshiftfrom(h.subhorizon, t)
    if ok && future_t > h.n
        return (future_t, false)
    end 
    return (future_t, ok)
end

update!(::ShortendHorizon, ::ProbTime) = nothing


# New horizon interface functions in 
# order to keep remote copies of self in sync
getchanges(::Horizon) = error()
setchanges(::Horizon, changes::Dict) = error()


getchanges(::SequentialHorizon) = Dict()
setchanges(::SequentialHorizon, changes::Dict) = nothing
