# Exercise the production Benders loop with a solver double that invalidates its
# solution on every model update. This test needs only Julia's standard library.
module BendersCutConsistencyTests

using Test
using Dates
using Distributed

module TuLiPa
mutable struct Cuts
    lower_bound::Float64
    probabilities::Vector{Float64}
    constants::Vector{Float64}
    slopes::Matrix{Float64}
    scenconstants::Matrix{Float64}
    scenslopes::Array{Float64, 3}
    numcuts::Int
    cutix::Int
end

mutable struct Prob
    constants::Vector{Float64}
    slopes::Matrix{Float64}
    solved_constants::Vector{Float64}
    solved_slopes::Matrix{Float64}
    bounds::Vector{Float64}
    objectives::Vector{Float64}
    solves::Int
    valid::Bool
    full_updates::Int
    scenario_states::Vector{Vector{Float64}}
end

getnumcuts(cuts) = cuts.numcuts
getcutix(cuts) = cuts.cutix
getmaxcuts(cuts) = length(cuts.constants)
getfuturecostvarid(cuts) = nothing
setwarmstart!(prob, enabled) = nothing

function solve!(prob)
    prob.solves += 1
    prob.solved_constants .= prob.constants
    prob.solved_slopes .= prob.slopes
    prob.valid = true
end

function getvarvalue(prob, id, ix)
    @assert prob.valid "Cannot retrieve values after modifying the model"
    return prob.bounds[prob.solves]
end

function getoutgoingstates!(prob, states)
    @assert prob.valid
    states[1] = 10.0 * prob.solves
end

function clearcuts!(cuts)
    fill!(cuts.constants, cuts.lower_bound)
    fill!(cuts.slopes, 0.0)
    fill!(cuts.scenconstants, cuts.lower_bound)
    fill!(cuts.scenslopes, 0.0)
    cuts.numcuts = 0
    cuts.cutix = 0
end

function updatecutparameters!(prob, cuts)
    cuts.cutix = mod1(cuts.cutix + 1, getmaxcuts(cuts))
    cuts.numcuts = max(cuts.numcuts, cuts.cutix)
    cuts.constants[cuts.cutix] = cuts.scenconstants[1, cuts.cutix]
    cuts.slopes[cuts.cutix, :] .= cuts.scenslopes[1, cuts.cutix, :]
end

function updatecuts!(prob, cuts)
    prob.constants .= cuts.constants
    prob.slopes .= cuts.slopes
    prob.valid = false
    prob.full_updates += 1
end

function updatelastcut!(prob, cuts)
    prob.constants[cuts.cutix] = cuts.constants[cuts.cutix]
    prob.slopes[cuts.cutix, :] .= cuts.slopes[cuts.cutix, :]
    prob.valid = false
end

function clearcuts!(prob, cuts)
    clearcuts!(cuts)
    updatecuts!(prob, cuts)
end
end

const SubsystemIx = Int
const MainTiming = :timing
const DATABASE = Ref{Any}()
get_local_db() = DATABASE[]
get_settings(db) = db.settings

function solve_sp(scenix, subix, states)
    prob = get_local_db().mp[subix].prob
    push!(prob.scenario_states, copy(states))
    objective = prob.objectives[prob.solves + 1]
    return (scenix, objective, [-4.0], 7.0 + prob.solves)
end

function solve_sp_with_startreservoirs(scenix, subix, states)
    @assert get_local_db().mp[subix].prob.solves == 0 "Solved master states must be used"
    return solve_sp(scenix, subix, states)
end

# Load the actual loop, not a duplicate implementation of its control flow.
source = Meta.parseall(read(joinpath(@__DIR__, "..", "src", "prob_stoch.jl"), String))
definitions = map(source.args) do expr
    # Documented functions are wrapped in a Core.@doc macro call by the parser.
    expr isa Expr && expr.head == :macrocall ? last(expr.args) : expr
end
definition = only(filter(definitions) do expr
    expr isa Expr && expr.head == :function &&
        expr.args[1] isa Expr && expr.args[1].head == :call &&
        expr.args[1].args[1] == :solve_benders
end)
Core.eval(@__MODULE__, definition)

function setup(; numcuts=1, maxcuts=4, lower_bound=-100.0,
               bounds=[100.0], objectives=[100.0, 100.0])
    constants = fill(lower_bound, maxcuts)
    slopes = zeros(maxcuts, 1)
    for i in 1:numcuts
        constants[i] = -20.0 + i
        slopes[i, 1] = -Float64(i)
    end
    cuts = TuLiPa.Cuts(lower_bound, [1.0], copy(constants), copy(slopes),
        reshape(copy(constants), 1, maxcuts), reshape(copy(slopes), 1, maxcuts, 1),
        numcuts, numcuts)
    prob = TuLiPa.Prob(copy(constants), copy(slopes), copy(constants), copy(slopes),
        bounds, objectives, 0, false, 0, Vector{Float64}[])
    mp = (prob=prob, cuts=cuts, states=Dict(1 => 0.0), div=Dict(MainTiming => zeros(5)))
    DATABASE[] = (mp=Dict(1 => mp), dist_sp=[(1, 1, myid())],
        settings=Dict("problems" => Dict("stochastic" => Dict("reltol" => 0.001))))
    return mp
end

function check_consistent(mp)
    @test mp.prob.valid
    @test mp.cuts.constants == mp.prob.constants == mp.prob.solved_constants
    @test mp.cuts.slopes == mp.prob.slopes == mp.prob.solved_slopes
end

@testset "Benders cut consistency" begin
    @testset "Convergence on first reused master preserves all cut data" begin
        mp = setup(numcuts=2)
        before = deepcopy(mp.cuts)
        solve_benders(2, 1)
        check_consistent(mp)
        for field in fieldnames(TuLiPa.Cuts)
            @test getfield(mp.cuts, field) == getfield(before, field)
        end
        @test mp.prob.solves == 1
        @test mp.prob.full_updates == 0
        @test mp.prob.scenario_states == [[10.0]]
    end

    @testset "Continuing warm reuse resets every row" begin
        mp = setup(numcuts=2, bounds=[0.0, 100.0], objectives=fill(100.0, 3))
        solve_benders(2, 1)
        check_consistent(mp)
        @test mp.prob.solves == 2
        @test mp.prob.full_updates == 1
        @test mp.cuts.numcuts == mp.cuts.cutix == 1
        @test mp.cuts.constants == [8.0, -100.0, -100.0, -100.0]
        @test mp.prob.scenario_states == [[10.0], [20.0]]
    end

    @testset "Cold start must solve a master even with zero initial gap" begin
        mp = setup(numcuts=0, lower_bound=0.0, bounds=[0.0], objectives=[0.0, 0.0])
        solve_benders(1, 1)
        check_consistent(mp)
        @test mp.prob.solves == 1
        @test mp.cuts.numcuts == 1
        @test mp.cuts.constants[1] == 7.0
    end

    @testset "No-progress reset removes obsolete solver rows" begin
        mp = setup(bounds=[0.0, 0.0, 100.0], objectives=fill(100.0, 4))
        solve_benders(1, 1)
        check_consistent(mp)
        @test mp.prob.full_updates == 1
        @test mp.cuts.numcuts == 1
        @test mp.cuts.constants == [9.0, -100.0, -100.0, -100.0]
    end

    @testset "Circular replacement is committed only when continuing" begin
        mp = setup(numcuts=2, maxcuts=2, bounds=[0.0, 100.0], objectives=fill(100.0, 3))
        solve_benders(1, 1)
        check_consistent(mp)
        @test mp.cuts.cutix == 1
        @test mp.cuts.constants == [8.0, -18.0]
    end

    @testset "Iteration limit discards candidate cut" begin
        mp = setup(bounds=zeros(15), objectives=fill(100.0, 16))
        @test_logs (:warn, r"Benders did not converge") solve_benders(1, 1)
        check_consistent(mp)
        @test mp.prob.solves == 15
        @test mp.cuts.constants[mp.cuts.cutix] == 21.0
    end

    @testset "Convergence at iteration limit does not warn" begin
        mp = setup(bounds=vcat(zeros(14), 100.0), objectives=fill(100.0, 16))
        @test_logs solve_benders(1, 1)
        check_consistent(mp)
        @test mp.prob.solves == 15
    end
end

end