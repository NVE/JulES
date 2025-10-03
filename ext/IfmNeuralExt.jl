"""
This code is copied from https://github.com/marv-in/HydroNODE
and has BSD 3-Clause License
"""
# TODO: Update License file
# TODO: Write license info in each source file

module IfmNeuralExt

using DiffEqFlux
using SciMLSensitivity
using Optimization
using OptimizationOptimisers
using OptimizationBBO
using Zygote
using Lux
using Random
using CSV
using JulES

function JulES.initialize_NN_model()
    rng = Random.default_rng()
    NNmodel = Lux.Chain(Lux.Dense(4, 32, tanh), Lux.Dense(32,32, leakyrelu), Lux.Dense(32,32, leakyrelu), Lux.Dense(32,32, leakyrelu), 
    Lux.Dense(32,32, leakyrelu), Lux.Dense(32,5))
    p_NN_init, st_NN_init = Lux.setup(rng, NNmodel)
    NN_in_fct(x, p) = NNmodel(x,p,st_NN_init)[1]
    p_NN_init = ComponentArray(p_NN_init)
    return NN_in_fct, p_NN_init
end

function JulES.NeuralODE_M100(p, norm_S0, norm_S1, norm_P, norm_T, itp_Lday, itp_P, itp_T, t_out, ann; S_init = [0.0, 0.0])
    function NeuralODE_M100_core!(dS,S,p,t)
        Lday = itp_Lday(t)
        P    = itp_P(t)
        T    = itp_T(t)
        g = ann([norm_S0(S[1]), norm_S1(S[2]), norm_P(P), norm_T(T)],p)
        melting = relu(step_fct(S[1])*sinh(g[3]))
        dS[1] = relu(sinh(g[4])*step_fct(-T)) - melting
        dS[2] = relu(sinh(g[5])) + melting - step_fct(S[2])*Lday*exp(g[1])- step_fct(S[2])*exp(g[2])
    end
    prob = ODEProblem(NeuralODE_M100_core!, S_init, Float64.((t_out[1], maximum(t_out))), p)
    # sol = solve(prob, BS3(), dt=1.0, saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
    sol = solve(prob, BS3(), dt=1.0, saveat=t_out, reltol=1e-3, abstol=1e-3)
    P_interp = norm_P.(itp_P.(t_out))
    T_interp = norm_T.(itp_T.(t_out))
    S0_ = norm_S0.(sol[1,:])
    S1_ = norm_S1.(sol[2,:])
    Qout_ =  exp.(ann(permutedims([S0_ S1_ P_interp T_interp]),p)[2,:])
    return Qout_, sol
end

end