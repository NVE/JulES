"""
This code is copied from https://github.com/marv-in/HydroNODE
and has BSD 3-Clause License
"""
# TODO: Update License file
# TODO: Write license info in each source file

function initialize_NN_model()
    rng = Random.default_rng()
    NNmodel = Lux.Chain(Lux.Dense(4, 32, tanh), Lux.Dense(32,32, leakyrelu), Lux.Dense(32,32, leakyrelu), Lux.Dense(32,32, leakyrelu), 
    Lux.Dense(32,32, leakyrelu), Lux.Dense(32,5))
    p_NN_init, st_NN_init = Lux.setup(rng, NNmodel)
    NN_in_fct(x, p) = NNmodel(x,p,st_NN_init)[1]
    p_NN_init = ComponentArray(p_NN_init)
    return NN_in_fct, p_NN_init
end

step_fct(x) = (tanh(5.0*x) + 1.0)*0.5
Ps(P, T, Tmin) = step_fct(Tmin-T)*P
Pr(P, T, Tmin) = step_fct(T-Tmin)*P
M(S0, T, Df, Tmax) = step_fct(T-Tmax)*step_fct(S0)*minimum([S0, Df*(T-Tmax)])
PET(T, Lday) = 29.8 * Lday * 0.611 * exp((17.3*T)/(T+237.3)) / (T + 273.2)
ET(S1, T, Lday, Smax) = step_fct(S1)*step_fct(S1-Smax)*PET(T,Lday) + step_fct(S1)*step_fct(Smax-S1)*PET(T,Lday)*(S1/Smax)
Qb(S1,f,Smax,Qmax) = step_fct(S1)*step_fct(S1-Smax)*Qmax + step_fct(S1)*step_fct(Smax-S1)*Qmax*exp(-f*(Smax-S1))
Qs(S1, Smax) = step_fct(S1)*step_fct(S1-Smax)*(S1-Smax)

function basic_bucket_incl_states(p_, itp_Lday, itp_P, itp_T, t_out)
    function exp_hydro_optim_states!(dS,S,ps,t)
        f, Smax, Qmax, Df, Tmax, Tmin = ps
        Lday = itp_Lday(t)
        P    = itp_P(t)
        T    = itp_T(t)
        Q_out = Qb(S[2],f,Smax,Qmax) + Qs(S[2], Smax)
        dS[1] = Ps(P, T, Tmin) - M(S[1], T, Df, Tmax)
        dS[2] = Pr(P, T, Tmin) + M(S[1], T, Df, Tmax) - ET(S[2], T, Lday, Smax) - Q_out
    end

    prob = ODEProblem(exp_hydro_optim_states!, p_[1:2], Float64.((t_out[1], maximum(t_out))))
    # sol = solve(prob, BS3(), u0=p_[1:2], p=p_[3:end], saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=ForwardDiffSensitivity())
    sol = solve(prob, BS3(), u0=p_[1:2], p=p_[3:end], saveat=t_out, reltol=1e-3, abstol=1e-3)
    Qb_ = Qb.(sol[2,:], p_[3], p_[4], p_[5])
    Qs_ = Qs.(sol[2,:], p_[4])
    Qout_ = Qb_.+Qs_
    return Qout_, sol
end

function NeuralODE_M100(p, norm_S0, norm_S1, norm_P, norm_T, itp_Lday, itp_P, itp_T, t_out, ann; S_init = [0.0, 0.0])
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