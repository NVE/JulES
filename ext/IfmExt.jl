"""
This code is copied from https://github.com/marv-in/HydroNODE
and has BSD 3-Clause License
"""
# TODO: Update License file
# TODO: Write license info in each source file

module IfmExt

using OrdinaryDiffEq
using JulES

step_fct(x) = (tanh(5.0*x) + 1.0)*0.5
Ps(P, T, Tmin) = step_fct(Tmin-T)*P
Pr(P, T, Tmin) = step_fct(T-Tmin)*P
M(S0, T, Df, Tmax) = step_fct(T-Tmax)*step_fct(S0)*minimum([S0, Df*(T-Tmax)])
PET(T, Lday) = 29.8 * Lday * 0.611 * exp((17.3*T)/(T+237.3)) / (T + 273.2)
ET(S1, T, Lday, Smax) = step_fct(S1)*step_fct(S1-Smax)*PET(T,Lday) + step_fct(S1)*step_fct(Smax-S1)*PET(T,Lday)*(S1/Smax)
Qb(S1,f,Smax,Qmax) = step_fct(S1)*step_fct(S1-Smax)*Qmax + step_fct(S1)*step_fct(Smax-S1)*Qmax*exp(-f*(Smax-S1))
Qs(S1, Smax) = step_fct(S1)*step_fct(S1-Smax)*(S1-Smax)

function JulES.basic_bucket_incl_states(p_, itp_Lday, itp_P, itp_T, t_out)
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

end