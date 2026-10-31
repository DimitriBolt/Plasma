"""Temperature extrema, hot intervals, energy and applicability diagnostics."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from .bubble_thermo import (Parameters,Gas,volume,adiabatic_temperature,adiabatic_pressure,
    gas_pressure,gas_density)
from .heat_transfer import heat_loss,relaxation_time
from .numerics import thermal_extrema,hot_intervals,roots_from_grid
from .one_way import temperature,temperature_derivative

THRESHOLDS=(1e3,1e4,1e5,1e6)


def threshold_metrics(Tfun,Rfun,grid):
    rows=[];gx,gw=leggauss(20)
    for threshold in THRESHOLDS:
        intervals=hot_intervals(Tfun,grid,threshold)
        exposure=0.
        for a,b in intervals:
            # Integrate per grid subinterval, not across unresolved internal structure.
            edges=np.unique(np.r_[a,grid[(grid>a)&(grid<b)],b])
            mid=(edges[1:]+edges[:-1])/2;half=np.diff(edges)/2
            nodes=mid[:,None]+half[:,None]*gx
            radii=Rfun(nodes.ravel()).reshape(nodes.shape)
            exposure+=np.sum(volume(radii)*gw*half[:,None])
        rows.append(dict(threshold_K=threshold,duration_s=sum(b-a for a,b in intervals),
            hot_volume_time_m3_s=float(exposure),intervals=intervals))
    return rows


def adiabatic_metrics(path,gas=Gas(),details=False):
    p=Parameters();s=path.summary
    tf=lambda t:adiabatic_temperature(path(t)[0],gas,p)
    rt=path.radius_extrema_times();vals=tf(rt);j=np.argmax(vals)
    # Independent result uses AUTHORITATIVE refined Rmin, not the interpolated argmin.
    independent=float(adiabatic_temperature(s['Rmin_m'],gas,p))
    tpeak=float(rt[j]);Tpeak=float(vals[j]);Rpeak=float(path(tpeak)[0])
    common=dict(mechanics=path.model,thermodynamics='A_adiabatic',Pa_bar=path.amplitude,
        gas=gas.name,gamma=gas.gamma,T0_K=gas.T0,kg_W_m_K=gas.kg,eta=0.,
        completed=path.completed,extrema_scope=s['extrema_scope'],last_time_s=path.end,
        observed_Tmax_K=Tpeak,global_Tmax_K=Tpeak if path.completed else np.nan,
        t_Tmax_s=tpeak,peak_cycle=int(tpeak*p.f)+1,R_at_Tmax_m=Rpeak,
        Tmax_from_refined_Rmin_K=independent,refined_vs_interpolant_relative_error=abs(Tpeak/independent-1),
        raw_saved_peak_relative_deficit=1-float(adiabatic_temperature(path.export.R_m.min(),gas,p))/independent,
        pg_at_Tmax_Pa=float(gas_pressure(Rpeak,Tpeak,gas,p)),V_at_Tmax_m3=float(volume(Rpeak)),
        density_at_Tmax_kg_m3=float(gas_density(Rpeak,gas,p)),
        max_wall_Mach=s['max_wall_Mach'],pressure_bulk_modulus_proxy=s['max_liquid_pressure_change_over_bulk_modulus'],
        numerical_status='inherited_saved_mechanics_plus_interpolation_check',
        applicability=('formal_RP_diagnostic' if path.model=='RP' and s['max_wall_Mach']>1 else
                       'formal_extreme_scenario' if path.amplitude>=1.2 else 'conditional_model_based_estimate'))
    if details:
        grid=np.unique(np.r_[path.t,rt])
        hot=threshold_metrics(tf,lambda t:path(t)[0],grid)
        for row in hot:common[f"duration_above_{int(row['threshold_K'])}K_s"]=row['duration_s']
        return common,hot
    return common


def solution_functions(sol):
    if sol.kind=='fully_coupled':
        Rfun=lambda t:sol.sol(t)[0]
        Ufun=lambda t:sol.sol(t)[1]
        Tfun=lambda t:sol.sol(t)[2]
        dTfun=lambda t:np.asarray([sol.fun(x,sol.sol(x))[2] for x in np.atleast_1d(t)]).reshape(np.shape(t))
    else:
        Rfun=lambda t:sol.path(t)[0]
        Ufun=lambda t:sol.path(t)[1]
        Tfun=lambda t:temperature(sol,t)
        dTfun=lambda t:temperature_derivative(sol,t)
    return Rfun,Ufun,Tfun,dTfun


def solution_metrics(sol,details=True):
    gas,p=sol.gas,sol.p
    Rfun,Ufun,Tfun,dTfun=solution_functions(sol)
    grid=np.unique(np.r_[sol.t,(sol.t[1:]+sol.t[:-1])/2])
    if sol.kind!='fully_coupled':grid=np.unique(np.r_[grid,sol.path.t[sol.path.t<=sol.t[-1]]])
    Tpeak,tpeak,tex=thermal_extrema(Tfun,dTfun,grid)
    rroots=roots_from_grid(Ufun,grid)
    grid=np.unique(np.r_[grid,tex,rroots])
    T=Tfun(grid);R=Rfun(grid);U=Ufun(grid);pg=gas_pressure(R,T,gas,p)
    y=sol.sol(grid); W,Q,S=(y[3],y[4],y[5]) if sol.kind=='fully_coupled' else (y[1],y[2],y[3])
    residual=T/gas.T0-1+W+Q-S
    scale=max(1,float(np.max(np.abs(T/gas.T0-1))),float(np.max(abs(W))),float(np.max(abs(Q))))
    pL=pg+p.pv-2*p.sigma/R-4*p.mu*U/R
    Rpeak=float(Rfun(tpeak))
    # Refine reference minimum from velocity roots and endpoints.
    rcandidates=np.unique(np.r_[0,sol.t[-1],rroots]);rmin_time=float(rcandidates[np.argmin(Rfun(rcandidates))])
    tau=relaxation_time(Rpeak,gas,sol.eta,p)
    mechanical_full=sol.kind=='fully_coupled' or sol.path.completed
    completed=sol.completed and mechanical_full
    mech='KM' if sol.kind=='fully_coupled' else sol.path.model
    m=dict(mechanics=mech,thermodynamics='B_'+sol.kind,Pa_bar=sol.amplitude,
        gas=gas.name,gamma=gas.gamma,T0_K=gas.T0,kg_W_m_K=gas.kg,eta=sol.eta,
        completed=completed,extrema_scope='global_on_requested_interval' if completed else 'observed_interval_censored',
        last_time_s=sol.t[-1],observed_Tmax_K=Tpeak,global_Tmax_K=Tpeak if completed else np.nan,
        t_Tmax_s=tpeak,peak_cycle=int(tpeak*p.f)+1,R_at_Tmax_m=Rpeak,
        t_Rmin_s=rmin_time,Tmax_minus_Rmin_time_s=tpeak-rmin_time,
        pg_at_Tmax_Pa=float(gas_pressure(Rpeak,Tpeak,gas,p)),max_pg_Pa=float(pg.max()),
        V_at_Tmax_m3=float(volume(Rpeak)),density_at_Tmax_kg_m3=float(gas_density(Rpeak,gas,p)),
        max_wall_Mach=float(np.max(abs(U))/p.c),
        pressure_bulk_modulus_proxy=float(np.max(abs(pL-p.p0))/(p.rho*p.c**2)),
        tau_thermal_at_peak_s=tau,tau_thermal_at_R0_s=relaxation_time(p.R0,gas,sol.eta,p),
        signed_heat_loss_J=float(Q[-1]*gas.E0(p)),work_by_gas_J=float(W[-1]*gas.E0(p)),
        source_energy_J=float(S[-1]*gas.E0(p)),energy_residual_max_J=float(np.max(abs(residual))*gas.E0(p)),
        energy_residual_scaled=float(np.max(abs(residual))/scale),
        solver_status=sol.status,solver_message=sol.message,accepted_steps=len(sol.t)-1,nfev=sol.nfev,
        numerical_status='requires_saved_convergence_comparison',
        applicability='formal_extreme_scenario' if sol.amplitude>=1.2 else 'conditional_model_based_estimate')
    if sol.kind!='fully_coupled':
        oldpg=p.pg0*(p.R0/R)**(3*p.kappa)
        m['max_abs_pressure_relative_mismatch']=float(np.max(abs(pg/oldpg-1)))
        m['pressure_mismatch_at_Tmax']=float(gas_pressure(Rpeak,Tpeak,gas,p)/(p.pg0*(p.R0/Rpeak)**(3*p.kappa))-1)
    hot=threshold_metrics(Tfun,Rfun,grid) if details else []
    for row in hot:m[f"duration_above_{int(row['threshold_K'])}K_s"]=row['duration_s']
    return m,hot,grid


def save_series(file,sol,grid):
    Rfun,Ufun,Tfun,_=solution_functions(sol);R,U,T=Rfun(grid),Ufun(grid),Tfun(grid)
    y=sol.sol(grid);W,Q,S=(y[3],y[4],y[5]) if sol.kind=='fully_coupled' else (y[1],y[2],y[3])
    np.savez_compressed(file,t_s=grid,R_m=R,U_m_s=U,T_K=T,pg_Pa=gas_pressure(R,T,sol.gas,sol.p),
        V_m3=volume(R),density_kg_m3=gas_density(R,sol.gas,sol.p),
        Qloss_W=heat_loss(R,T,sol.gas,sol.eta),work_J=W*sol.gas.E0(sol.p),
        signed_heat_loss_J=Q*sol.gas.E0(sol.p),source_energy_J=S*sol.gas.E0(sol.p))
