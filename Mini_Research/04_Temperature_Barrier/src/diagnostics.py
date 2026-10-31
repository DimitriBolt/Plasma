"""Refined thermal extrema and signed energy ledgers on observed intervals."""
import math
import numpy as np
from scipy.optimize import brentq
from .upstream import Parameters
from .thermodynamics import T0, RU,volume

def peak_conditioning(series,relative_tolerance=1e-8):
    """Conditioning of argmax time, without changing the selected maximum.

    The declared near-tie band is a diagnostic convention, not an error bound.
    All candidate times and actual values remain in the physical series.
    """
    times=np.unique(np.r_[series['t_s'][0],series['t_s'][-1],series['temperature_stationary_times_s']])
    ids=np.searchsorted(series['t_s'],times);values=series['T_K'][ids]
    peak=float(np.max(series['T_K']));band=peak*relative_tolerance
    near=times[values>=peak-band]
    return dict(near_peak_relative_tolerance=relative_tolerance,near_peak_tolerance_K=band,
                near_peak_candidate_count=len(near),near_peak_time_first_s=float(near.min()),
                near_peak_time_last_s=float(near.max()))

def roots(fun,grid,vals=None):
    if vals is None: vals=np.array([fun(t) for t in grid])
    ids=np.flatnonzero(vals[:-1]*vals[1:]<0)
    rr=[brentq(fun,grid[i],grid[i+1],xtol=1e-24,rtol=4*np.finfo(float).eps) for i in ids]
    return np.unique(np.r_[rr,grid[vals==0]])

def evaluate(sol):
    mix=sol.mix;p=Parameters();coupled=sol.mode=='coupled'
    state=lambda t:sol.state(t,sol.sol(t))
    temp=lambda t:float(state(t)[2])
    speed=lambda t:float(state(t)[1])
    dT=lambda t:mix.temperature_rate(*state(t),sol.eta,trial=True)
    grid=np.unique(np.r_[sol.t,*[sol.t[:-1]+w*np.diff(sol.t) for w in (.25,.5,.75)]])
    if not coupled:
        grid=np.unique(np.r_[grid,sol.path.t[sol.path.t<=sol.t[-1]]])
    R,U,T=state(grid);pg=mix.nt*RU*T/volume(R)
    rates=(-pg*4*np.pi*R**2*U-4*np.pi*sol.eta*mix.kg*R*(T-T0))/(mix.nt*mix.cv_array(T))
    tex=roots(dT,grid,rates);rex=roots(speed,grid,U)
    candidates=np.unique(np.r_[0.,sol.t[-1],tex]);Ts=np.asarray(state(candidates)[2])
    j=Ts.argmax();tp=float(candidates[j]);peak=float(Ts[j])
    grid=np.unique(np.r_[grid,tex,rex,tp])
    R,U,T=state(grid);V=volume(R);pg=mix.nt*RU*T/V
    du=mix.du_array(T);yy=sol.sol(grid);ix=3 if coupled else 1
    W,Q,C=yy[ix:ix+3]*mix.Eref
    resid=du+W+Q
    scale=max(mix.Eref,float(np.max(abs(du))),float(np.max(abs(W))),float(np.max(abs(Q))))
    ip=int(np.flatnonzero(grid==tp)[0]);Rp=float(R[ip]);Tp=float(T[ip]);pp=float(pg[ip])
    rc=np.unique(np.r_[0.,sol.t[-1],rex]);rv=np.asarray(state(rc)[0]);tr=float(rc[rv.argmin()])
    rmax=[t for t in rex if 0<t<tp and speed(max(0,t-1e-13))>0 and speed(min(sol.t[-1],t+1e-13))<0]
    tc=float(max(rmax)) if rmax else 0.
    yc=sol.sol(tc);yp=sol.sol(tp);Uc=mix.du(temp(tc));Cp=float((yp[ix+2]-yc[ix+2])*mix.Eref)
    pL=pg+mix.p_constant-2*p.sigma/R-4*p.mu*U/R
    pressure_old=p.pg0*(p.R0/R)**(3*p.kappa)+p.pv
    mismatch=(pg+mix.p_constant)/pressure_old-1
    # Temperature half-height duration around the selected peak (relative to T0).
    half=T0+(peak-T0)/2;cross=roots(lambda t:temp(t)-half,grid,T-half)
    left=cross[cross<tp];right=cross[cross>tp]
    width=float(right[0]-left[-1]) if len(left) and len(right) else math.nan
    complete=sol.completed and (coupled or sol.path.completed)
    summary=dict(Pa_bar=sol.amplitude,gas=mix.gas,mode=sol.mode,cv_mode=mix.cv_mode,
        vapor_factor=mix.vapor_factor,vapor_mole_fraction=mix.fractions.get('H2O',0.),
        vapor_convention=mix.convention,n_total_mol=mix.nt,n_H2O_mol=mix.nt*mix.fractions.get('H2O',0.),
        kg_W_m_K=mix.kg,eta=sol.eta,Tmin_domain_K=mix.Tmin,Tmax_domain_K=mix.Tmax,
        completed=complete,extrema_scope='five_cycle_global' if complete else 'observed_interval_censored',
        observed_Tmax_K=peak,global_Tmax_K=peak if complete else math.nan,
        start_s=0.,end_s=float(sol.t[-1]),boundary_temperature_K=float(T[-1]),
        t_Tmax_s=tp,peak_cycle=min(5,int(tp*p.f)+1),R_at_Tmax_m=Rp,p_at_Tmax_Pa=pp,
        p_mechanical_at_Tmax_Pa=pp+mix.p_constant,density_at_Tmax_kg_m3=mix.nt*mix.mass_mol/float(volume(Rp)),
        number_density_at_Tmax_m3=mix.nt*6.02214076e23/float(volume(Rp)),
        max_p_Pa=float(pg.max()),max_density_kg_m3=float((mix.nt*mix.mass_mol/V).max()),
        max_wall_Mach=float(np.max(abs(U))/p.c),pressure_bulk_modulus_proxy=float(np.max(abs(pL-p.p0))/(p.rho*p.c**2)),
        t_Rmin_s=tr,Tmax_minus_Rmin_s=tp-tr,Rmin_m=float(rv.min()),
        gamma_at_peak=1+RU/mix.cv(Tp),thermal_acoustic_crossing_peak_s=Rp/math.sqrt((1+RU/mix.cv(Tp))*RU*Tp/mix.mass_mol),
        temperature_half_height_width_s=width,collapse_start_s=tc,collapse_to_Tmax_s=tp-tc,
        delta_U_peak_J=float(du[ip]),work_by_gas_peak_J=float(W[ip]),signed_heat_loss_peak_J=float(Q[ip]),
        gross_compression_peak_J=float(C[ip]),collapse_compression_J=Cp,collapse_delta_U_J=mix.du(Tp)-Uc,
        delta_U_end_J=float(du[-1]),work_by_gas_end_J=float(W[-1]),signed_heat_loss_end_J=float(Q[-1]),
        gross_compression_end_J=float(C[-1]),chemical_energy_change_J=0.,source_energy_J=0.,
        energy_residual_max_J=float(np.max(abs(resid))),energy_residual_scaled=float(np.max(abs(resid))/scale),
        solver_status=sol.status,solver_message=sol.message,accepted_steps=len(sol.t)-1,nfev=sol.nfev,
        max_abs_pressure_relative_mismatch=float(np.max(abs(mismatch))) if not coupled else math.nan,
        pressure_relative_mismatch_at_peak=float(mismatch[ip]) if not coupled else math.nan,
        applicability='formal_extreme_scenario' if sol.amplitude==1.3 else 'conditional_model_based_estimate')
    series=dict(t_s=grid,R_m=R,U_m_s=U,T_K=T,p_Pa=pg,delta_U_J=du,work_by_gas_J=W,
        signed_heat_loss_J=Q,gross_compression_J=C,energy_residual_J=resid,
        pressure_relative_mismatch=mismatch if not coupled else np.full_like(T,np.nan),
        density_kg_m3=mix.nt*mix.mass_mol/V,temperature_stationary_times_s=tex)
    summary.update(peak_conditioning(series))
    return summary,series
