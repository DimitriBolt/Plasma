"""Bounded numerical verification; never regenerates MR1 or MR2 outputs."""
from dataclasses import asdict
import json
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from .bubble_thermo import ROOT,Parameters,Gas,adiabatic_temperature,gas_pressure,mechanics
from .saved_mechanics import SavedMechanics
from .coupled_model import solve_coupled
from .one_way import solve_one_way
from .diagnostics import solution_metrics,solution_functions,adiabatic_metrics
from .heat_transfer import heat_loss
from .numerics import Settings
from .run_experiments import json_write,AMPLITUDES

OUT=ROOT/'results'


def independent_energy(sol):
    """Gauss integration on each accepted interval, independent of ODE ledgers."""
    Rf,Uf,Tf,_=solution_functions(sol);g,p=sol.gas,sol.p
    gx,gw=leggauss(16)
    mid=(sol.t[1:]+sol.t[:-1])/2;half=np.diff(sol.t)/2
    t=mid[:,None]+half[:,None]*gx
    R,U,T=Rf(t.ravel()).reshape(t.shape),Uf(t.ravel()).reshape(t.shape),Tf(t.ravel()).reshape(t.shape)
    work=gas_pressure(R,T,g,p)*4*np.pi*R**2*U
    loss=heat_loss(R,T,g,sol.eta)
    W=np.r_[0,np.cumsum(np.sum(work*gw*half[:,None],axis=1))]
    Q=np.r_[0,np.cumsum(np.sum(loss*gw*half[:,None],axis=1))]
    residual=g.capacity(p)*(Tf(sol.t)-g.T0)+W+Q
    scale=max(g.E0(p),np.max(abs(W)),np.max(abs(Q)),np.max(abs(g.capacity(p)*(Tf(sol.t)-g.T0))))
    return dict(gauss_energy_max_residual_J=float(np.max(abs(residual))),
        gauss_energy_scaled_residual=float(np.max(abs(residual))/scale),
        outward_heat_J=float(np.sum(np.maximum(loss,0)*gw*half[:,None])),
        inward_heat_J=float(-np.sum(np.minimum(loss,0)*gw*half[:,None])))


def differences(m,ref):
    d=dict(relative_Tmax_error=abs(m['observed_Tmax_K']/ref['observed_Tmax_K']-1),
           abs_peak_time_error_s=abs(m['t_Tmax_s']-ref['t_Tmax_s']))
    for q in (1000,10000,100000,1000000):
        name=f'duration_above_{q}K_s';a=m.get(name,np.nan);b=ref.get(name,np.nan)
        d[f'abs_duration_{q}K_error_s']=abs(a-b)
        d[f'relative_duration_{q}K_error']=abs(a/b-1) if b else (0. if a==0 else np.nan)
    return d


def numerical_convergence():
    production=pd.read_csv(OUT/'thermal_summary.csv');rows=[];recovery=[]
    for amplitude in AMPLITUDES:
        for eta in (0.,1.):
            prod=production[(production.Pa_bar==amplitude)&(production.gas=='air_like')&
                (production.eta==eta)&(production.thermodynamics=='B_fully_coupled')].iloc[0].to_dict()
            ref=solve_coupled(amplitude,eta=eta,settings=Settings('DOP853',2e-11,1200))
            rm,_,_=solution_metrics(ref);energy=independent_energy(ref)
            configurations=[('production',prod,Settings(),None)]
            if eta==1:
                coarse=solve_coupled(amplitude,eta=eta,settings=Settings('DOP853',1e-7,300))
                cm,_,_=solution_metrics(coarse);configurations.append(('coarse',cm,coarse.settings,coarse))
            other=solve_coupled(amplitude,eta=eta,settings=Settings('Radau',2e-9,600))
            om,_,_=solution_metrics(other);configurations.append(('independent_Radau',om,other.settings,other))
            configurations.append(('tight_reference',rm,ref.settings,ref))
            for label,m,settings,sol in configurations:
                row=dict(configuration=label,**asdict(settings),**m,**differences(m,rm))
                if sol is not None:row.update(independent_energy(sol))
                rows.append(row)
                print('convergence',amplitude,eta,label,row['relative_Tmax_error'],m['completed'],flush=True)
            if eta==0:
                path=SavedMechanics('KM',amplitude)
                Rf,Uf,Tf,_=solution_functions(ref)
                rt=np.unique(np.r_[ref.t,path.t]);rt=rt[rt<=ref.t[-1]]
                temp_error=float(np.max(abs(Tf(rt)/adiabatic_temperature(Rf(rt))-1)))
                # Compare the coupled run to the ACCEPTED trajectories, no new mechanical solve.
                R0,U0=path(rt)
                recovery.append(dict(Pa_bar=amplitude,completed=ref.completed,
                    max_temperature_invariant_relative_error=temp_error,
                    max_common_time_radius_error_over_R0=float(np.max(abs(Rf(rt)-R0))/Parameters().R0),
                    max_common_time_velocity_error_m_s=float(np.max(abs(Uf(rt)-U0))),
                    relative_Tmax_vs_saved_mechanics=abs(rm['observed_Tmax_K']/adiabatic_metrics(path)['observed_Tmax_K']-1),
                    **energy))
            pd.DataFrame(rows).to_csv(OUT/'solver_convergence.csv',index=False)
            pd.DataFrame(recovery).to_csv(OUT/'adiabatic_recovery.csv',index=False)
    # Every transport/composition sensitivity receives a tolerance/step check.
    for amplitude in AMPLITUDES:
        for gas,eta in ((Gas(),.1),(Gas(),10.),(Gas(gamma=5/3,kg=.0177,molar_mass=.039948,name='argon_like'),1.)):
            ref=solve_coupled(amplitude,gas,eta,settings=Settings('DOP853',2e-11,1200));rm,_,_=solution_metrics(ref)
            pm=production[(production.Pa_bar==amplitude)&(production.gas==gas.name)&(production.eta==eta)&
                         (production.thermodynamics=='B_fully_coupled')].iloc[0].to_dict()
            for label,m,sol in [('production',pm,None),('tight_reference',rm,ref)]:
                row=dict(**m,configuration=label,**differences(m,rm))
                if sol is not None:row.update(independent_energy(sol))
                rows.append(row)
            print('sensitivity tolerance',amplitude,gas.name,eta,differences(pm,rm)['relative_Tmax_error'],flush=True)
            pd.DataFrame(rows).to_csv(OUT/'solver_convergence.csv',index=False)


def input_resolution():
    production=pd.read_csv(OUT/'thermal_summary.csv');rows=[];interpolation=[]
    for amplitude in AMPLITUDES:
        full=SavedMechanics('KM',amplitude)
        prod=production[(production.Pa_bar==amplitude)&(production.gas=='air_like')&
                       (production.eta==1)&(production.thermodynamics=='B_one_way_post_processing')].iloc[0].to_dict()
        ref=solve_one_way(full,eta=1,settings=Settings('DOP853',2e-11,1200));rm,_,_=solution_metrics(ref)
        rows.append(dict(**prod,grid='full_export_production',**differences(prod,rm)))
        rows.append(dict(**rm,grid='full_export_tight',**differences(rm,rm),**independent_energy(ref)))
        for gas,eta in ((Gas(),0.),(Gas(),.1),(Gas(),10.),(Gas(gamma=5/3,kg=.0177,molar_mass=.039948,name='argon_like'),1.)):
            sol=solve_one_way(full,gas,eta,settings=Settings('DOP853',2e-11,1200));m,_,_=solution_metrics(sol)
            pm=production[(production.Pa_bar==amplitude)&(production.gas==gas.name)&(production.eta==eta)&
                         (production.thermodynamics=='B_one_way_post_processing')].iloc[0].to_dict()
            rows.append(dict(**pm,grid='full_export_production',**differences(pm,m)))
            rows.append(dict(**m,grid='full_export_tight',**differences(m,m),**independent_energy(sol)))
        for grid,stride in (('export',2),('export',4),('accepted',1)):
            path=SavedMechanics('KM',amplitude,stride,grid)
            sol=solve_one_way(path,eta=1,settings=Settings('DOP853',2e-11,1200))
            m,_,_=solution_metrics(sol)
            rows.append(dict(**m,grid=f'{grid}_stride{stride}',**differences(m,rm),**independent_energy(sol)))
            # Withheld export locations quantify R/U interpolation (U is derivative of R).
            t=full.export.t_s.to_numpy();R,U=path(t)
            interpolation.append(dict(Pa_bar=amplitude,grid=f'{grid}_stride{stride}',nodes=len(path.t),
                max_R_error_over_R0=float(np.max(abs(R-full.export.R_m))/Parameters().R0),
                max_U_error_m_s=float(np.max(abs(U-full.export.U_m_s))),
                relative_Rmin_temperature_error=adiabatic_metrics(path)['refined_vs_interpolant_relative_error']))
            print('input grid',amplitude,grid,stride,rows[-1]['relative_Tmax_error'],flush=True)
        # Adiabatic hot widths: extrema are retained at each resolution.
        for model in ('KM','RP'):
            fullA=SavedMechanics(model,amplitude);am,hot=adiabatic_metrics(fullA,details=True)
            for grid,stride in (('export',2),('export',4),('accepted',1)):
                path=SavedMechanics(model,amplitude,stride,grid)
                m,_=adiabatic_metrics(path,details=True)
                rows.append(dict(**m,grid=f'{grid}_stride{stride}',**differences(m,am)))
        pd.DataFrame(rows).to_csv(OUT/'input_grid_convergence.csv',index=False)
        pd.DataFrame(interpolation).to_csv(OUT/'interpolation_checks.csv',index=False)


def finalize_numerical_status():
    b=pd.read_csv(OUT/'thermal_summary.csv',float_precision='round_trip')
    c=pd.read_csv(OUT/'solver_convergence.csv');q=pd.read_csv(OUT/'input_grid_convergence.csv')
    assert len(c)==39 and c.completed.all() and len(q)==57
    for i,r in b.iterrows():
        independent=r.gas=='air_like' and r.eta in (0.,1.) and r.thermodynamics=='B_fully_coupled'
        status='tolerance_step_and_independent_Radau_checked' if independent else 'tolerance_and_step_checked'
        table=c if r.thermodynamics=='B_fully_coupled' else q
        mask=(table.Pa_bar==r.Pa_bar)&(table.gas==r.gas)&(table.eta==r.eta)
        if r.thermodynamics=='B_fully_coupled':mask &= table.configuration=='production'
        else:mask &= (table.grid=='full_export_production')&(table.thermodynamics==r.thermodynamics)
        error=float(table[mask].iloc[0].relative_Tmax_error)
        b.loc[i,'numerical_status']=status;b.loc[i,'numerical_peak_difference_relative']=error
        file=OUT/(r.run_id+'.json');meta=json.loads(file.read_text())
        meta['summary']['numerical_status']=status
        meta['summary']['numerical_peak_difference_relative']=error
        json_write(file,meta)
    b.to_csv(OUT/'thermal_summary.csv',index=False)


def main():
    numerical_convergence();input_resolution();finalize_numerical_status()


if __name__=='__main__':main()
