"""Reproduce the bounded new calculations; read accepted results without writes."""
import hashlib,json,platform,sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq
from scipy.integrate import quad
from .upstream import ROOT,REPO,thermo
from .mechanics import solve,Tait,mr2

OUT=ROOT/'results'

def inventory():
    old=json.loads((OUT/'accepted_tree_sha256.json').read_text())
    changed=[f for f,h in old.items() if not (REPO/f).exists() or hashlib.sha256((REPO/f).read_bytes()).hexdigest()!=h]
    if changed:raise AssertionError('Accepted input changes: '+str(changed))
    return {'accepted_files_checked':len(old),'changed':changed}

def mechanics():
    rows=[]
    for a in [1.,1.15,1.30]:
        for name in ['KM','Gilmore']:
            for method,rtol,steps,label in [('DOP853',1e-9,800,'production'),('DOP853',1e-11,1600,'strict'),('Radau',2e-10,1200,'independent')]:
                row,sol,t=solve(a,name,method,rtol,steps);row['control']=label;rows.append(row)
                if label=='strict':
                    np.savez_compressed(OUT/f'mechanics_{name}_{a:.2f}.npz',t_s=t,R_m=sol.sol(t)[0],U_m_s=sol.sol(t)[1])
                print('mechanics',a,name,label,row['Rmin_m'],flush=True)
    df=pd.DataFrame(rows)
    controls=[]
    for (a,model),g in df.groupby(['Pa_bar','model']):
        ref=g[g.control=='strict'].iloc[0]
        for _,r in g.iterrows():
            controls.append(dict(Pa_bar=a,model=model,control=r.control,
                relative_Rmin=abs(r.Rmin_m/ref.Rmin_m-1),relative_speed=abs(r.max_speed_m_s/ref.max_speed_m_s-1),
                event_time_difference_s=r.t_Rmin_s-ref.t_Rmin_s))
    df.to_csv(OUT/'mechanics.csv',index=False,float_format='%.17g')
    pd.DataFrame(controls).to_csv(OUT/'convergence.csv',index=False,float_format='%.17g')

def evidence():
    mr7=ROOT.parent/'07_Integrated_Model/results'
    summary=pd.read_csv(mr7/'summary.csv',float_precision='round_trip').set_index('run_id')
    diag=pd.read_csv(mr7/'physical_diagnostics.csv',float_precision='round_trip')
    diag.to_csv(OUT/'inherited_physical_diagnostics.csv',index=False,float_format='%.17g')
    ledger=[];ideal=[];pairs=[];sources=[];kinematics=[]
    for name in ['base_1.00','base_1.15','base_1.30','stable_base_1.00_F10_L3_q3','stable_base_1.15_F10_L3_q3','stable_extra_high_1.30_BDF_base_F10_L3_q3']:
        path=mr7/'references'/f'{name}.npz';z=np.load(path)
        sources.append(dict(path=str(path.relative_to(REPO)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        ug,ew,ea,w,qg,qa,qb,ab=[z[k] for k in ['U_gas_change_J','E_water_J','E_Au_J','work_by_gas_J','Q_l_to_g_J','Q_Au_to_l_J','Q_bath_J','E_absorbed_J']]
        residuals=[ug+w-qg,ew-qa+qg+qb,ea-ab+qa,ug+ew+ea+w+qb-ab]
        scale=max(np.max(abs(ug)),np.max(abs(w)),np.max(abs(qg)),1e-30)
        row=dict(run_id=name,end_s=float(z['t_s'][-1]),observed_Tmax_K=float(np.max(z['T_g_K'])),
           E_incident_J=float(z['E_incident_J'][-1]),E_absorbed_J=float(ab[-1]),
           E_Au_J=float(ea[-1]),E_liquid_J=float(ew[-1]),delta_U_gas_J=float(ug[-1]),
           W_by_gas_J=float(w[-1]),Q_l_to_g_J=float(qg[-1]),Q_Au_to_l_J=float(qa[-1]),Q_bath_J=float(qb[-1]))
        for k,v in zip(['gas','liquid','Au','total'],residuals):row[k+'_max_residual_J']=float(np.max(abs(v)));row[k+'_scaled_residual']=float(np.max(abs(v))/scale)
        row['signed_gas_power_min_W']=float(np.min(z['Q_l_to_g_W']));row['signed_gas_power_max_W']=float(np.max(z['Q_l_to_g_W']))
        ledger.append(row)
        if name.startswith('base_'):
            # Kinematic mismatch, not a force law: required constant-gap speed
            # minus the incompressible passive-water velocity at the center.
            R=z['R_m'];U=z['U_m_s'];rnp=R+120e-9
            slip=U*(1-(R/rnp)**2)
            kinematics.append(dict(run_id=name,max_required_slip_m_s=float(np.max(abs(slip))),
                max_particle_Re=float(np.max(2*20e-9*998.2*abs(slip)/.001002)),
                Stokes_relaxation_reference_s=2*19320*(20e-9)**2/(9*.001002),
                meaning='Required slip for prescribed confinement relative to accepted incompressible flow; no retaining force derived.'))
        if 'paired_delta_T_K' in z:
            # Independent nonlinear sensible-energy increment, same common prefix.
            mix=thermo.Mixture();tb=z['T_g_K']-z['paired_delta_T_K']
            delta=np.array([d*quad(lambda s:mix.capacity(float(b+s*d)),0,1,epsabs=1e-30)[0]
                       for b,d in zip(tb,z['paired_delta_T_K'])])
            dl=z['paired_delta_ledgers_J'];r=delta+dl[0]-dl[1]
            pairs.append(dict(run_id=name,E_absorbed_J=float(dl[4,-1]),
                delta_U_gas_end_J=float(delta[-1]),delta_work_by_gas_end_J=float(dl[0,-1]),
                delta_Q_l_to_g_end_J=float(dl[1,-1]),delta_Q_Au_to_l_end_J=float(dl[2,-1]),
                delta_Q_bath_end_J=float(dl[3,-1]),
                independent_gas_increment_residual_over_absorption=float(np.max(abs(r))/dl[4,-1]),
                stored_paired_total_residual_over_absorption=float(np.max(abs(z['paired_total_residual_J']))/dl[4,-1]),
                interpretation='Same-prefix nonlinear increment; net gas flux includes mechanical feedback, not tagged photon delivery.'))
    pd.DataFrame(ledger).to_csv(OUT/'energy_audit.csv',index=False,float_format='%.17g')
    pd.DataFrame(pairs).to_csv(OUT/'paired_energy_audit.csv',index=False,float_format='%.17g')
    pd.DataFrame(kinematics).to_csv(OUT/'particle_kinematics.csv',index=False,float_format='%.17g')
    for gas,vapor,temps in [('air',1,[553.0230367271835,2798.2812622841393]),('air',0,[6263.578953999877]),('argon',0,[22808.715597374918])]:
        mix=thermo.Mixture(gas=gas,vapor_factor=vapor)
        for T in temps:
            for F in [1.,10.,100.,1000.]:
                E=3.4337324474065e-15*F
                umax=mix.du(mix.Tmax) if np.isfinite(mix.Tmax) else np.inf
                target=mix.du(T)+E
                top=mix.Tmax if np.isfinite(mix.Tmax) else T+E/mix.capacity(T)*2+1
                Tf=brentq(lambda x:mix.du(x)-target,T,top,xtol=1e-10) if target<umax else np.nan
                ideal.append(dict(gas=gas,vapor_factor=vapor,baseline_T_K=T,fluence_J_m2=F,
                    E_absorbed_J=E,ideal_all_absorbed_at_fixed_volume_delta_T_K=Tf-T,
                    meaning='Defined fixed-volume, frozen-composition caloric allowance only; not a dynamic or universal upper bound.'))
    pd.DataFrame(ideal).to_csv(OUT/'ideal_caloric_allowance.csv',index=False,float_format='%.17g')
    p=mr2.Parameters();alpha=.6/(1000*4187);a=20e-9;cap=4*np.pi*a**3/3*19320*129
    scales=dict(T0_K=293.15,alpha_water_m2_s=alpha,C_Au_J_K=cap,
         tau_Au_interface_s=cap/(4*np.pi*a*a*1e8),tau_Au_internal_s=a*a/(317/(19320*129)),
         tau_diffusion_gap100nm_s=(100e-9)**2/alpha,tau_diffusion_center120nm_s=(120e-9)**2/alpha,
         tau_diffusion_gap1nm_s=(1e-9)**2/alpha,tau_diffusion_gap500nm_s=(500e-9)**2/alpha,
         laser_duration_s=10e-9,electron_phonon_reference_s=1.7e-12,
         incident_beam_radius_m=100e-6,local_incident_energy_F10_J=np.pi*(100e-6)**2*10,
         absorbed_one_particle_F10_J=3.4337324474065e-14,
         isolated_Au_rise_F10_K=23.19,radial_MR7_Au_rise_F10_K=6.81048558667788,
         adiabatic_Au_sensible_rise_F10_K=3.4337324474065e-14/cap,
         bulk_Au_melting_reference_K=1337.33,sensible_to_bulk_melting_J=cap*(1337.33-293.15),
         Tait_exponent=Tait().exponent,Tait_B_Pa=Tait().B,
         Tait_calibration='B=rho0*c0^2/n-p0 to match accepted ambient sound speed; n=7.15 reference, no temperature fit.')
    (OUT/'scales.json').write_text(json.dumps(scales,indent=2)+'\n')
    ext=[]
    for fname in ['adiabatic_summary','thermal_summary']:
        df=pd.read_csv(ROOT.parent/f'03_In_Bubble_Thermodynamics/results/{fname}.csv',float_precision='round_trip')
        for scope,g in [('observed',df),('complete',df[df.completed])]:
            r=g.loc[g.observed_Tmax_K.idxmax()].to_dict();r.update(source=f'MR3/{fname}.csv',selection=scope)
            ext.append(r)
    (OUT/'formal_extrema.json').write_text(json.dumps(ext,indent=2,default=str,allow_nan=True)+'\n')
    for stage in sorted(ROOT.parent.glob('0[1-7]_*')):
        f=next((stage/'report').glob('*.pdf'));sources.append(dict(path=str(f.relative_to(REPO)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    for f in [ROOT.parent/'03_In_Bubble_Thermodynamics/results/adiabatic_summary.csv',
              ROOT.parent/'03_In_Bubble_Thermodynamics/results/thermal_summary.csv',
              ROOT.parent/'04_Temperature_Barrier/results/species_coefficients.json',
              ROOT.parent/'05_Plasmonic_Heating/data/Johnson.yml',
              ROOT.parent/'05_Plasmonic_Heating/results/report_values.tex',
              ROOT.parent/'06_AuNP_Bubble_Coupling/results/report_values.tex',
              mr7/'stable_pairs.csv',mr7/'physical_diagnostics.csv',mr7/'transport_timing.csv']:
        sources.append(dict(path=str(f.relative_to(REPO)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    (OUT/'source_inventory.json').write_text(json.dumps(sources,indent=2)+'\n')

def main():
    inventory();mechanics();evidence()
    (OUT/'environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,platform=platform.platform()),indent=2)+'\n')
    (OUT/'preservation.json').write_text(json.dumps(inventory(),indent=2)+'\n')

if __name__=='__main__':main()
