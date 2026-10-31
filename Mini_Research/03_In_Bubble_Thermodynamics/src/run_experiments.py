"""Read accepted MR2 outputs; write ONLY the MR3 reproducibility package."""
import argparse
from dataclasses import asdict,replace
import json
import platform
import time
import numpy as np
import pandas as pd
import scipy
from .bubble_thermo import (ROOT,MR2,Gas,Parameters,MECHANICAL_COMMIT,MR1_COMMIT,
    adiabatic_temperature,adiabatic_pressure,volume,gas_density)
from .saved_mechanics import audit_inputs,SavedMechanics
from .diagnostics import adiabatic_metrics,solution_metrics,save_series
from .coupled_model import solve_coupled
from .one_way import solve_one_way
from .numerics import Settings

OUT=ROOT/'results'
AMPLITUDES=(1.,1.15,1.3)


def json_write(file,obj):
    def clean(x):
        if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(tuple,list)):return [clean(v) for v in x]
        if isinstance(x,(np.bool_,)):return bool(x)
        if isinstance(x,(np.integer,)):return int(x)
        if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
        return x
    file.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')


def adiabatic_stage():
    rows=[];events=[];hotrows=[];sensitivity=[]
    admission=pd.read_csv(OUT/'input_admission.csv')
    for a in admission.to_dict('records'):
        path=SavedMechanics(a['model'],a['Pa_bar'])
        selected=a['Pa_bar'] in AMPLITUDES
        if selected:
            m,hot=adiabatic_metrics(path,details=True)
            for h in hot:
                hotrows.append({**{k:m[k] for k in ['mechanics','Pa_bar','completed','extrema_scope']},
                                **{k:v for k,v in h.items() if k!='intervals'}})
        else:m=adiabatic_metrics(path)
        m['upstream_next_stage_use']=a['next_stage_use'];m['upstream_numerical_status']=a['numerical_status']
        rows.append(m)
        if selected:
            T=adiabatic_temperature(path.export.R_m)
            np.savez_compressed(OUT/f'A_{path.key}.npz',t_s=path.export.t_s,T_K=T,
                pg_Pa=adiabatic_pressure(path.export.R_m),V_m3=volume(path.export.R_m),
                density_kg_m3=gas_density(path.export.R_m))
            for e in path.events:
                if e['kind'] not in ('minimum','maximum'):continue
                events.append(dict(mechanics=a['model'],Pa_bar=a['Pa_bar'],**e,
                    cycle=int(e['time_s']*Parameters().f)+1,T_K=float(adiabatic_temperature(e['R_m'])),
                    first_strong=(e['kind']=='minimum' and e['time_s']==path.summary['first_strong_time_s']),
                    global_deepest=(path.completed and e['kind']=='minimum' and e['R_m']==path.summary['Rmin_m'])))
        if selected and a['model']=='KM':
            for gamma in (1.2,1.4,1.55,5/3):
                for T0 in (273.15,293.15,313.15):
                    gas=replace(Gas(),gamma=gamma,T0=T0)
                    r=adiabatic_metrics(path,gas);r['fixed_mechanics']=True
                    r['fixed_water_properties']=True;r['n_mol']=gas.n();r['Cv_J_mol_K']=gas.Cv
                    r['pressure_mismatch_at_peak']=float(adiabatic_pressure(r['R_at_Tmax_m'],gas)/(
                        Parameters().pg0*(Parameters().R0/r['R_at_Tmax_m'])**(3*1.4))-1)
                    sensitivity.append(r)
        print('A',path.key,m['observed_Tmax_K'],m['completed'],flush=True)
    pd.DataFrame(rows).to_csv(OUT/'adiabatic_summary.csv',index=False)
    pd.DataFrame(events).to_csv(OUT/'collapse_temperatures.csv',index=False)
    pd.DataFrame(hotrows).to_csv(OUT/'adiabatic_hot_durations.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/'adiabatic_sensitivity.csv',index=False)
    # Explicit selections preserve the upstream first-strong/cycle/global distinction.
    cyc=pd.read_csv(MR2/'results/cycle_comparison.csv')
    cyc=cyc[cyc.Pa_bar.isin(AMPLITUDES)].copy()
    for model in ('RP','KM'):
        cyc[model+'_T_cycle_deepest_K']=adiabatic_temperature(cyc[model+'_Rmin_m'])
    cyc.to_csv(OUT/'RP_KM_cycle_temperature_comparison.csv',index=False)


def thermal_stage(force=False):
    rows=[];intervals=[]
    for amplitude in AMPLITUDES:
        path=SavedMechanics('KM',amplitude)
        gases=[(Gas(),eta) for eta in (0.,.1,1.,10.)]
        gases.append((Gas(gamma=5/3,kg=.0177,molar_mass=.039948,name='argon_like'),1.))
        for gas,eta in gases:
            for mode in ('one_way','coupled'):
                key=f'{mode}_KM_{amplitude:.3f}bar_{gas.name}_eta{eta:g}'
                if not force and (OUT/(key+'.json')).exists():
                    cached=json.loads((OUT/(key+'.json')).read_text())
                    if mode=='coupled' or cached.get('interpolation_version')=='round_trip_v1':
                        rows.append(cached['summary'])
                        for h in cached['hot_intervals']:
                            intervals.append(dict(run_id=key,**{k:v for k,v in h.items() if k!='intervals'}))
                        print('reused',key,flush=True)
                        continue
                start=time.monotonic()
                sol=solve_one_way(path,gas,eta) if mode=='one_way' else solve_coupled(amplitude,gas,eta)
                m,hot,grid=solution_metrics(sol);m['run_id']=key
                m['runtime_s']=time.monotonic()-start
                save_series(OUT/(key+'.npz'),sol,grid)
                np.savez_compressed(OUT/(key+'_accepted.npz'),t_s=sol.t,state=sol.y)
                json_write(OUT/(key+'.json'),dict(summary=m,gas=asdict(gas),eta=eta,
                    settings=asdict(sol.settings),source_commit=MECHANICAL_COMMIT,hot_intervals=hot,
                    interpolation_version='round_trip_v1'))
                rows.append(m)
                for h in hot:
                    intervals.append(dict(run_id=key,**{k:v for k,v in h.items() if k!='intervals'}))
                print(key,m['observed_Tmax_K'],'energy',m['energy_residual_scaled'],
                      'completed',m['completed'],'seconds',round(m['runtime_s'],1),flush=True)
                pd.DataFrame(rows).to_csv(OUT/'thermal_summary.csv',index=False)
    pd.DataFrame(rows).to_csv(OUT/'thermal_summary.csv',index=False)
    pd.DataFrame(intervals).to_csv(OUT/'thermal_hot_durations.csv',index=False)


def metadata():
    p=Parameters();g=Gas()
    json_write(OUT/'parameter_provenance.json',dict(source_commit=MECHANICAL_COMMIT,MR1_commit=MR1_COMMIT,
        water=asdict(p),gas=asdict(g),pg0_Pa=p.pg0,n_mol=g.n(),Cv_J_mol_K=g.Cv,
        E0_J=g.E0(),conductivity_source='Huber, Thermal Conductivity of Gases, NIST table, 300 K',
        conductivity_url='https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=926438',
        conductivity_note='Air 0.0264 W/(m K) at 100 kPa; argon 0.0177 W/(m K), dilute limit. '
            '300 K values held constant even though T0=293.15 K. Phenomenological sensitivity, not high-T properties.',
        gas_inventory_note='n fixed by equilibrium for each T0; same n for air and argon at equal initial state; '
            'molar mass affects mass density, not the number of moles.',
        eta_values=[0,.1,1,10],eta_note='Uncalibrated decade sensitivity around geometry estimate eta=1.',
        source_power_W=0,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        pandas=pd.__version__,interpolation='New C1 cubic Hermite R(t), nodal slopes=stored U; U(t)=dR/dt; '
            'event-aware CSV primary, accepted nodes and thinned exports are resolution checks; no extrapolation.'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['audit','adiabatic','thermal','metadata','all'],default='all')
    parser.add_argument('--force',action='store_true',help='Recompute cached thermal runs')
    args=parser.parse_args();OUT.mkdir(exist_ok=True)
    if args.stage in ('audit','all'):audit_inputs()
    if args.stage in ('adiabatic','all'):adiabatic_stage()
    if args.stage in ('thermal','all'):thermal_stage(args.force)
    if args.stage in ('metadata','all'):metadata()


if __name__=='__main__':main()
