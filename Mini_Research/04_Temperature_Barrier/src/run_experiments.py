"""All new runs write only MR4; dependencies are imported read-only."""
import argparse
from dataclasses import asdict
import json
import math
import platform
import time
import numpy as np
import pandas as pd
import scipy
from .upstream import ROOT,SavedMechanics,old,BASE
from .thermodynamics import Mixture
from .model import solve,Settings
from .diagnostics import evaluate,peak_conditioning

OUT=ROOT/'results'
AMPLITUDES=(1.,1.15,1.3)
CONFIGS={
 'M0':dict(cv_mode='mr3',vapor_factor=0,convention='legacy'),
 'M1_const':dict(cv_mode='frozen',vapor_factor=0,convention='legacy'),
 'M1_variable':dict(cv_mode='variable',vapor_factor=0,convention='legacy'),
 'M2_dry':dict(cv_mode='variable',vapor_factor=0,convention='mixture'),
 'M2_sat':dict(cv_mode='variable',vapor_factor=1,convention='mixture'),
 'M2_sat_const':dict(cv_mode='frozen',vapor_factor=1,convention='mixture'),
 'M2_v5':dict(cv_mode='variable',vapor_factor=5,convention='mixture'),
 'M2_v10':dict(cv_mode='variable',vapor_factor=10,convention='mixture')}

def clean(x):
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [clean(v) for v in x]
    if isinstance(x,(np.bool_,)): return bool(x)
    if isinstance(x,np.integer): return int(x)
    if isinstance(x,(float,np.floating)): return float(x) if math.isfinite(x) else None
    return x

def write_json(path,obj): path.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')

def design():
    for amp in AMPLITUDES:
        for gas in ('air','argon'):
            for mode in ('one_way','coupled'):
                for config in list(CONFIGS)[:6]: yield amp,gas,mode,config,1.
                if amp==1.15 and gas=='air':
                    for eta in (.5,2.): yield amp,gas,mode,'M2_sat',eta
                if amp==1.3:
                    for config in ('M2_v5','M2_v10'): yield amp,gas,mode,config,1.
                    for eta in (.5,2.): yield amp,gas,mode,'M2_sat',eta

def run_key(amp,gas,mode,config,eta): return f'{config}_{gas}_{amp:.2f}bar_{mode}_eta{eta:g}'

def run_case(amp,gas,mode,config,eta,settings=Settings(),path=None):
    mix=Mixture(gas=gas,**CONFIGS[config]);start=time.monotonic()
    if mode=='one_way' and path is None: path=SavedMechanics('KM',amp)
    sol=solve(amp,mix,eta,mode,settings,path);m,z=evaluate(sol)
    m.update(config=config,run_id=run_key(amp,gas,mode,config,eta),runtime_s=time.monotonic()-start,
             method=settings.method,rtol=settings.rtol,steps_per_period=settings.steps_per_period)
    return m,z,sol

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--force',action='store_true')
    parser.add_argument('--only',default='');args=parser.parse_args();rows=[];paths={}
    (OUT/'trajectories').mkdir(exist_ok=True)
    for amp,gas,mode,config,eta in design():
        key=run_key(amp,gas,mode,config,eta);file=OUT/'trajectories'/f'{key}.json'
        if args.only and args.only not in key: continue
        if file.exists() and not args.force:
            m=json.loads(file.read_text())
            if m.get('Tmax_domain_K') is None: m['Tmax_domain_K']=math.inf
            with np.load(file.with_suffix('.npz')) as z: m.update(peak_conditioning(z))
            write_json(file,m)
            rows.append(m);print('cached',key,flush=True);continue
        if mode=='one_way' and amp not in paths: paths[amp]=SavedMechanics('KM',amp)
        m,z,sol=run_case(amp,gas,mode,config,eta,path=paths.get(amp))
        write_json(file,m);np.savez_compressed(OUT/'trajectories'/f'{key}.npz',**z)
        np.savez_compressed(OUT/'trajectories'/f'{key}_accepted.npz',t_s=sol.t,state=sol.y)
        rows.append(m);print(key,f"T={m['observed_Tmax_K']:.5f}",m['completed'],m['solver_message'],f"E={m['energy_residual_scaled']:.2g}",round(m['runtime_s'],1),flush=True)
        pd.DataFrame(rows).to_csv(OUT/'summary_partial.csv',index=False)
    if args.only: return
    df=pd.DataFrame(rows);df.to_csv(OUT/'summary.csv',index=False)
    reg=[];prior=pd.read_csv(ROOT.parent/'03_In_Bubble_Thermodynamics/results/thermal_summary.csv',float_precision='round_trip')
    for m in df[df.config=='M0'].to_dict('records'):
        sel=prior[(prior.Pa_bar==m['Pa_bar'])&(prior.gas==m['gas']+'_like')&(prior.eta==1)&
                  (prior.thermodynamics==('B_fully_coupled' if m['mode']=='coupled' else 'B_one_way_post_processing'))].iloc[0]
        reg.append(dict(run_id=m['run_id'],MR3_Tmax_K=sel.observed_Tmax_K,MR4_M0_Tmax_K=m['observed_Tmax_K'],
                        relative_difference=m['observed_Tmax_K']/sel.observed_Tmax_K-1,
                        time_difference_s=m['t_Tmax_s']-sel.t_Tmax_s))
    pd.DataFrame(reg).to_csv(OUT/'mr3_regression.csv',index=False)
    write_json(OUT/'run_metadata.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        pandas=pd.__version__,base=BASE,settings=asdict(Settings()),cases=len(rows),chemistry='energy assessment only',
        coefficient_stopping='Terminal boundaries; trial stages use analytic continuations only for event bracketing',
        one_way_mechanics='MR3 SavedMechanics: MR2 export with C1 Hermite radius and exact derivative velocity',
        heat_transport='Fixed MR3 conductivities, even for mixtures; uncalibrated controlled comparison'))
    (OUT/'summary_partial.csv').unlink(missing_ok=True)

if __name__=='__main__': main()
