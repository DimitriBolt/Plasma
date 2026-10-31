"""Bounded numerical checks, independent method and offline quadrature."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import argparse
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from .upstream import ROOT,SavedMechanics
from .model import Settings
from .run_experiments import run_case,write_json

OUT=ROOT/'results'

def quadrature(sol,tp):
    mix=sol.mix
    records=[]
    for end,label in ((tp,'Tmax'),(sol.t[-1],'end')):
        edges=np.unique(np.r_[sol.t[sol.t<end],end]);mid=(edges[:-1]+edges[1:])/2;half=np.diff(edges)/2
        for order in (8,16):
            nodes,weights=leggauss(order);ts=mid[:,None]+half[:,None]*nodes
            R,U,T=sol.state(ts.ravel(),sol.sol(ts.ravel()));pg=mix.nt*8.31446261815324*T/(4*np.pi*R**3/3)
            work=(pg*4*np.pi*R**2*U).reshape(ts.shape)
            loss=(4*np.pi*sol.eta*mix.kg*R*(T-293.15)).reshape(ts.shape)
            W=float(np.sum(work*weights*half[:,None]));Q=float(np.sum(loss*weights*half[:,None]));du=mix.du(float(sol.state(end,sol.sol(end))[2]))
            records.append(dict(endpoint=label,order=order,delta_U_J=du,quadrature_work_J=W,
                                quadrature_heat_J=Q,residual_J=du+W+Q,
                                scaled_residual=abs(du+W+Q)/max(mix.Eref,abs(du),abs(W),abs(Q))))
    return records

def check_case(task):
    base,kind=task;settings=(Settings(rtol=2e-11,steps_per_period=1200) if kind!='Radau' else Settings(method='Radau',rtol=2e-10,steps_per_period=900))
    path=None
    if base['mode']=='one_way':
        path=SavedMechanics('KM',base['Pa_bar'],stride=2 if kind=='input_stride2' else 1,
                            grid='accepted' if kind=='input_accepted' else 'export')
    m,z,sol=run_case(base['Pa_bar'],base['gas'],base['mode'],base['config'],base['eta'],settings,path)
    record=dict(run_id=base['run_id'],check=kind,method=settings.method,rtol=settings.rtol,
                completed=m['completed'],solver_message=m['solver_message'],end_s=m['end_s'],
                production_Tmax_K=base['observed_Tmax_K'],check_Tmax_K=m['observed_Tmax_K'],
                peak_relative_difference=m['observed_Tmax_K']/base['observed_Tmax_K']-1,
                peak_time_difference_s=m['t_Tmax_s']-base['t_Tmax_s'],
                end_time_difference_s=m['end_s']-base['end_s'],
                energy_residual_scaled=m['energy_residual_scaled'],
                extrema_scope=m['extrema_scope'])
    quad=[]
    if kind=='strict':
        for q in quadrature(sol,m['t_Tmax_s']): quad.append(dict(run_id=base['run_id'],**q))
    return record,quad

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--force',action='store_true');parser.add_argument('--available',action='store_true');args=parser.parse_args()
    if args.available:
        import json
        df=pd.DataFrame([json.loads(p.read_text()) for p in sorted((OUT/'trajectories').glob('*.json'))])
    else: df=pd.read_csv(OUT/'summary.csv',float_precision='round_trip')
    main=df[(df.config.isin(['M1_variable','M2_sat']))&(df.eta==1)]
    tasks=[(r,'strict') for r in main.to_dict('records')]
    tasks += [(r,'Radau') for r in main[main.config=='M2_sat'].to_dict('records')]
    tasks += [(r,kind) for r in main[(main.Pa_bar==1.3)&(main['mode']=='one_way')].to_dict('records') for kind in ('input_stride2','input_accepted')]
    # M0 checks distinguish exact equation recovery from the accuracy of new integration.
    tasks += [(r,'strict') for r in df[df.config=='M0'].to_dict('records')]
    rows=[];quads=[];(OUT/'validation_cases').mkdir(exist_ok=True);pending=[]
    for task in tasks:
        base,kind=task;file=OUT/'validation_cases'/(base['run_id']+'_'+kind+'.json')
        if file.exists() and not args.force:
            data=__import__('json').loads(file.read_text());rows.append(data['summary']);quads.extend(data['quadrature'])
        else: pending.append(task)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for task,(row,quad) in zip(pending,pool.map(check_case,pending)):
            rows.append(row);quads.extend(quad)
            write_json(OUT/'validation_cases'/(row['run_id']+'_'+row['check']+'.json'),dict(summary=row,quadrature=quad))
            print(row['check'],row['run_id'],row['peak_relative_difference'],flush=True)
    pd.DataFrame(rows).to_csv(OUT/'convergence.csv',index=False)
    pd.DataFrame(quads).to_csv(OUT/'independent_quadrature.csv',index=False)
    time_selection_audit()

def time_selection_audit():
    """Keep large argmax-time changes visible when successive peaks nearly tie."""
    df=pd.read_csv(OUT/'summary.csv',float_precision='round_trip')
    conv=pd.read_csv(OUT/'convergence.csv',float_precision='round_trip');rows=[]
    for r in conv.to_dict('records'):
        base=df[df.run_id==r['run_id']].iloc[0];dt=r['peak_time_difference_s']
        checktime=base.t_Tmax_s+dt
        neartie=(base.near_peak_candidate_count>1 and base.near_peak_time_first_s-1e-9<=checktime<=base.near_peak_time_last_s+1e-9)
        rows.append(dict(run_id=r['run_id'],check=r['check'],peak_time_difference_s=dt,
                         production_peak_s=base.t_Tmax_s,check_peak_s=checktime,
                         production_near_peak_candidates=base.near_peak_candidate_count,
                         near_peak_first_s=base.near_peak_time_first_s,near_peak_last_s=base.near_peak_time_last_s,
                         near_peak_tolerance_relative=base.near_peak_relative_tolerance,
                         interpretation='near_degenerate_cycle_selection' if abs(dt)>1e-9 and neartie else 'same_event_numerical_comparison'))
    pd.DataFrame(rows).to_csv(OUT/'peak_time_conditioning.csv',index=False)

if __name__=='__main__': main()
