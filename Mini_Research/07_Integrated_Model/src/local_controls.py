"""Unconfounded radius and wet-gas composition controls for the finite local design."""
from dataclasses import replace,asdict
import pandas as pd
from .upstream import ROOT
from .model import Model
from .particle import Particle
from .experiments import DEFAULT,record,table
from .diagnostics import evaluate

def main():
    cases=pd.read_csv(ROOT/'results/local_sensitivity.csv',float_precision='round_trip').to_dict('records')
    runs=pd.read_csv(ROOT/'results/sensitivity_runs.csv',float_precision='round_trip').to_dict('records')
    for name,cfg,p in [('a40_fixed_gap',DEFAULT,Particle(radius_m=40e-9,gap_m=100e-9)),
                       ('wetargon',replace(DEFAULT,gas='argon'),Particle())]:
        m=Model(cfg,p);b,br,_=record('sens_'+name+'_base',m,quad=False);runs.append(br)
        s=m.sensitivity(b,br['tc_s']-5e-6);r,z=evaluate(s,'sens_'+name+'_tangent',do_quad=False)
        cases.append(dict(case=name,**asdict(cfg),**p.metadata(),dTmax_dF_K_m2_J=r['dTmax_dF_K_m2_J'],
          observed_Tmax_K=br['observed_Tmax_K'],t_Tmax_s=br['t_Tmax_s'],status=b.status,
          incremental_balance_over_unit_absorbed=r['incremental_balance_over_unit_absorbed']))
        table('local_sensitivity',cases);table('sensitivity_runs',runs)

if __name__=='__main__':main()
