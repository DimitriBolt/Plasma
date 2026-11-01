"""Moving-grid sensitivity at strong collapse, including the weak late-pulse tail."""
from dataclasses import replace
from .model import Model
from .particle import Particle
from .experiments import DEFAULT,record,table
from .diagnostics import evaluate

def main():
    rows=[];cfg=replace(DEFAULT,amplitude_bar=1.15,cells=96,formal_material=True)
    m=Model(cfg,Particle());b,r,z=record('high_mesh96_1.15_base',m)
    for offset in (-5e-6,-100e-9,0.):
        s=m.sensitivity(b,r['tc_s']+offset)
        rr,_=evaluate(s,f'high_mesh96_1.15_{offset/1e-9:g}ns',do_quad=False)
        rows.append(dict(Pa_bar=1.15,cells=96,offset_s=offset,dTmax_dF_K_m2_J=rr['dTmax_dF_K_m2_J'],
           T_baseline_K=r['observed_Tmax_K'],t_baseline_peak_s=r['t_Tmax_s'],
           incremental_balance_over_unit_absorbed=rr['incremental_balance_over_unit_absorbed']))
        table('high_mesh_sensitivity',rows);print(rows[-1],flush=True)

if __name__=='__main__':main()
