"""Separated domain extension, repaired photon-ledger tolerance and high-Pa controls."""
from dataclasses import replace
import numpy as np
import pandas as pd
from .upstream import ROOT
from .model import Config,Model
from .particle import Particle
from .transport import MaterialGrid
from .experiments import DEFAULT,record,table
from .diagnostics import evaluate,pair

def main():
    rows=[]
    old=MaterialGrid(DEFAULT.cells,DEFAULT.padding_m)
    h=np.log1p(DEFAULT.padding_m/5e-9)/DEFAULT.cells
    tail=5e-9*np.expm1(np.arange(DEFAULT.cells+1,100)*h)
    tail=tail[tail<50e-6];tail=np.r_[tail,50e-6]
    faces=np.r_[old.faces,((old.R0+tail)**3-old.R0**3)/3]
    design=[('outer50_fixed_inner',replace(DEFAULT,cells=len(faces)-1,padding_m=50e-6,grid_faces_m3=tuple(faces))),
       ('strict_incident_scaled',replace(DEFAULT,rtol=2e-10,atol=1e-12,max_step_s=5e-8,incident_atol=1e-6)),
       ('high_1.15_BDF',replace(DEFAULT,amplitude_bar=1.15,formal_material=True)),
       ('high_1.15_Radau',replace(DEFAULT,amplitude_bar=1.15,formal_material=True,method='Radau',rtol=2e-9,atol=1e-11)),
       ('high_1.30_BDF',replace(DEFAULT,amplitude_bar=1.3,formal_material=True)),
       ('high_1.30_Radau',replace(DEFAULT,amplitude_bar=1.3,formal_material=True,method='Radau',rtol=2e-9,atol=1e-11)),
       ('high_1.30_mesh96',replace(DEFAULT,amplitude_bar=1.3,formal_material=True,cells=96))]
    for name,cfg in design:
        m=Model(cfg,Particle());b,r,z=record('extra_'+name+'_base',m,quad=True);r['check']=name
        if cfg.amplitude_bar==1:
            center=r['tc_s']-5e-6;s=m.sensitivity(b,center)
            rr,_=evaluate(s,'extra_'+name+'_tangent',do_quad=False)
            r['dTmax_dF_K_m2_J']=rr['dTmax_dF_K_m2_J']
            if name=='strict_incident_scaled':
                h,hr,_=record('extra_'+name+'_laser',m,30,center)
                r.update(pair(b,h));r['finite_difference_slope']=r['global_delta_Tmax_K']/30
        rows.append(r);table('extra_controls',rows)

if __name__=='__main__':main()
