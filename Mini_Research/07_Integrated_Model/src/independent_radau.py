"""Full nonlinear independent finite-difference check at a resolved weak source."""
from dataclasses import replace
import pandas as pd
from .model import Model
from .particle import Particle
from .experiments import DEFAULT,record,table
from .diagnostics import pair,evaluate

def main():
    cfg=replace(DEFAULT,method='Radau',rtol=2e-11,atol=1e-13,incident_atol=1e-6)
    m=Model(cfg,Particle());b,br,_=record('independent_Radau_base',m)
    h,hr,_=record('independent_Radau_laser',m,10,br['tc_s']-5e-6)
    table('independent_Radau',[dict(**pair(b,h),baseline_id=br['run_id'],laser_id=hr['run_id'],
         baseline_total_residual_scaled=br['total_residual_scaled'],laser_total_residual_scaled=hr['total_residual_scaled'])])
    s=m.sensitivity(b,br['tc_s']-5e-6);r,z=evaluate(s,'independent_Radau_tangent',do_quad=False)
    table('independent_Radau_sensitivity',[r])

if __name__=='__main__':main()
