"""Bounded fresh one-cycle replay against the archived five-cycle baseline."""
from dataclasses import replace
import json
import numpy as np
from .experiments import DEFAULT
from .model import Model
from .particle import Particle
from .stable_pairs import archived
from .upstream import ROOT,T0

def main():
    b,_=archived('base_1.00');m=Model(replace(DEFAULT,cycles=1.),Particle());s=m.solve()
    tt=np.linspace(0,.95/m.p.f,1001);a=b.sol(tt);v=s.sol(tt)
    result=dict(completed=bool(s.completed),interval_s=[0,float(tt[-1])],samples=len(tt),
        R_max_difference_m=float(np.max(abs(a[0]-v[0]))*m.p.R0),
        Tg_max_difference_K=float(np.max(abs(a[2]-v[2]))*T0),
        water_max_difference_K=float(np.max(abs(a[3:m.au]-v[3:m.au]))*T0),
        scope='fresh one-cycle no-laser replay; not a rerun of all sweeps')
    assert s.completed and result['Tg_max_difference_K']<.01 and result['R_max_difference_m']<1e-9
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
