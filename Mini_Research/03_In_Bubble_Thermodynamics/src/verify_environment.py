"""Independent locked-environment reproduction of the principal thermal cases."""
import platform
import sys
import numpy as np
import scipy
import pandas as pd
import matplotlib
from .bubble_thermo import ROOT
from .coupled_model import solve_coupled
from .one_way import solve_one_way
from .saved_mechanics import SavedMechanics
from .diagnostics import solution_metrics
from .run_experiments import json_write
from .validate import differences


def main():
    production=pd.read_csv(ROOT/'results/thermal_summary.csv');rows=[]
    for amp,mode in [(1.,'coupled'),(1.15,'coupled'),(1.3,'coupled'),(1.3,'one_way')]:
        sol=solve_coupled(amp,eta=1) if mode=='coupled' else solve_one_way(SavedMechanics('KM',amp),eta=1)
        m,_,_=solution_metrics(sol)
        label='B_fully_coupled' if mode=='coupled' else 'B_one_way_post_processing'
        ref=production[(production.Pa_bar==amp)&(production.gas=='air_like')&
            (production.eta==1)&(production.thermodynamics==label)].iloc[0].to_dict()
        row=dict(Pa_bar=amp,mode=mode,completed=m['completed'],**differences(m,ref))
        rows.append(row);print(row,flush=True)
    result=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        pandas=pd.__version__,matplotlib=matplotlib.__version__,isolated_environment=sys.prefix!=sys.base_prefix,
        lock_file='requirements-lock.txt',compared_to_production_python='3.14.4',cases=rows,
        passed=all(r['completed'] and r['relative_Tmax_error']<1e-6 for r in rows))
    json_write(ROOT/'results/clean_environment_verification.json',result)
    if not result['passed']:raise SystemExit('Independent environment check failed')


if __name__=='__main__':main()
