"""Extended mathematical-limit and independent clean-environment checks."""
from dataclasses import asdict, replace
from importlib.metadata import version
import json
import platform
import numpy as np
from .model import Parameters, Settings, solve, summarize, sample_times
from .experiments import (RESULTS, CANONICAL, read_csv, write_csv, save_json)


def extended_limit():
    cfg=Settings(rtol=1e-11,atol_R=1e-17,atol_U=1e-9,steps_per_period=2000)
    p=Parameters(); ref=solve(1.3e5,'RP',p,cfg); sr=summarize(ref); rows=[]
    for factor in (1e9,1e11,1e13,1e15):
        pc=replace(p,c=p.c*factor); sol=solve(1.3e5,'KM',pc,cfg); s=summarize(sol)
        t=np.unique(np.r_[sample_times(ref),sample_times(sol)]); a,b=ref.sol(t),sol.sol(t)
        rows.append(dict(Pa_bar=1.3,c_factor=factor,c_m_s=pc.c,completed=s['completed'],
            max_R_difference_over_R0=float(abs(a[0]-b[0]).max()/p.R0),
            max_U_difference_over_reference_max=float(abs(a[1]-b[1]).max()/sr['max_abs_U_m_s']),
            Rmin_m=s['Rmin_m'],RP_Rmin_m=sr['Rmin_m'],
            relative_Rmin_error=abs(s['Rmin_m']-sr['Rmin_m'])/sr['Rmin_m'],
            max_abs_U_m_s=s['max_abs_U_m_s'],
            relative_max_speed_error=abs(s['max_abs_U_m_s']-sr['max_abs_U_m_s'])/sr['max_abs_U_m_s'],
            collapse_time_shift_s=s['collapse_time_s']-sr['collapse_time_s']))
        write_csv(RESULTS/'c_limit_extended.csv',rows)
        print(f'Extended c limit x{factor:g}: radius trajectory error {rows[-1]["max_R_difference_over_R0"]:.3g}',flush=True)


def reproduce():
    original=read_csv(RESULTS/'amplitude_summary.csv') if (RESULTS/'amplitude_summary.csv').exists() else [
        json.loads(path.read_text())['summary'] for path in (RESULTS/'trajectories').glob('*.json')]
    rows=[]
    cfg=Settings(rtol=1e-11,atol_R=1e-17,atol_U=1e-9,steps_per_period=1600)
    for bar in CANONICAL:
        for model in ('RP','KM'):
            saved=next(r for r in original if r['model']==model and float(r['Pa_bar'])==bar)
            sol=solve(bar*1e5,model,settings=cfg); s=summarize(sol)
            for k in ('Rmax_m','Rmin_m','max_abs_U_m_s','collapse_time_s'):
                ref=float(saved[k]) if saved[k] is not None else np.nan; value=s[k]
                rows.append(dict(model=model,Pa_bar=bar,metric=k,production=ref,reproduced=value,
                                 absolute_error=abs(value-ref),
                                 relative_error=abs(value-ref)/abs(ref) if ref else np.nan,
                                 completed=s['completed']))
            print(f'Clean repetition {model} {bar:g}',flush=True)
    write_csv(RESULTS/'clean_reproduction.csv',rows)
    save_json(RESULTS/'clean_environment.json',dict(python=platform.python_version(),
        versions={x:version(x) for x in ('numpy','scipy','matplotlib','pytest')},
        settings=asdict(cfg),description='Independent settings for all six original amplitudes'))


def strongest_rp():
    original=read_csv(RESULTS/'amplitude_summary.csv')
    complete=[r for r in original if r['model']=='RP' and r['completed']=='True']
    row=max(complete,key=lambda r:float(r['max_abs_U_m_s']))
    bar=float(row['Pa_bar']); rows=[]
    for cfg in (Settings(rtol=1e-11,atol_R=1e-17,atol_U=1e-9,steps_per_period=2000),Settings(method='Radau')):
        sol=solve(bar*1e5,'RP',settings=cfg); s=summarize(sol)
        rows.append(dict(method=cfg.method,rtol=cfg.rtol,atol_R_m=cfg.atol_R,
                         atol_U_m_s=cfg.atol_U,max_step_s=Parameters().T/cfg.steps_per_period,**s))
        print(f'Strongest complete RP {bar:g} bar {cfg.method}: complete={s["completed"]}',flush=True)
    write_csv(RESULTS/'strongest_RP_validation.csv',rows)


if __name__=='__main__':
    extended_limit(); reproduce(); strongest_rp()
