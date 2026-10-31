"""Run the saved-parameter sweep, event comparisons and convergence audit."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
import subprocess
from dataclasses import asdict, replace
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import numpy as np
from .model import (ROOT, RP_PATH, BASE_COMMIT, Parameters, Settings, rp, solve,
                    summarize, sample_times, gas_pressure, denominator)

RESULTS = ROOT / 'results'
METRICS = ('Rmax_m','Rmin_m','Umax_positive_m_s','Umax_inward_m_s',
           'max_abs_U_m_s','max_wall_Mach','collapse_time_s','collapse_duration_s',
           'Rmax_over_Rmin','Vmax_over_Vmin')
CANONICAL = (0., .2, .5, .8, 1., 1.3)


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError('Cannot write a headerless empty CSV: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save_json(path, data):
    def cast(x):
        if isinstance(x, dict): return {k:cast(v) for k,v in x.items()}
        if isinstance(x, (list,tuple)): return [cast(v) for v in x]
        if isinstance(x, np.generic): return cast(x.item())
        if isinstance(x,float) and not np.isfinite(x): return None
        return x
    path.write_text(json.dumps(cast(data),indent=2,allow_nan=False)+'\n')


def run_key(model, bar):
    return f'{model}_{bar:.3f}bar'


def export(sol, s):
    key = run_key(sol.model, sol.Pa/1e5)
    t = sample_times(sol)
    R,U = sol.sol(t)
    p=sol.parameters
    # gzip mtime=0 makes the numeric compressed output deterministic.
    with (RESULTS/'trajectories'/f'{key}.csv.gz').open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as z:
            with io.TextIOWrapper(z,encoding='ascii',newline='') as f:
                f.write('t_s,R_m,U_m_s,p_infty_Pa,pg_Pa,D_m,wall_Mach\n')
                np.savetxt(f,np.column_stack((t,R,U,rp.external_pressure(t,sol.Pa,p),
                           gas_pressure(R,p),denominator(R,U,p),abs(U)/p.c)),
                           delimiter=',',fmt='%.17g')
    np.savez_compressed(RESULTS/'trajectories'/f'{key}_accepted.npz',
                        t_s=sol.t,R_m=sol.y[0],U_m_s=sol.y[1])
    save_json(RESULTS/'trajectories'/f'{key}.json',dict(
        summary=s,parameters=asdict(p),settings=asdict(sol.settings),events=sol.events,
        interpolation='CSV samples evaluated on accepted-step dense polynomials; no reintegration',
        columns='CSV column suffixes are SI units; wall_Mach is dimensionless',
        gas_law='pg0*(R0/R)**(3*kappa); pg0=p0+2*sigma/R0-pv',
        numerical_warning='Numerical convergence does not establish physical validity'))


def sweep(force=False):
    rows, events = [], []
    def run(bar,model):
        key=run_key(model,bar)
        path=RESULTS/'trajectories'/f'{key}.json'
        if path.exists() and not force:
            old=json.loads(path.read_text()); s=old['summary']
            s={k:np.nan if v is None else v for k,v in s.items()}
            ev=old['events']
        else:
            sol=solve(bar*1e5,model)
            s=summarize(sol); ev=sol.events
            export(sol,s)
            print(f"{model} {bar:.3f} bar: Rmin={s['Rmin_m']*1e6:.6g} um, "
                  f"M={s['max_wall_Mach']:.6g}, complete={s['completed']}",flush=True)
        rows.append(s)
        events.extend(dict(model=model,Pa_bar=bar,**e) for e in ev)
        return s
    coarse=np.round(np.arange(1.,1.30001,.05),3)
    pilot=[]
    for bar in sorted(set(CANONICAL)|set(coarse)):
        for model in ('RP','KM'):
            s=run(bar,model)
            if bar>=1.: pilot.append(s)
    write_csv(RESULTS/'preliminary_scan.csv',pilot)
    # Actual pilot behavior selects refinement: a log-radius slope >4 per bar
    # (~22% change per .05 bar) is an operational mesh rule, not a bifurcation.
    intervals=[]
    for a,b in zip(coarse[:-1],coarse[1:]):
        slopes=[]
        for model in ('RP','KM'):
            ra=next(r for r in pilot if r['model']==model and r['Pa_bar']==a)
            rb=next(r for r in pilot if r['model']==model and r['Pa_bar']==b)
            slopes.append(abs(np.log(rb['Rmin_m']/ra['Rmin_m']))/(b-a))
        if max(slopes)>4:
            intervals.append((float(a),float(b)))
    added=sorted({round(float(x),3) for a,b in intervals for x in np.arange(a,b+.0001,.005)}
                 -{r['Pa_bar'] for r in rows})
    for bar in added:
        for model in ('RP','KM'): run(bar,model)
    rows.sort(key=lambda r:(r['Pa_bar'],r['model']))
    events.sort(key=lambda r:(r['Pa_bar'],r['model'],r['time_s']))
    write_csv(RESULTS/'amplitude_summary.csv',rows)
    write_csv(RESULTS/'extrema.csv',events)
    save_json(RESULTS/'scan_design.json',dict(
        pilot_step_bar=.05,range_bar=[1.,1.3],refined_step_bar=.005,
        refinement_rule='max_model abs(delta log Rmin)/delta Pa_bar >4',
        refined_intervals_bar=intervals,amplitudes_bar=sorted({r['Pa_bar'] for r in rows}),
        interpretation='Finite five-cycle mesh; no bifurcation or unique critical pressure claimed'))
    comparisons=[]
    for bar in sorted({r['Pa_bar'] for r in rows}):
        a=next(r for r in rows if r['Pa_bar']==bar and r['model']=='RP')
        b=next(r for r in rows if r['Pa_bar']==bar and r['model']=='KM')
        valid=a['completed'] and b['completed']
        row=dict(Pa_bar=bar,RP_completed=a['completed'],KM_completed=b['completed'])
        for k in METRICS:
            row['RP_'+k]=a[k]; row['KM_'+k]=b[k]
            row['KM_minus_RP_'+k]=b[k]-a[k] if valid else np.nan
            row['relative_difference_'+k]=(b[k]-a[k])/abs(a[k]) if a[k] and valid else np.nan
        comparisons.append(row)
    write_csv(RESULTS/'comparison.csv',comparisons)
    # Corresponding primary event per forcing cycle, independently of which
    # cycle provides the global deepest minimum for each model.
    pairs=[]
    for bar in sorted({r['Pa_bar'] for r in rows}):
        for cycle in range(1,6):
            row=dict(Pa_bar=bar,cycle=cycle)
            selected={}
            for model in ('RP','KM'):
                summ=next(r for r in rows if r['model']==model and r['Pa_bar']==bar)
                # A cycle crossing the stopping time is censored. Earlier full
                # cycles may still be compared, even on an incomplete run.
                cycle_complete=summ['last_time_s'] >= cycle*Parameters().T
                row[model+'_cycle_complete']=bool(cycle_complete)
                mins=[e for e in events if cycle_complete and e['model']==model and e['Pa_bar']==bar
                      and e['kind']=='minimum' and int(e['time_s']/Parameters().T)+1==cycle]
                e=min(mins,key=lambda e:e['R_m']) if mins else None
                maxs=[x for x in events if x['model']==model and x['Pa_bar']==bar
                      and x['kind']=='maximum' and e and x['time_s']<e['time_s']]
                prev=maxs[-1] if maxs else None
                row[model+'_time_s']=e['time_s'] if e else np.nan
                row[model+'_Rmin_m']=e['R_m'] if e else np.nan
                row[model+'_preceding_Rmax_m']=prev['R_m'] if prev else np.nan
                row[model+'_duration_s']=e['time_s']-prev['time_s'] if prev else np.nan
                selected[model]=e
            row['KM_minus_RP_time_s']=row['KM_time_s']-row['RP_time_s']
            row['relative_Rmin_difference']=(row['KM_Rmin_m']-row['RP_Rmin_m'])/row['RP_Rmin_m']
            pairs.append(row)
    write_csv(RESULTS/'cycle_comparison.csv',pairs)
    # Compare full CSV precision with the accepted reference values.
    baseline=[]
    for old in read_csv(ROOT.parent/'01_Bubble_Dynamics/results/amplitude_summary.csv'):
        current=next(r for r in rows if r['model']=='RP' and r['Pa_bar']==float(old['Pa_bar']))
        for k in ('Rmax_m','Rmin_m','max_abs_U_m_s','collapse_time_s'):
            x=float(old[k]); y=current[k]
            baseline.append(dict(Pa_bar=current['Pa_bar'],metric=k,accepted=x,reproduced=y,
                                 absolute_error=abs(y-x),relative_error=abs(y-x)/abs(x) if x else np.nan))
    write_csv(RESULTS/'rp_reproduction.csv',baseline)
    return rows


def convergence():
    rows=[]
    configs=[Settings(rtol=1e-8,atol_R=1e-14,atol_U=1e-6,steps_per_period=300),
             Settings(),Settings(rtol=1e-11,atol_R=1e-17,atol_U=1e-9,steps_per_period=2000),
             Settings(method='Radau')]
    for model in ('RP','KM'):
        for bar in (1.,1.1,1.15,1.3):
            group=[]
            for cfg in configs:
                sol=solve(bar*1e5,model,settings=cfg)
                s=summarize(sol)
                row=dict(method=cfg.method,rtol=cfg.rtol,atol_R_m=cfg.atol_R,
                         atol_U_m_s=cfg.atol_U,max_step_s=Parameters().T/cfg.steps_per_period,**s)
                group.append(row)
                print(f'Convergence {model} {bar} {cfg.method} {cfg.rtol:g}',flush=True)
            ref=group[2]
            for row in group:
                for k in ('Rmax_m','Rmin_m','max_abs_U_m_s','collapse_time_s','collapse_duration_s'):
                    row['relative_error_'+k]=abs(row[k]-ref[k])/abs(ref[k])
                row['collapse_time_error_s']=row['collapse_time_s']-ref['collapse_time_s']
            rows.extend(group)
            write_csv(RESULTS/'convergence.csv',rows)


def c_limit():
    rows=[]
    cfg=Settings(rtol=1e-11,atol_R=1e-17,atol_U=1e-9,steps_per_period=2000)
    for bar in (1.,1.3):
        p=Parameters(); ref=solve(bar*1e5,'RP',p,cfg); sr=summarize(ref)
        for factor in (1.,10.,100.,1000.,1e5,1e7):
            pc=replace(p,c=p.c*factor); sol=solve(bar*1e5,'KM',pc,cfg); s=summarize(sol)
            # Union of both dense event-aware grids prevents a false convergence
            # conclusion from a common uniform grid missing the RP narrow core.
            t=np.unique(np.r_[sample_times(ref),sample_times(sol)])
            a,b=ref.sol(t),sol.sol(t)
            row=dict(Pa_bar=bar,c_factor=factor,c_m_s=pc.c,completed=s['completed'],
                     max_R_difference_over_R0=float(abs(a[0]-b[0]).max()/p.R0),
                     max_U_difference_over_reference_max=float(abs(a[1]-b[1]).max()/sr['max_abs_U_m_s']),
                     Rmin_m=s['Rmin_m'],RP_Rmin_m=sr['Rmin_m'],
                     relative_Rmin_error=abs(s['Rmin_m']-sr['Rmin_m'])/sr['Rmin_m'],
                     max_abs_U_m_s=s['max_abs_U_m_s'],
                     relative_max_speed_error=abs(s['max_abs_U_m_s']-sr['max_abs_U_m_s'])/sr['max_abs_U_m_s'],
                     collapse_time_shift_s=s['collapse_time_s']-sr['collapse_time_s'])
            rows.append(row)
            print(f"c limit {bar} x{factor:g}: Rmin rel={row['relative_Rmin_error']:.3g}",flush=True)
            write_csv(RESULTS/'c_limit.csv',rows)


def metadata():
    source_root=ROOT.parents[1]
    files=[RP_PATH,source_root/'Plasma-plan.tex',*[source_root/'Library'/n for n in
           ('035218_1_online.pdf','nature03361.pdf','PhysRevLett.104.244301.pdf')]]
    for f in files:
        if not f.is_file() or not f.stat().st_size: raise FileNotFoundError(f)
    def sha(f): return hashlib.sha256(f.read_bytes()).hexdigest()
    assert sha(RP_PATH)=='f8c1baaf5a4fc891383b7203b67c026908f4bf005d5418e273b640784e161188'
    try: head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    except subprocess.CalledProcessError: head=None
    save_json(RESULTS/'run_metadata.json',dict(
        created_utc=datetime.now(timezone.utc).isoformat(),study_started_client_date='2026-10-30',
        source_repository_commit=head,accepted_RP_commit=BASE_COMMIT,
        parameters=asdict(Parameters()),initial_liquid_temperature_K=293.15,
        initial_state=[Parameters().R0,0.],settings=asdict(Settings()),
        python=platform.python_version(),platform=platform.platform(),
        libraries={x:version(x) for x in ('numpy','scipy','matplotlib','pytest')},
        source_sha256={str(f.relative_to(source_root)):sha(f) for f in files},
        code_sha256={str(f.relative_to(ROOT)):sha(f) for f in sorted((ROOT/'src').glob('*.py'))},
        scope='Mechanical radius dynamics only; temperature and reactions are not evaluated'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=('all','sweep','convergence','c-limit','metadata'),default='all')
    parser.add_argument('--force',action='store_true')
    args=parser.parse_args()
    RESULTS.mkdir(exist_ok=True); (RESULTS/'trajectories').mkdir(exist_ok=True)
    if args.stage in ('all','sweep'): sweep(args.force)
    if args.stage in ('all','convergence'): convergence()
    if args.stage in ('all','c-limit'): c_limit()
    metadata()


if __name__=='__main__': main()
