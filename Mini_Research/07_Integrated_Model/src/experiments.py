"""Finite declared experiments; archive every run and censored status, no fitting."""
from dataclasses import asdict,replace
import argparse,json,time,platform,subprocess,sys
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from .model import Model,Config
from .particle import Particle
from .diagnostics import evaluate,extrema,pair
from .upstream import ROOT,Parameters,T0,thermo,mr4,load

RESULT=ROOT/'results'
DEFAULT=Config(cells=48,rtol=1e-7,atol=1e-10)


def table(name,rows):
    pd.DataFrame(rows).to_csv(RESULT/(name+'.csv'),index=False,float_format='%.17g')

def record(name,model,fluence=0.,center=0.,duration=10e-9,path=None,quad=True):
    start=time.time();s=model.solve(fluence,center,duration,path=path)
    r,z=evaluate(s,name,do_quad=quad);r['run_id']=name;r['runtime_s']=time.time()-start
    print(name,s.status,'Tmax',r['observed_Tmax_K'],'steps',len(s.t),'seconds',r['runtime_s'],flush=True)
    return s,r,z

class SavedPath:
    def __init__(self,base):
        self.end=base.end
        y=base.sol(base.t);R=y[0]*base.model.p.R0;U=y[1]*base.model.p.c
        self.r=CubicHermiteSpline(base.t,R,U,extrapolate=False);self.u=self.r.derivative()
    def __call__(self,t):
        if np.any(np.asarray(t)<0) or np.any(np.asarray(t)>self.end+1e-19):raise ValueError('No path extrapolation')
        return self.r(t),self.u(t)


def main_runs():
    summaries=[];pairs=[];slopes=[];pathchecks=[]
    for pa in (1.,1.15,1.3):
        cfg=replace(DEFAULT,amplitude_bar=pa,formal_material=(pa==1.3))
        p=Particle();model=Model(cfg,p);tag=f'{pa:.2f}'
        n0,rr,_=record('N0_'+tag,Model(cfg,replace(p,count=0)));summaries.append(rr);table('summary',summaries)
        base,br,bz=record('base_'+tag,model);summaries.append(br);table('summary',summaries)
        tc=br['tc_s'];center=max(1e-6,tc-5e-6)
        if not base.completed and center+5e-9>=base.end:center=base.end-1e-6
        for F in (10.,100.):
            mm=Model(replace(cfg,formal_material=True) if F==100. else cfg,p)
            heated,hr,_=record(f'laser_F{F:g}_'+tag,mm,F,center)
            summaries.append(hr);pairs.append(dict(run_id=hr['run_id'],baseline_id=br['run_id'],level=3,**pair(base,heated)))
            table('summary',summaries);table('pairs',pairs)
        for offset in (-5e-6,-1e-6,-100e-9,0.,100e-9,1e-6):
            c=tc+offset
            if c-5e-9<0 or c+5e-9>base.end:continue
            name=f'tangent_{tag}_{offset/1e-9:g}ns'
            start=time.time();s=model.sensitivity(base,c)
            r,z=evaluate(s,name,do_quad=False)
            row=dict(run_id=name,Pa_bar=pa,level=3,offset_s=offset,center_s=c,
                     baseline_id=br['run_id'],completed=base.completed,
                     baseline_Tmax_K=br['observed_Tmax_K'],baseline_t_Tmax_s=br['t_Tmax_s'],
                     dTmax_dF_K_m2_J=r['dTmax_dF_K_m2_J'],
                     incremental_balance_over_unit_absorbed=r['incremental_balance_over_unit_absorbed'])
            slopes.append(row);table('timing_sensitivity',slopes)
            print(name,'slope',row['dTmax_dF_K_m2_J'],'seconds',time.time()-start,flush=True)
        # Same no-laser mechanics, C1 interpolation checked against retained dense solution.
        path=SavedPath(base);t=base.t[:-1]+.5*np.diff(base.t);bb=base.sol(t)
        pathchecks.append(dict(Pa_bar=pa,R_max_relative=np.max(abs(path(t)[0]-bb[0]*model.p.R0))/model.p.R0,
                 U_max_m_s=np.max(abs(path(t)[1]-bb[1]*model.p.c))))
        if base.completed:
            one,b2,_=record('oneway_base_'+tag,model,path=path);summaries.append(b2)
            hs,h2,_=record('oneway_laser_'+tag,model,100,center,path=path);summaries.append(h2)
            pairs.append(dict(run_id=h2['run_id'],baseline_id=b2['run_id'],level=2,**pair(one,hs)))
            st=model.sensitivity(one,center);r,z=evaluate(st,'oneway_tangent_'+tag,do_quad=False)
            slopes.append(dict(run_id='oneway_tangent_'+tag,Pa_bar=pa,level=2,offset_s=center-tc,center_s=center,
                baseline_id=b2['run_id'],completed=one.completed,baseline_Tmax_K=b2['observed_Tmax_K'],
                baseline_t_Tmax_s=b2['t_Tmax_s'],dTmax_dF_K_m2_J=r['dTmax_dF_K_m2_J'],
                incremental_balance_over_unit_absorbed=r['incremental_balance_over_unit_absorbed']))
        if pa==1.3:
            for gas in ('air','argon'):
                drycfg=replace(cfg,gas=gas,vapor_factor=0.,formal_material=True)
                dm=Model(drycfg,p)
                ds,dr,_=record('dry_'+gas+'_base',dm);summaries.append(dr)
                dc=min(dr['tc_s']-1e-6,ds.end-1e-6)
                dh,dhr,_=record('dry_'+gas+'_laser',dm,100,dc);summaries.append(dhr)
                pairs.append(dict(run_id=dhr['run_id'],baseline_id=dr['run_id'],level=3,**pair(ds,dh)))
            _,rr,_=record('material_guard_1.30',Model(replace(cfg,formal_material=False),p));summaries.append(rr)
        table('summary',summaries);table('pairs',pairs);table('timing_sensitivity',slopes);table('mechanics_interpolation',pathchecks)


def sensitivities():
    rows=[];pairs=[]
    cfg=DEFAULT;p=Particle()
    design=[('passive_near',cfg,replace(p,motion='passive')),
            ('passive_far',cfg,replace(p,motion='passive',gap_m=1e-6)),
            ('fixed',cfg,replace(p,motion='fixed')),
            ('N8',cfg,replace(p,count=8)),
            ('a40',cfg,replace(p,radius_m=40e-9,gap_m=200e-9)),
            ('gap500',cfg,replace(p,gap_m=500e-9)),
            ('lambda600',cfg,replace(p,wavelength_m=600e-9)),
            ('G50',cfg,replace(p,G_W_m2K=5e7)),
            ('eta_half',replace(cfg,eta=.5),p),
            ('k_half',replace(cfg,k_water=.3),p),
            ('vapor5',replace(cfg,vapor_factor=5.),p),
            ('dryargon',replace(cfg,gas='argon',vapor_factor=0.,formal_material=True),p),
            ('duration100',cfg,p)]
    for name,c,part in design:
        m=Model(c,part);b,br,_=record('sens_'+name+'_base',m,quad=False);rows.append(br)
        center=max(.5e-6,min(br['tc_s']-5e-6,b.end-.1e-6))
        if not b.completed:
            table('sensitivity_runs',rows);continue
        duration=100e-9 if name=='duration100' else 10e-9
        s=m.sensitivity(b,center,duration);r,_=evaluate(s,'sens_'+name+'_tangent',do_quad=False)
        pairs.append(dict(case=name,**asdict(c),**part.metadata(),dTmax_dF_K_m2_J=r['dTmax_dF_K_m2_J'],
             observed_Tmax_K=br['observed_Tmax_K'],t_Tmax_s=br['t_Tmax_s'],status=b.status,
             incremental_balance_over_unit_absorbed=r['incremental_balance_over_unit_absorbed']))
        if name in ('passive_far','N8','duration100'):
            h,hr,_=record('sens_'+name+'_laser',m,10,center,duration);rows.append(hr)
            pairs[-1].update(pair(b,h))
        table('sensitivity_runs',rows);table('local_sensitivity',pairs)
    # Contact cases at stronger amplitude remain explicit, with zero optical source before contact.
    for mode in ('passive','fixed'):
        _,r,_=record('contact_'+mode+'_1.15',Model(replace(cfg,amplitude_bar=1.15),replace(p,motion=mode)),quad=False)
        rows.append(r)
    table('sensitivity_runs',rows)


def convergence():
    rows=[]
    design=[('coarse',replace(DEFAULT,cells=24)),('production',DEFAULT),
      ('fine',replace(DEFAULT,cells=96)),('outer50',replace(DEFAULT,padding_m=50e-6)),
      ('strict',replace(DEFAULT,rtol=2e-10,atol=1e-12,max_step_s=5e-8)),
      ('Radau',replace(DEFAULT,method='Radau',rtol=2e-9,atol=1e-11)),
      ('step_half',replace(DEFAULT,max_step_s=5e-8))]
    solutions={}
    for name,cfg in design:
        m=Model(cfg,Particle());b,br,bz=record('conv_'+name+'_base',m,quad=False)
        center=br['tc_s']-5e-6;s=m.sensitivity(b,center);r,rz=evaluate(s,'conv_'+name+'_tangent',do_quad=False)
        row=dict(check=name,**br,dTmax_dF_K_m2_J=r['dTmax_dF_K_m2_J'],
                 incremental_balance_over_unit_absorbed=r['incremental_balance_over_unit_absorbed'])
        if name in ('production','strict','Radau'):
            h,hr,_=record('conv_'+name+'_laser',m,30,center,quad=False)
            row.update(pair(b,h));row['finite_difference_slope']=row['global_delta_Tmax_K']/30
        if 'production' in solutions:
            ref=solutions['production'];tt=np.linspace(0,min(ref.end,b.end),1300)
            R,U,T,tl,_=m.state(tt,b.sol(tt));pr=ref.model;Rb,Ub,Tb,tlb,_=pr.state(tt,ref.sol(tt))
            row['R_field_maxnorm_m']=float(np.max(abs(R-Rb)));row['Tg_field_maxnorm_K']=float(np.max(abs(T-Tb)))
            # Compare liquid in the common material coordinate, not physical radii at different states.
            interp=np.array([np.interp(pr.grid.centers,m.grid.centers,tl[:,j]) for j in range(len(tt))]).T
            row['liquid_field_maxnorm_K']=float(np.max(abs(interp-tlb)))
        solutions[name]=b;rows.append(row);table('convergence',rows)
    # Regress original clamped-water ODE against the published MR4 summary.
    m=Model(replace(DEFAULT,clamp_water=True),Particle(count=0));s=m.solve();e=extrema(s)
    df=pd.read_csv(ROOT.parent/'04_Temperature_Barrier/results/summary.csv')
    r=df[df.run_id=='M2_sat_air_1.00bar_coupled_eta1'].iloc[0]
    table('MR4_regression',[dict(MR7_Tmax_K=e['observed_Tmax_K'],MR4_Tmax_K=r.observed_Tmax_K,
         relative_difference=e['observed_Tmax_K']/r.observed_Tmax_K-1,MR7_t_s=e['t_Tmax_s'],MR4_t_s=r.t_Tmax_s)])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('group',choices=['main','sensitivity','convergence']);args=parser.parse_args()
    {'main':main_runs,'sensitivity':sensitivities,'convergence':convergence}[args.group]()

if __name__=='__main__':main()
