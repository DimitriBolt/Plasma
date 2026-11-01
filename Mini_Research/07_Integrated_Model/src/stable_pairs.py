"""Stable nonlinear paired differences on a checked archived baseline.

f(b+d)-f(b) is evaluated as integral_0^1 J(b+s*d)d ds. Three Gauss nodes
avoid subtractive cancellation. Five nodes provide an independent quadrature
check. This is a nonlinear increment, distinct from an infinitesimal tangent.
The archived no-laser base is newly reconstructed by C1 Hermite interpolation.
Its defect is shared by the pair, and its interpolation errors are recorded.
"""
from dataclasses import fields
from types import SimpleNamespace
import json,time,argparse
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
from .upstream import ROOT,T0
from .model import Model,Config,TS
from .particle import Particle
from .diagnostics import evaluate,extrema,pair
from .data_access import data_path


def archived(name):
    meta=json.loads(data_path(name,'json').read_text())
    model=Model(Config(**meta['config']),Particle(**{k:v for k,v in meta['particle'].items() if k in {f.name for f in fields(Particle)}}))
    z=np.load(data_path(name));t=z['accepted_t_s'];y=z['accepted_y'][:model.dim]
    deriv=np.column_stack([model.rhs(tt/TS,yy)/TS for tt,yy in zip(t,y.T)])
    interp=CubicHermiteSpline(t,y.T,deriv.T,axis=0,extrapolate=False)
    def dense(x):
        if np.any(np.asarray(x)<0) or np.any(np.asarray(x)>t[-1]+1e-19):raise ValueError('No extrapolation')
        return np.moveaxis(interp(x),-1,0)
    rr=meta['summary']
    base=SimpleNamespace(model=model,sol=dense,t=t,completed=rr['completed'],status=rr['status'],end=t[-1],
         requested_end=rr['requested_end_s'],fluence=0.,center=0.,duration=1e-8,path=None,tangent=False,nfev=0)
    tt=z['t_s'];delta=dense(tt)[2]*T0-z['T_g_K'];e=extrema(base)
    check=dict(base_id=name,max_interpolation_T_K=float(np.max(abs(delta))),
       peak_difference_K=e['observed_Tmax_K']-rr['observed_Tmax_K'],
       peak_time_difference_s=e['t_Tmax_s']-rr['t_Tmax_s'])
    return base,check


def solve_difference(base,F,center,duration=10e-9,order=3,oneway=False,method=None):
    m=base.model;onset=center-duration/2;offset=center+duration/2;end=base.end
    gx,gw=leggauss(order);gx=(gx+1)/2;gw=gw/2
    edges=sorted(set([0.,end]+[x for x in (onset,offset) if 0<x<end]))
    segments=[];y0=np.zeros(m.dim);status=base.status
    for a,b in zip(edges[:-1],edges[1:]):
        on=onset<(a+b)/2<offset
        P=m.abs_per_fluence*F/duration if on else 0.
        inc=np.pi*m.particle.beam_radius_m**2*F/duration if on else 0.
        def fun(t,delta):
            yy=base.sol(t*TS)
            dd=delta.copy()
            if oneway:dd[:2]=0.
            result=np.zeros(m.dim)
            for s,w in zip(gx,gw):
                point=yy+s*dd
                result+=w*m.rhs(t,point.astype(complex)+1j*1e-20*dd).imag/1e-20
            if oneway:result[:2]=0.
            if m.particle.count:result[m.au]+=P/(m.particle.count*m.particle.capacity*T0)*TS
            result[m.ledger+4:m.ledger+6]+=np.array([P,inc])/m.mix.Eref*TS
            return result
        def jac(t,y):
            # d[f(b+d)-f(b)]/dd = J(b+d); no numerical second derivative.
            point=(base.sol(t*TS)+y).astype(complex)
            columns=[]
            for j in range(m.dim):
                d=np.zeros(m.dim,complex);d[j]=1j*1e-24
                columns.append(m.rhs(t,point+d).imag/1e-24)
            J=np.array(columns).T
            if oneway:J[:2]=0.;J[:,:2]=0.
            return J
        events=[]
        for g in m.guards():
            def event(t,d,g=g):return g(t,base.sol(t*TS)+d)
            event.terminal=True;event.direction=g.direction;event.name=g.name;events.append(event)
        atols=np.full(m.dim,1e-14);atols[:3]=1e-15;atols[m.au]=1e-11
        atols[m.ledger:]=1e-12;atols[-1]=1e-6;atols*=max(F,1.)
        sol=solve_ivp(fun,(a/TS,b/TS),y0,method=method or m.cfg.method,jac=jac,dense_output=True,
                  rtol=min(m.cfg.rtol,2e-8),atol=atols,max_step=m.cfg.max_step_s/TS,events=events)
        segments.append(sol);y0=sol.y[:,-1]
        if sol.status!=0:
            status=next((g.name for g,tt in zip(events,sol.t_events) if len(tt)),sol.message);break
    end=segments[-1].t[-1]*TS
    if not base.completed and abs(end-base.end)<1e-18 and segments[-1].status==0:
        status='reference_censored_end: '+base.status
    def difference(t):
        arr=np.asarray(t);flat=arr.ravel();out=np.empty((m.dim,len(flat)))
        if np.any(flat<0) or np.any(flat>end+1e-19):raise ValueError('No extrapolation')
        ends=[s.t[-1]*TS for s in segments];ids=np.minimum(np.searchsorted(ends,flat),len(ends)-1)
        for j,s in enumerate(segments):
            mask=ids==j
            if mask.any():out[:,mask]=s.sol(flat[mask]/TS)
        return out[:,0] if arr.ndim==0 else out.reshape((m.dim,)+arr.shape)
    def dense(t):return base.sol(t)+difference(t)
    complete=base.completed and status=='completed' and abs(end-base.end)<1e-18
    s=SimpleNamespace(model=m,sol=dense,t=np.unique(np.r_[base.t[base.t<end],np.concatenate([s.t*TS for s in segments])]),
        completed=complete,status=status,end=end,requested_end=base.requested_end,fluence=F,center=center,
        duration=duration,path=None,tangent=False,nfev=sum(s.nfev for s in segments),difference=difference)
    return s


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--names',nargs='+',default=['base_1.00','base_1.15','base_1.30']);args=parser.parse_args()
    rows=pd.read_csv(ROOT/'results/stable_pairs.csv',float_precision='round_trip').to_dict('records') if (ROOT/'results/stable_pairs.csv').exists() else []
    checks=pd.read_csv(ROOT/'results/archive_interpolation.csv',float_precision='round_trip').to_dict('records') if (ROOT/'results/archive_interpolation.csv').exists() else []
    for name in args.names:
        b,check=archived(name);checks.append(check);m=b.model;e=extrema(b);tc=e['tc_s']
        design=[(False,10,3),(False,100,3),(True,10,3),(False,10,5)] if m.cfg.amplitude_bar==1 else [(False,10,3),(True,10,3)]
        for one,F,order in design:
            if F==100:object.__setattr__(m.cfg,'formal_material',True)
            center=tc-5e-6
            if center<0 or center+5e-9>b.end:continue
            label=f'stable_{name}_F{F}_L{2 if one else 3}_q{order}'
            start=time.time();s=solve_difference(b,F,center,order=order,oneway=one)
            r,z=evaluate(s,label,do_quad=True);r.update(pair(b,s));r.update(run_id=label,base_id=name,order=order,level=2 if one else 3)
            # Stable temperature increment is saved directly, independent of absolute-state subtraction.
            t=z['t_s'];d=s.difference(t)
            z['paired_delta_T_K']=d[2]*T0;z['paired_delta_R_m']=d[0]*m.p.R0
            z['paired_delta_Tl_K']=d[3:m.au]*T0;z['paired_delta_ledgers_J']=d[m.ledger:m.ledger+6]*m.mix.Eref
            np.savez_compressed(ROOT/'results/raw'/(label+'.npz'),**z)
            rows=[x for x in rows if x['run_id']!=label]+[r]
            pd.DataFrame(rows).to_csv(ROOT/'results/stable_pairs.csv',index=False,float_format='%.17g')
            pd.DataFrame(checks).to_csv(ROOT/'results/archive_interpolation.csv',index=False,float_format='%.17g')
            print(label,s.status,'delta',r['global_delta_Tmax_K'],'seconds',time.time()-start,flush=True)

if __name__=='__main__':main()
