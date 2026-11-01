"""Dense-output events, independent flux quadrature, SI archives and pair scope."""
from dataclasses import asdict
import json
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from .upstream import ROOT,T0
from .model import TS
from .transport import center_radius

def roots(fun,t,values=None):
    v=np.array([fun(x) for x in t]) if values is None else np.asarray(values)
    ids=np.flatnonzero(v[:-1]*v[1:]<0)
    r=[brentq(fun,t[j],t[j+1],xtol=1e-23,rtol=1e-14) for j in ids]
    return np.unique(np.r_[r,t[v==0]])

def extrema(s,end=None):
    m=s.model;end=s.end if end is None else end
    nodes=s.t[s.t<=end];nodes=np.unique(np.r_[nodes,end])
    t=np.unique(np.r_[nodes,nodes[:-1]+.25*np.diff(nodes),nodes[:-1]+.5*np.diff(nodes),nodes[:-1]+.75*np.diff(nodes)])
    state=lambda x:m.state(x,s.sol(x)[:m.dim],s.path)
    def drate(x):
        return m.rhs(x/TS,s.sol(x)[:m.dim],path=s.path)[2]
    tex=roots(drate,t)
    rex=roots(lambda x:state(x)[1],t,state(t)[1])
    cand=np.unique(np.r_[0.,end,tex]);temp=state(cand)[2]
    tp=float(cand[np.argmax(temp)]);peak=float(np.max(temp))
    rc=np.unique(np.r_[0.,end,rex]);rval=state(rc)[0];tr=float(rc[np.argmin(rval)])
    minima=[float(x) for x in rex if x>0 and state(max(0,x-1e-12))[1]<0 and state(min(end,x+1e-12))[1]>0]
    strong=[x for x in minima if state(x)[0]<m.p.R0/2]
    first=[x for x in minima if x<=1/m.p.f]
    tc=strong[0] if strong else (min(first,key=lambda x:state(x)[0]) if first else tr)
    return dict(observed_Tmax_K=peak,global_Tmax_K=peak if s.completed and end==s.end else np.nan,
        t_Tmax_s=tp,Rmin_m=float(rval.min()),t_Rmin_s=tr,tc_s=tc,
        selected_event='first_strong_minimum' if strong else ('cycle_1_deepest_minimum' if first else 'observed_interval_endpoint_reference'),
        first_strong_s=strong[0] if strong else np.nan,
        cycle_1_deepest_s=min(first,key=lambda x:state(x)[0]) if first else np.nan,
        temperature_stationary_times_s=tex,radius_stationary_times_s=rex,
        peak_candidates_s=cand,peak_candidates_K=temp)

def quadrature(s,order=8):
    m=s.model;gx,gw=leggauss(order)
    totals=np.zeros(6);pos=np.zeros(3);neg=np.zeros(3)
    # Independently evaluate physical channels on dense polynomials, not ledger differences.
    for a,b in zip(s.t[:-1],s.t[1:]):
        tt=(a+b)/2+(b-a)/2*gx
        for t,w in zip(tt,gw):
            on=abs(t-s.center)<s.duration/2
            pa=m.abs_per_fluence*s.fluence/s.duration if on else 0.
            inc=np.pi*m.particle.beam_radius_m**2*s.fluence/s.duration if on else 0.
            qg,qa,qb,work,*_=m.powers(t,s.sol(t)[:m.dim],pa,inc,s.path)
            weight=(b-a)/2*w
            totals+=weight*np.array([work,qg,qa,qb,pa,inc])
            flows=np.array([qg,qa,qb]);pos+=weight*np.maximum(flows,0);neg+=weight*np.minimum(flows,0)
    return totals,pos,neg

def evaluate(s,save_name=None,do_quad=True):
    m=s.model;e=extrema(s)
    t=np.unique(np.r_[s.t,np.linspace(0,s.end,1400),e['temperature_stationary_times_s'],e['radius_stationary_times_s']])
    y=s.sol(t);R,U,T,tl,ta=m.state(t,y[:m.dim],s.path)
    rnp=center_radius(R,m.particle.radius_m,m.particle.gap_m,m.particle.motion)
    gap=rnp-R-m.particle.radius_m
    powers=np.array([m.powers(x,y[:m.dim,j],path=s.path)[:4] for j,x in enumerate(t)])
    qg,qa,qb,work=powers.T
    W,Qg,Qa,Qb,Eabs,Einc=y[m.ledger:m.ledger+6]*m.mix.Eref
    Ug=m.mix.du_array(T)
    El=m.grid.capacity@tl;Ea=m.particle.count*m.particle.capacity*ta
    rg=Ug+W-Qg;rl=El-Qa+Qg+Qb;ra=Ea-Eabs+Qa
    rt=Ug+El+Ea+W+Qb-Eabs
    scale=max(m.mix.Eref,float(np.max(np.abs(W))),float(np.max(np.abs(Eabs))))
    H=4*np.pi*m.cfg.eta*m.mix.kg*R
    tw=T+np.divide(qg,H,out=np.zeros_like(qg),where=H!=0)
    tw=np.where(H!=0,tw,T0 if m.cfg.clamp_water else T0+tl[0])
    pp=m.mix.nt*8.31446261815324*T/(4*np.pi*R**3/3)
    # Ring normal to beam, centered on bubble: no projected shadowing if separation>2a.
    sep=2*rnp*np.sin(np.pi/m.particle.count) if m.particle.count>1 else np.full_like(R,np.inf)
    optics_ok=np.min(gap)>=5*m.particle.radius_m*(1-1e-10)
    reasons=[]
    if m.cfg.amplitude_bar==1.3:reasons.append('1.3 bar mandated formal; strong compressibility')
    if np.max(np.abs(U))/m.p.c>.1:reasons.append('wall Mach above 0.1 diagnostic screen')
    if m.cfg.amplitude_bar*1e5>m.p.p0:reasons.append('negative far-field absolute pressure during forcing')
    if m.particle.count and m.particle.motion!='passive':reasons.append('prescribed confinement or fixed-center force unverified')
    if not optics_ok and m.particle.count:reasons.append('homogeneous isolated-host near-wall screen fails')
    if np.max(T)>2000:reasons.append('fixed composition and constant gas transport at high temperature')
    if np.max(ta)>30:reasons.append('Au exceeds moderate-rise constant-property screen')
    if np.max(tl)+T0>373.15:reasons.append('water exceeds bulk boiling reference; no phase-change model')
    if s.status!='completed':reasons.append(s.status)
    reasons.append('spherical source averaging; local hotspot unresolved')
    row={k:v for k,v in e.items() if not isinstance(v,np.ndarray)}
    row.update(Pa_bar=m.cfg.amplitude_bar,gas=m.cfg.gas,vapor_factor=m.cfg.vapor_factor,
        completed=s.completed,status=s.status,end_s=s.end,requested_end_s=s.requested_end,
        boundary_T_K=float(T[-1]),N=m.particle.count,motion=m.particle.motion,
        a_m=m.particle.radius_m,d0_m=m.particle.gap_m,fluence_J_m2=s.fluence,
        pulse_center_s=s.center,duration_s=s.duration,wavelength_m=m.particle.wavelength_m,
        max_wall_Mach=float(np.max(np.abs(U))/m.p.c),min_gap_m=float(np.min(gap)),
        min_ring_separation_over_2a=float(np.min(sep)/(2*m.particle.radius_m)),
        max_gold_K=float(T0+np.max(ta)),max_spherical_water_K=float(T0+np.max(tl)),
        gas_domain_K=str((m.mix.Tmin,m.mix.Tmax)),optical_host_screen=bool(optics_ok),
        pulse_equilibration_screen=bool(s.duration>=10*max(1.7e-12,m.particle.radius_m**2/(317/(19320*129)))),
        gold_moderate_screen=bool(np.max(ta)<=30),water_373K_reference_screen=bool(T0+np.max(tl)<=373.15),
        Au_bulk_melting_reference_screen=bool(T0+np.max(ta)<1337.33),contact_screen=bool(np.min(gap)>-1e-18),
        photon_allowance_ratio=m.particle.metadata()['N_sigma_ext_over_beam_area'],
        thermal_energy_scale_J=scale,gas_residual_scaled=float(np.max(abs(rg))/scale),
        water_residual_scaled=float(np.max(abs(rl))/scale),Au_residual_scaled=float(np.max(abs(ra))/scale),
        total_residual_scaled=float(np.max(abs(rt))/scale),nfev=s.nfev,accepted_steps=len(s.t)-1,
        E_incident_J=float(Einc[-1]),E_absorbed_J=float(Eabs[-1]),E_Au_J=float(Ea[-1]),E_liquid_J=float(El[-1]),
        delta_U_gas_J=float(Ug[-1]),W_by_gas_J=float(W[-1]),Q_l_to_g_J=float(Qg[-1]),
        Q_Au_to_l_J=float(Qa[-1]),Q_bath_J=float(Qb[-1]),
        classification='formal_extreme' if m.cfg.amplitude_bar==1.3 or np.max(T)>2000 else
             ('conditional_prescribed_geometry' if m.particle.count and m.particle.motion!='passive' else 'conditional_reduced_model'),
        limitations='; '.join(reasons),tangent=s.tangent,cells=m.cfg.cells,padding_m=m.cfg.padding_m,
        method=m.cfg.method,rtol=m.cfg.rtol,max_step_s=m.cfg.max_step_s)
    fields=dict(t_s=t,R_m=R,U_m_s=U,T_g_K=T,T_l_rise_K=tl,T_Au_rise_K=ta,
        p_g_Pa=pp,T_wall_K=tw,r_np_m=rnp,gap_m=gap,Q_l_to_g_W=qg,Q_Au_to_l_W=qa,
        Q_bath_W=qb,work_by_gas_W=work,U_gas_change_J=Ug,E_water_J=El,E_Au_J=Ea,
        work_by_gas_J=W,Q_l_to_g_J=Qg,Q_Au_to_l_J=Qa,Q_bath_J=Qb,E_absorbed_J=Eabs,
        E_incident_J=Einc,gas_residual_J=rg,water_residual_J=rl,Au_residual_J=ra,total_residual_J=rt,
        xi_faces_m3=m.grid.faces,xi_centers_m3=m.grid.centers,capacity_water_J_K=m.grid.capacity,
        accepted_t_s=s.t,accepted_y=s.sol(s.t),temperature_stationary_times_s=e['temperature_stationary_times_s'],
        radius_stationary_times_s=e['radius_stationary_times_s'],peak_candidates_s=e['peak_candidates_s'],
        peak_candidates_K=e['peak_candidates_K'])
    if s.tangent:
        delta=y[m.dim:];dt=delta[2]*T0
        fields.update(dT_dF_K_m2_J=dt,dR_dF_m3_J=delta[0]*m.p.R0,dU_dF=delta[1]*m.p.c,
            dTl_dF=delta[3:m.au]*T0,dTAu_dF=delta[m.au]*T0,
            dledgers_dF_J=delta[m.ledger:m.ledger+6]*m.mix.Eref)
        yp=s.sol(e['t_Tmax_s']);row['dTmax_dF_K_m2_J']=float(yp[m.dim+2]*T0)
        deltaUg=m.mix.nt*m.mix.cv_array(T)*dt
        deltaEl=m.grid.capacity@(delta[3:m.au]*T0)
        deltaEa=m.particle.count*m.particle.capacity*(delta[m.au]*T0)
        dw,dqg,dqa,dqb,deabs,deinc=delta[m.ledger:m.ledger+6]*m.mix.Eref
        fields['incremental_total_residual_J_per_F']=deltaUg+deltaEl+deltaEa+dw+dqb-deabs
        row['incremental_balance_over_unit_absorbed']=float(np.max(abs(fields['incremental_total_residual_J_per_F']))/max(m.abs_per_fluence,1e-30))
    if do_quad:
        quad,pos,neg=quadrature(s)
        for name,val,pv,nv in zip(('l_to_g','Au_to_l','bath'),quad[1:4],pos,neg):
            row['quadrature_'+name+'_J']=float(val);row[name+'_positive_J']=float(pv);row[name+'_negative_J']=float(nv)
        row['independent_total_residual_scaled']=float(abs(Ug[-1]+El[-1]+Ea[-1]+quad[0]+quad[3]-quad[4])/scale)
        row['quadrature_ledger_max_error_scaled']=float(np.max(abs(quad-np.array([W[-1],Qg[-1],Qa[-1],Qb[-1],Eabs[-1],Einc[-1]])))/max(scale,abs(quad[5])))
    if save_name:
        (ROOT/'results/raw').mkdir(parents=True,exist_ok=True)
        path=ROOT/'results/raw'/save_name
        np.savez_compressed(str(path)+'.npz',**fields)
        def clean(x):
            if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
            if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
            if isinstance(x,np.generic):return x.item()
            return x
        path.parent.joinpath(path.name+'.json').write_text(json.dumps(clean(dict(summary=row,config=asdict(m.cfg),particle=m.particle.metadata())),indent=2)+'\n')
    return row,fields

def pair(base,heated):
    a,b=extrema(base),extrema(heated);end=min(base.end,heated.end)
    aa,bb=extrema(base,end),extrema(heated,end)
    completed=base.completed and heated.completed and base.requested_end==heated.requested_end
    return dict(completed_pair=completed,comparison_end_s=end,
        baseline_observed_Tmax_K=a['observed_Tmax_K'],laser_observed_Tmax_K=b['observed_Tmax_K'],
        baseline_t_Tmax_s=a['t_Tmax_s'],laser_t_Tmax_s=b['t_Tmax_s'],
        global_delta_Tmax_K=b['global_Tmax_K']-a['global_Tmax_K'] if completed else np.nan,
        common_prefix_delta_Tmax_K=bb['observed_Tmax_K']-aa['observed_Tmax_K'],
        event_shift_s=b['tc_s']-a['tc_s'],pulse_offset_s=heated.center-a['tc_s'],
        scope='five_cycle_global' if completed else 'common_prefix_only_unknown_tail')
