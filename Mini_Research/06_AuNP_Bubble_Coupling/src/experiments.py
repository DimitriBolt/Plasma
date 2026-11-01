"""Rebuild finite MR6 experiments; no MR1--MR5 file is written."""
from dataclasses import asdict,replace
from pathlib import Path
import json,sys,hashlib,subprocess,platform
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.optimize import brentq,minimize_scalar
from .upstream import ROOT,REPO,BASE,optics,thermo,saved,Parameters
from .coupling import (Geometry,FrozenHeatModel,ThermalHistory,particle_trajectory,
                       gaussian_source,integrated_parts,MODEL_ID)

RESULT=ROOT/'results'
def write_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def csv(name,rows):
    pd.DataFrame(rows).to_csv(RESULT/(name+'.csv'),index=False,float_format='%.17g')

def audit_inputs():
    """Git blob comparison plus SHA-256 of every tracked pre-MR6 input."""
    inv=RESULT/'input_inventory.json'
    if (REPO/'.git').exists():
        records=subprocess.check_output(['git','ls-tree','-rz','--full-tree',BASE],cwd=REPO).split(b'\0')
        files=[]
        for item in filter(None,records):
            head,name=item.split(b'\t',1);rel=name.decode();p=REPO/rel
            expected=subprocess.check_output(['git','cat-file','blob',head.split()[2].decode()],cwd=REPO)
            actual=p.read_bytes()
            if actual!=expected: raise RuntimeError('Accepted upstream bytes changed: '+rel)
            files.append(dict(path=rel,bytes=len(actual),sha256=hashlib.sha256(actual).hexdigest()))
        write_json(inv,dict(accepted_snapshot=BASE,all_match=True,files=files))
    else:
        data=json.loads(inv.read_text())
        for rec in data['files']:
            if hashlib.sha256((REPO/rec['path']).read_bytes()).hexdigest()!=rec['sha256']:
                raise RuntimeError('Exported input hash mismatch: '+rec['path'])

def times(duration=10e-9,end=2e-6):
    return np.unique(np.r_[0.,np.geomspace(min(1e-13,duration/1000),end,620),
                           np.linspace(0,duration,50),duration+np.geomspace(max(1e-13,duration*1e-5),end-duration,120)])

def response(geometry=Geometry(gap_m=100e-9),duration=10e-9,h=None,props=None,ref=1.,end=2e-6,field=False,**kwargs):
    if h is None: h=.0264/geometry.radius_m
    model=FrozenHeatModel(geometry,h,props,refinement=ref,**kwargs)
    out=model.pulse(1e-15,duration,times(duration,end),save_field=field)
    return model,out

def summary(model,out,energy=1e-15):
    t=out['time_s'];peak=int(np.argmax(out['q_wall_laser_increment_W']))
    row=dict(a_m=model.geometry.particle_radius_m,d_m=model.geometry.gap_m,R_m=model.geometry.radius_m,
             h_W_m2K=model.h,G_W_m2K=model.G,k_W_mK=model.props.k_water_W_mK,
             energy_absorbed_J=energy,gas_energy_J=out['gas_energy_J'][-1],
             gas_fraction=out['gas_energy_J'][-1]/energy,
             max_gold_rise_K=float(max(out['gold_rise_K'])),
             max_water_rise_K=float(max(out['nearest_water_rise_K'])),
             t_Au_peak_s=t[np.argmax(out['gold_rise_K'])],t_water_peak_s=t[np.argmax(out['nearest_water_rise_K'])],
             t_q_peak_s=t[peak] if np.max(out['q_wall_laser_increment_W'])>0 else np.nan,
             q_peak_bracket_width_s=t[min(peak+1,len(t)-1)]-t[max(0,peak-1)],
             q_peak_at_window_end=bool(peak==len(t)-1),q_peak_scope='observed_finite_response_window',
             max_q_W=out['q_wall_laser_increment_W'][peak],
             max_balance_relative=max(abs(out['balance_residual_J']))/energy,
             minimum_water_rise_K=min(out['minimum_water_rise_K']),nodes=len(model.xy),
             continuum_gap_screen=bool(model.geometry.gap_m>=2e-9),
             homogeneous_water_optics_screen=bool(model.geometry.gap_m>=5*model.geometry.particle_radius_m))
    return row

def grid_and_sensitivities():
    rows=[];base=Geometry(gap_m=100e-9);cache={}
    for a in [20e-9,40e-9,80e-9,100e-9]:
        for ratio in [1.,5.,20.]:
            g=replace(base,particle_radius_m=a,gap_m=ratio*a)
            for duration in [1e-9,10e-9,100e-9]:
                model,out=response(g,duration)
                s=summary(model,out)
                for wavelength in [450e-9,506e-9,526e-9,532e-9,600e-9,800e-9]:
                    section=optics.optical_cross_sections(a,wavelength,1.33)
                    F=min(1000.,1e10*duration,30/s['max_gold_rise_K']*1e-15/float(section.absorption_m2))
                    E=float(optics.pulse_energy(F,section.absorption_m2));factor=E/1e-15
                    moderate=duration>=10*max(optics.thermal_scales(a)['gold_diffusion_s'],1.7e-12) and ratio>=5
                    rows.append({**s,'duration_s':duration,'gap_over_a':ratio,'wavelength_m':wavelength,
                                 'sigma_abs_m2':float(section.absorption_m2),'fluence_J_m2':F,'intensity_W_m2':F/duration,
                                 'incident_local_beam_J':F*np.pi*(25e-6)**2,
                                 'energy_absorbed_J':E,'gas_energy_J':s['gas_energy_J']*factor,
                                 'max_gold_rise_K':s['max_gold_rise_K']*factor,
                                 'max_water_rise_K':s['max_water_rise_K']*factor,'moderate_static_candidate':moderate,
                                 'status':'conditional_static_design' if moderate else 'near_boundary_optics_unvalidated'})
                if a==20e-9 and ratio==5 and duration==10e-9:
                    cross=optics.optical_cross_sections(a,532e-9,1.33);E=float(cross.absorption_m2)
                    # Store a physical F=1 J/m^2 reference and the entire field.
                    model,physical=response(g,duration,field=True)
                    for k in physical:
                        if k not in ('time_s','xy_m','triangles','boundary_edges','boundary_markers','wall_node_area_m2','h_W_m2K'):
                            physical[k]=physical[k]*(E/1e-15)
                    np.savez_compressed(RESULT/'reference.npz',**physical)
    csv('finite_design',rows)
    gaps=[]
    for d in [1e-9,2e-9,5e-9,10e-9,20e-9,50e-9,100e-9,200e-9,500e-9,1000e-9]:
        model,out=response(replace(base,gap_m=d));gaps.append(summary(model,out))
    csv('gap_scan',gaps)
    sensitivities=[]
    for param,vals in [('h',[0.,1e3,4.4e3,1e4,1e5,1e6,1e8,1e10]),('k',[.3,.6,1.2]),('G',[0.,5e7,1e8,2e8,1e12])]:
        for v in vals:
            kw={};props=None
            if param=='h':kw['h']=v
            if param=='G':kw['interface_W_m2K']=v
            if param=='k':props=replace(optics.ThermalProperties(),k_water_W_mK=v)
            model,out=response(base,props=props,**kw)
            sensitivities.append({**summary(model,out),'varied_parameter':param,'value_SI':v})
    csv('transport_sensitivity',sensitivities)
    formal=[]
    for duration in [1e-13,1e-9,100e-9]:
        model,out=response(replace(base,gap_m=1e-9),duration)
        s=summary(model,out)
        for F in [0.,.1,1.,10.,100.,1000.]:
            E=float(optics.pulse_energy(F,optics.optical_cross_sections(20e-9,526e-9,1.33).absorption_m2));f=E/1e-15
            flags=['near_interface_optics','sub_2nm_continuum','unresolved_capture']
            if duration<10*max(1.7e-12,optics.thermal_scales(20e-9)['gold_diffusion_s']): flags.append('one_temperature_failure')
            if s['max_gold_rise_K']*f>30:flags.append('moderate_screen_exceeded')
            if s['max_gold_rise_K']*f+293.15>1337.33:flags.append('bulk_melting_diagnostic')
            if s['max_water_rise_K']*f+293.15>373.15:flags.append('atmospheric_boiling_diagnostic')
            formal.append(dict(a_m=20e-9,d_m=1e-9,duration_s=duration,fluence_J_m2=F,absorbed_J=E,
                               gas_J=s['gas_energy_J']*f,formal_Au_rise_K=s['max_gold_rise_K']*f,
                               formal_water_rise_K=s['max_water_rise_K']*f,flags=';'.join(flags)))
    csv('formal_regimes',formal)
    # Dilute, nonshadowing equatorial ring; superposition is conditional on separation.
    multi=[];alpha=.6/(1000*4187);R=base.radius_m;rc=R+base.particle_radius_m+base.gap_m
    E=float(optics.optical_cross_sections(20e-9,532e-9).absorption_m2)
    ref=np.load(RESULT/'reference.npz');single=ref['gas_energy_J'][-1]
    ext=float(optics.optical_cross_sections(20e-9,532e-9).extinction_m2);beam=np.pi*(25e-6)**2
    for N in [0,1,2,4,8]:
        sep=2*rc*np.sin(np.pi/N) if N>1 else np.inf
        overlap=0 if N<=1 else (N-1)*np.exp(-(sep-40e-9)**2/(4*alpha*2e-6))
        multi.append(dict(N_particles=N,local_absorbed_J=N*E,incident_beam_J=beam,
                          beam_extinction_fraction=N*ext/beam,gas_energy_J=N*single,
                          center_separation_m=sep if np.isfinite(sep) else np.nan,
                          nearest_overlap_exponential_bound=overlap,
                          nonoverlap_screen=bool(sep>max(200e-9,40e-9+4*np.sqrt(alpha*2e-6))),
                          ring_envelope_number_density_m3=N/(4*np.pi*rc**2*40e-9),
                          status='conditional_separate_patches'))
    csv('multiparticle',multi)
    data=pd.DataFrame(rows);admitted=data[data.moderate_static_candidate]
    def winner(frame,column):
        top=frame[column].max()
        tied=frame[frame[column]>=top*(1-1e-12)]
        return tied.sort_values('sigma_abs_m2',ascending=False).iloc[0].to_dict()
    write_json(RESULT/'design_extrema.json',dict(
        moderate_max_absorbed=winner(admitted,'energy_absorbed_J'),
        moderate_max_gas_energy=winner(admitted,'gas_energy_J'),
        moderate_max_fraction=winner(admitted,'gas_fraction'),
        all_grid_max_gas_energy=winner(data,'gas_energy_J')))

def timing_maps():
    ref=np.load(RESULT/'reference.npz');t=ref['time_s'];q=ref['q_wall_laser_increment_W']
    source=PchipInterpolator(t,q,extrapolate=False);E=float(optics.optical_cross_sections(20e-9,532e-9).absorption_m2)
    rows=[];window=(-10e-9,10e-9);duration=10e-9
    for offset in [-2000e-9,-1000e-9,-500e-9,-200e-9,-100e-9,-50e-9,-20e-9,-10e-9,0.,10e-9,20e-9,50e-9,100e-9,500e-9]:
        start=offset-duration/2
        def integral(a,b):
            a=max(a-start,0.);b=min(b-start,t[-1])
            return float(source.integrate(a,b)) if b>a else 0.
        rows.append(dict(delta_t_s=offset,pulse_start_relative_tc_s=start,
                         gas_by_tc_J=integral(-3e-6,0.),gas_late_window_J=integral(*window),
                         efficiency_by_tc=integral(-3e-6,0.)/E,efficiency_late_window=integral(*window)/E,
                         t_Au_peak_relative_tc_s=start+t[np.argmax(ref['gold_rise_K'])],
                         t_water_peak_relative_tc_s=start+t[np.argmax(ref['nearest_water_rise_K'])],
                         t_q_peak_relative_tc_s=start+t[np.argmax(q)],
                         Au_lag_from_center_s=t[np.argmax(ref['gold_rise_K'])]-duration/2,
                         water_lag_from_center_s=t[np.argmax(ref['nearest_water_rise_K'])]-duration/2,
                         q_lag_from_center_s=t[np.argmax(q)]-duration/2))
    csv('timing_scan',rows)
    # Gap x timing uses unit-energy shape recomputed for each gap (not an arbitrary lag).
    matrix=[]
    for d in [5e-9,20e-9,100e-9,500e-9,1000e-9]:
        model,o=response(Geometry(gap_m=d));f=PchipInterpolator(o['time_s'],o['q_wall_laser_increment_W'])
        for offset in [-1000e-9,-500e-9,-200e-9,-100e-9,-50e-9,-20e-9,0.,20e-9,100e-9]:
            start=offset-5e-9;a=max(-10e-9-start,0);b=min(10e-9-start,o['time_s'][-1])
            matrix.append(dict(d_m=d,delta_t_s=offset,late_window_fraction=float(f.integrate(a,b))/1e-15 if b>a else 0.))
    csv('gap_timing_map',matrix)

def mechanical_events():
    admission=pd.read_csv(ROOT.parent/'02_Compressible_Bubble_Dynamics/results/admission.csv')
    rows=[];kin=[];snap=[]
    for amp in [1.,1.15,1.3]:
        path=saved.SavedMechanics('KM',amp)
        arow=admission[(admission.model=='KM')&np.isclose(admission.Pa_bar,amp)].iloc[0]
        assert bool(arow.completed) and path.completed and path.end>=5/Parameters().f*(1-1e-14)
        ev=[e for e in path.events if e['kind']=='minimum' and e['time_s']<1/Parameters().f]
        strong=[e for e in ev if e['R_m']<Parameters().R0/2]
        event=min(strong,key=lambda e:e['time_s']) if strong else min(ev,key=lambda e:e['R_m'])
        tc=event['time_s'];Rc=event['R_m'];prev=max(e['time_s'] for e in path.events if e['kind']=='maximum' and e['time_s']<tc)
        t=np.unique(np.r_[np.linspace(prev,tc+200e-9,2001),path.t[(path.t>=prev)&(path.t<=tc+200e-9)],tc])
        RR,UU=path(t);half=2*Rc
        left=brentq(lambda x:float(path(x)[0])-half,prev,tc) if float(path(prev)[0])>half else np.nan
        # Right crossing is optional, and never fabricated beyond the saved interval.
        tr=np.linspace(tc,tc+2e-6,2000);rg=path(tr)[0]-half
        candidates=np.where(rg[1:]*rg[:-1]<0)[0]
        right=brentq(lambda x:float(path(x)[0])-half,tr[candidates[0]],tr[candidates[0]+1]) if len(candidates) else np.nan
        local=path.t[(path.t>tc-20e-9)&(path.t<tc+20e-9)];steps=np.diff(local)
        thermal_tab=pd.read_csv(ROOT.parent/'04_Temperature_Barrier/results/summary.csv')
        th=thermal_tab[(thermal_tab.Pa_bar==amp)&(thermal_tab.gas=='air')&(thermal_tab['mode']=='coupled')&(thermal_tab.config=='M2_sat')&(thermal_tab.eta==1)].iloc[0]
        rows.append(dict(Pa_bar=amp,tc_s=tc,Rc_m=Rc,preceding_maximum_s=prev,
                         event_rule='first_strong_R_less_R0_over_2' if strong else 'cycle1_deepest_no_strong_event',
                         radius_below_2Rc_width_s=right-left,local_sample_max_step_s=max(steps) if len(steps) else np.nan,
                         local_sample_min_step_s=min(steps) if len(steps) else np.nan,
                         max_wall_Mach=float(arow.max_wall_Mach),input_admission=arow.next_stage_use,
                         completed=True,requested_end_s=path.end,
                         separate_MR4_coupled_global_Tpeak_s=th.t_Tmax_s,
                         separate_MR4_coupled_thermal_width_s=th.temperature_half_height_width_s))
        for mode in ['fixed','passive','capture']:
            rr,dd,status=particle_trajectory(RR,20e-9,100e-9,mode)
            for i in range(0,len(t),max(1,len(t)//300)):
                kin.append(dict(Pa_bar=amp,mode=mode,t_s=t[i],t_minus_tc_s=t[i]-tc,R_m=RR[i],U_m_s=UU[i],
                                particle_center_m=rr[i],surface_gap_m=dd[i],status=status[i],
                                invariant_m3=rr[i]**3-RR[i]**3))
            for offset in [-200e-9,-20e-9,0.,20e-9,200e-9]:
                tl=tc+offset
                if tl<prev: continue
                Rl,Ul=map(float,path(tl));r_l,d_l,s_l=particle_trajectory(np.array([float(path(prev)[0]),Rl]),20e-9,100e-9,mode)
                prior_contact=bool(np.any(dd[t<=tl]<=0))
                if prior_contact or d_l[-1]<=0:
                    snap.append(dict(Pa_bar=amp,mode=mode,delta_t_s=offset,status='contact_stopped',R_frozen_m=Rl,d_m=d_l[-1]));continue
                g=Geometry(radius_m=Rl,gap_m=float(d_l[-1]));m,o=response(g)
                alpha=.6/(1000*4187);L=g.particle_radius_m+g.gap_m;td=L**2/alpha
                check=np.unique(np.r_[np.linspace(tl,min(tl+2e-6,path.end),300),tc,
                                      path.t[(path.t>=tl)&(path.t<=min(tl+2e-6,path.end))]])
                check=check[(check>=tl)&(check<=tl+2e-6)];rcheck,ucheck=path(check)
                Pe=max(abs(ucheck))*L/alpha
                change=max(abs(rcheck-Rl))/Rl
                rtraj, dtraj, _=particle_trajectory(np.r_[float(path(prev)[0]),rcheck],20e-9,100e-9,mode)
                gap_change=max(abs(dtraj[1:]-g.gap_m))/max(g.gap_m,20e-9)
                td_end=tl+td
                td_grid=np.unique(np.r_[np.linspace(tl,min(td_end,path.end),300),
                                        path.t[(path.t>=tl)&(path.t<=min(td_end,path.end))]])
                td_change=(max(abs(path(td_grid)[0]-Rl))/Rl if td_end<=path.end else np.nan)
                frozen_ok=Pe<.1 and change<.01 and gap_change<.1 and amp<1.3
                snap.append({**summary(m,o),'Pa_bar':amp,'mode':mode,'delta_t_s':offset,'R_frozen_m':Rl,
                             'diffusion_time_s':td,'Pe_max_response':Pe,'relative_radius_change_response':change,
                             'relative_radius_change_diffusion_time':td_change,'relative_gap_change_response':gap_change,
                             'status':'conditional_frozen' if frozen_ok else 'formal_frozen_motion_omitted'})
    csv('mechanical_events',rows);csv('kinematics',kin);csv('event_snapshots',snap)

def save_benchmarks():
    """Full finite-particle fields/time curves for multiple explicit geometries."""
    directory=RESULT/'fields';directory.mkdir(exist_ok=True)
    events=pd.read_csv(RESULT/'mechanical_events.csv')
    cases=[('near_contact',Geometry(particle_radius_m=100e-9,gap_m=1e-9),1e6,None),
           ('far_particle',Geometry(gap_m=500e-9),4400.,None)]
    for amp in [1.,1.15]:
        rec=events[events.Pa_bar==amp].iloc[0]
        cases.append((f'event_{amp:.2f}',Geometry(radius_m=rec.Rc_m,gap_m=100e-9),.0264/rec.Rc_m,rec.tc_s))
    keep={'time_s','xy_m','triangles','boundary_edges','boundary_markers','wall_node_area_m2','h_W_m2K'}
    for key,g,h,tc in cases:
        m,o=response(g,h=h,field=True)
        E=float(optics.optical_cross_sections(g.particle_radius_m,532e-9,1.33).absorption_m2)
        o={k:v if k in keep else v*(E/1e-15) for k,v in o.items()}
        meta=dict(geometry=asdict(g),duration_s=10e-9,fluence_J_m2=1.,wavelength_m=532e-9,
                  h_W_m2K=h,model_id=MODEL_ID,collapse_time_s=tc,
                  status='formal_frozen_event' if tc is not None else ('formal_near_boundary' if g.gap_m<5*g.particle_radius_m else 'conditional_static'))
        if tc is not None:o['absolute_time_s']=tc-25e-9+o['time_s']
        o['metadata_json']=np.array(json.dumps(meta,sort_keys=True))
        np.savez_compressed(directory/(key+'.npz'),**o)

def gas_diagnostic(amplitudes=(1.,1.15,1.3),prefix='',delta_max_step_ns=2.):
    """Common nonlinear baseline plus first-order laser sensitivity on fixed R(t).

    Integrating the tiny increment directly avoids subtracting two noisy large
    temperatures. This is the controlled small-increment branch, not a nonlinear
    or mechanically self-consistent laser experiment.
    """
    tab=pd.read_csv(ROOT.parent/'04_Temperature_Barrier/results/summary.csv');rows=[]
    for amp in amplitudes:
        selected=tab[(tab.Pa_bar==amp)&(tab.gas=='air')&(tab['mode']=='one_way')&(tab.cv_mode=='variable')&(tab.vapor_factor==1)&(tab.eta==1)]
        if len(selected)!=1: raise RuntimeError('Ambiguous published MR4 baseline')
        rec=selected.iloc[0];path=saved.SavedMechanics('KM',amp);mix=thermo.Mixture()
        event=pd.read_csv(RESULT/'mechanical_events.csv').query('Pa_bar==@amp').iloc[0];tc=event.tc_s
        with np.load(ROOT.parent/f'04_Temperature_Barrier/results/trajectories/{rec.run_id}.npz') as z:
            tbase=z['t_s'];Tbase=z['T_K']
        start=tc-3e-6;end=tc+2e-6
        if start>tbase[-1]: raise RuntimeError('Published baseline cannot initialize window')
        Tinitial=float(PchipInterpolator(tbase,Tbase)(start))
        def base_rhs(s,y):
            t=start+s*1e-9;R,U=map(float,path(t));T=y[0]
            work=mix.pressure(R,T)*4*np.pi*R*R*U;loss=mix.heat_loss(R,T,1.)
            return np.array([(-work-loss)/mix.capacity(T,trial=True),work/mix.Eref,loss/mix.Eref])*1e-9
        def upper_base(s,y):return mix.Tmax-y[0]
        def lower_base(s,y):return y[0]-mix.Tmin
        upper_base.terminal=lower_base.terminal=True;upper_base.direction=lower_base.direction=-1
        base=solve_ivp(base_rhs,(0,(end-start)/1e-9),[Tinitial,0.,0.],rtol=2e-11,
                       atol=[1e-9,1e-12,1e-12],method='DOP853',max_step=2.,dense_output=True,events=[upper_base,lower_base])
        if not base.success:raise RuntimeError(base.message)
        grid=np.unique(np.r_[np.linspace(0,base.t[-1],2500),base.t]);temperatures=base.sol(grid)[0];ib=np.argmax(temperatures)
        lo=grid[max(0,ib-1)];hi=grid[min(len(grid)-1,ib+1)]
        fit=minimize_scalar(lambda x:-float(base.sol(lo+(hi-lo)*x)[0]),bounds=(0,1),method='bounded',options={'xatol':1e-13})
        base_peak_s=lo+(hi-lo)*fit.x
        if base.status:base_peak_s=base.t[-1]
        g=Geometry(radius_m=float(path(tc)[0]),gap_m=100e-9)
        model,heat=response(g,end=5e-6);E=float(optics.optical_cross_sections(20e-9,532e-9).absorption_m2)
        qinterp=PchipInterpolator(heat['time_s'],heat['q_wall_laser_increment_W']*(E/1e-15))
        EB=heat['gas_energy_J'][-1]*E/1e-15
        for mechanism,offset in [('baseline',0.),('A_preheat',-2e-6),('B_wall',-20e-9),('B_after',200e-9),('C_B_energy',-5e-9),('C_ideal_upper',-5e-9)]:
            pulse_start=tc+offset-5e-9
            def source(t):
                if mechanism=='baseline':return 0.
                if mechanism.startswith('C_'):
                    return float(gaussian_source(t,EB if mechanism=='C_B_energy' else E,tc+offset,2e-9,(tc-20e-9,tc+10e-9)))
                age=t-pulse_start
                return float(qinterp(age)) if 0<=age<=heat['time_s'][-1] else 0.
            def rhs(s,y):
                t=start+s*1e-9;R,U=map(float,path(t));T=float(base.sol(s)[0]);delta=y[0]
                C=mix.capacity(T,trial=True);coef=mix._a(T)
                Cprime=mix.nt*thermo.RU*(-2*coef[0]/T**3-coef[1]/T**2+coef[3]+2*coef[4]*T+3*coef[5]*T*T+4*coef[6]*T**3)
                work_prime=3*mix.nt*thermo.RU*U/R;K=4*np.pi*mix.kg*R
                N=-work_prime*T-K*(T-thermo.T0);fprime=(-work_prime-K)/C-N*Cprime/C**2;q=source(t)
                return np.array([fprime*delta+q/C,work_prime*delta/mix.Eref,K*delta/mix.Eref,q/mix.Eref])*1e-9
            def upper(s,y):return mix.Tmax-float(base.sol(s)[0])-y[0]
            upper.terminal=True;upper.direction=-1
            sol=solve_ivp(rhs,(0,base.t[-1]),np.zeros(4),rtol=2e-10,atol=[1e-12,1e-18,1e-18,1e-18],
                          method='DOP853',max_step=delta_max_step_ns,dense_output=True,events=upper)
            if not sol.success:raise RuntimeError(sol.message)
            stop=start+sol.t[-1]*1e-9;censored=bool(base.status or sol.status)
            ss=np.unique(np.r_[np.linspace(0,sol.t[-1],2500),base.t[base.t<=sol.t[-1]],sol.t,
                               base_peak_s if base_peak_s<=sol.t[-1] else sol.t[-1]])
            tt=start+ss*1e-9;states=sol.sol(ss);baseline_states=base.sol(ss);Tb=baseline_states[0];delta=states[0];T=Tb+delta
            def minus_T(s):return -float(base.sol(s)[0]+sol.sol(s)[0])
            idx=np.argmax(T);lo=ss[max(0,idx-1)];hi=ss[min(len(ss)-1,idx+1)]
            fit=minimize_scalar(lambda x:minus_T(lo+(hi-lo)*x),bounds=(0,1),method='bounded',options={'xatol':1e-13})
            sp=lo+(hi-lo)*fit.x;tp=start+sp*1e-9;Tm=-fit.fun
            if censored:tp=stop;Tm=float(T[-1])
            R,U=path(tt);pg=np.array([mix.pressure(r,t) for r,t in zip(R,T)])
            old=Parameters().pg0*(Parameters().R0/R)**(3*Parameters().kappa)+Parameters().pv;mismatch=pg/old-1
            dU=mix.nt*mix.cv_array(Tb)*delta
            increment_residual=dU+(states[1]+states[2]-states[3])*mix.Eref
            residual=mix.du_array(Tb)-mix.du(Tinitial)+(baseline_states[1]+baseline_states[2])*mix.Eref+increment_residual
            baseline=-np.array([mix.heat_loss(r,t,1) for r,t in zip(R,T)])
            laser=np.array([source(t) for t in tt]);parts=integrated_parts(tt,baseline+laser)
            key=f'{prefix}gas_{amp:.2f}_{mechanism}'
            metadata=dict(model_id=MODEL_ID,time_unit='s',flux_unit='W',energy_unit='J',collapse_time_s=tc,
                          interval_s=[start,stop],frozen_radius_m=g.radius_m,gap_m=g.gap_m,particle_radius_m=20e-9,
                          fluence_J_m2=1.,pulse_center_s=tc+offset,baseline='signed MR4 closure; not absolute liquid baseline',
                          laser='zero gas-perturbation BVP; C is a diagnostic Gaussian',gas='first-order small-increment sensitivity',
                          status='formal_frozen_one_way_diagnostic',mechanism=mechanism)
            ThermalHistory(tt,laser,baseline,metadata).save(RESULT/'histories'/f'{key}_source.npz')
            np.savez_compressed(RESULT/'histories'/f'{key}.npz',time_s=tt,time_minus_tc_s=tt-tc,R_m=R,T_K=T,
                                baseline_T_K=Tb,delta_T_K=delta,p_Pa=pg,pressure_relative_mismatch=mismatch,
                                work_by_gas_J=(baseline_states[1]+states[1])*mix.Eref,
                                signed_loss_J=(baseline_states[2]+states[2])*mix.Eref,laser_energy_J=states[3]*mix.Eref,
                                energy_residual_J=residual,incremental_energy_residual_J=increment_residual)
            rows.append(dict(Pa_bar=amp,mechanism=mechanism,start_s=start,requested_end_s=end,end_s=stop,
                             completed_window=not censored,censoring_temperature_K=float(T[-1]) if censored else np.nan,
                             stop_reason='coefficient_boundary' if censored else 'window_complete',
                             observed_Tmax_K=Tm,t_Tmax_s=tp,tc_s=tc,Tmax_minus_tc_s=tp-tc,
                             max_delta_T_K=float(max(delta)),min_delta_T_K=float(min(delta)),
                             delta_T_at_baseline_peak_K=float(sol.sol(base_peak_s)[0]) if base_peak_s<=sol.t[-1] else np.nan,
                             laser_pressure_increment_relative_max=float(max(delta/Tb)),
                             local_absorbed_J=E,ideal_upper_J=E,computed_B_energy_J=EB,direct_source_energy_J=states[3,-1]*mix.Eref,
                             max_pressure_mismatch=float(max(abs(mismatch))),max_energy_residual_J=float(max(abs(residual))),
                             incremental_energy_residual_J=float(max(abs(increment_residual))),
                             published_baseline_run=rec.run_id,initial_temperature_K=Tinitial,**parts,
                             status='formal_frozen_first_order_one_way_diagnostic'))
    csv(prefix+'gas_diagnostic',rows)

def refine_gas():
    gas_diagnostic(amplitudes=(1.15,),prefix='refined_',delta_max_step_ns=.5)
    a=pd.read_csv(RESULT/'gas_diagnostic.csv');b=pd.read_csv(RESULT/'refined_gas_diagnostic.csv');rows=[]
    for _,r in b.iterrows():
        p=a[(a.Pa_bar==r.Pa_bar)&(a.mechanism==r.mechanism)].iloc[0]
        rows.append(dict(Pa_bar=r.Pa_bar,mechanism=r.mechanism,max_delta_T_K=p.max_delta_T_K,
                         refined_max_delta_T_K=r.max_delta_T_K,
                         relative_difference=(p.max_delta_T_K/r.max_delta_T_K-1 if r.max_delta_T_K else 0.),
                         stop_time_difference_s=p.end_s-r.end_s))
    csv('gas_increment_convergence',rows)

def main():
    audit_inputs()
    write_json(RESULT/'parameters.json',dict(model_id=MODEL_ID,thermal_properties=asdict(optics.ThermalProperties()),
        geometry=asdict(Geometry(gap_m=100e-9)),reference_duration_s=10e-9,reference_fluence_J_m2=1.,
        reference_wavelength_m=532e-9,refinement=1.,rtol=2e-7,reference_h_W_m2K=4400.,response_end_s=2e-6,
        finite_design={'a_nm':[20,40,80,100],'gap_over_a':[1,5,20],'duration_ns':[1,10,100],
                       'wavelength_nm':[450,506,526,532,600,800],'fluence_cap_J_m2':1000.,
                       'intensity_cap_W_m2':1e10,'gold_rise_design_cap_K':30.,'objective':'separate absorbed J, gas J at 2us, fraction'},
        timing_window_relative_tc_s=[-1e-8,1e-8],baseline_matching_relative_difference=-.0264/(.6+.0264)))
    write_json(RESULT/'environment.json',dict(python=sys.version,platform=platform.platform(),
        numpy=np.__version__,scipy=__import__('scipy').__version__,triangle=__import__('triangle').__version__))
    for name,fn in [('finite design',grid_and_sensitivities),('timing',timing_maps),('mechanics',mechanical_events),('fields',save_benchmarks),('gas',gas_diagnostic)]:
        print(name,flush=True);fn()
    print('done',flush=True)
if __name__=='__main__':
    import argparse
    argparse.ArgumentParser(description='Recompute the bounded MR6 experiments from accepted read-only inputs.').parse_args()
    main()
