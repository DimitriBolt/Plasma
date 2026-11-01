"""Saved-data curation, physical power replay and paired energy accounting."""
from dataclasses import fields
import json,hashlib
import numpy as np
import pandas as pd
from .upstream import ROOT,REPO,T0
from .model import Model,Config
from .particle import Particle
from .stable_pairs import archived
from .experiments import table

R=ROOT/'results'

def load_model(meta):
    return Model(Config(**meta['config']),Particle(**{k:v for k,v in meta['particle'].items() if k in {f.name for f in fields(Particle)}}))

def all_metadata():
    data={}
    stable={r['run_id']:r for r in pd.read_csv(R/'stable_pairs.csv',float_precision='round_trip').to_dict('records')} if (R/'stable_pairs.csv').exists() else {}
    for path in sorted((R/'raw').glob('*.json')):
        if path.stem.startswith('pilot_') or path.stem in ['N0_1','base_1','laser_F10_1','laser_F100_1']:continue
        meta=json.loads(path.read_text());s=meta['summary'];m=load_model(meta)
        if path.stem in stable:
            s.update({k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in stable[path.stem].items()})
            s['mechanical_mode']='prescribed' if s['level']==2 else 'coupled'
            meta['paired_formulation']='nonlinear directional Jacobian integral on checked C1 no-laser reference'
        elif path.stem.startswith('oneway_'):s['mechanical_mode']='prescribed'
        else:s['mechanical_mode']='coupled'
        s['run_id']=path.stem
        s['optical_host_screen']=s['min_gap_m']>=5*s['a_m']*(1-1e-10)
        s['continuum_gap_2nm_screen']=s['min_gap_m']>=2e-9*(1-1e-10)
        if s.get('cycle_1_deepest_s') is None and s.get('first_strong_s') is None:
            s['selected_event']='observed_interval_endpoint_reference'
        with np.load(path.with_suffix('.npz')) as z:
            V=4*np.pi*z['R_m']**3/3;rho=m.mix.nt*m.mix.mass_mol/V
            ip=int(np.argmax(z['T_g_K']));Tp=float(z['T_g_K'][ip]);Rp=float(z['R_m'][ip])
            s.update(max_pressure_Pa=float(z['p_g_Pa'].max()),max_gas_density_kg_m3=float(rho.max()),
                density_at_Tmax_kg_m3=float(rho[ip]),pressure_at_Tmax_Pa=float(z['p_g_Pa'][ip]),
                pressure_over_liquid_bulk_modulus=float(z['p_g_Pa'].max()/(m.p.rho*m.p.c**2)),
                gas_acoustic_crossing_at_peak_s=Rp/np.sqrt((1+8.31446261815324/m.mix.cv(Tp))*8.31446261815324*Tp/m.mix.mass_mol),
                beam_position_screen=bool(np.max(z['r_np_m'])+m.particle.radius_m<=m.particle.beam_radius_m))
        s['local_hotspot_resolved']=False if s['N'] else None
        s['retaining_mechanism_verified']=False if s['N'] and s['motion']!='passive' else None
        if s['optical_host_screen']:
            s['limitations']='; '.join(x for x in s['limitations'].split('; ') if x!='homogeneous isolated-host near-wall screen fails')
        if s['Pa_bar']==1.3:s['classification']='formal_extreme'
        elif not s['gold_moderate_screen'] or not s['continuum_gap_2nm_screen'] or not s['optical_host_screen']:
            s['classification']='formal_property_or_geometry_continuation'
        elif s['observed_Tmax_K']>2000:s['classification']='formal_fixed_composition_sensitivity'
        elif s['N'] and s['motion']!='passive':s['classification']='conditional_prescribed_geometry'
        else:s['classification']='conditional_reduced_model'
        meta['summary']=s
        path.write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
        data[path.stem]=meta
    if stable:
        # Keep the scientific pair table and per-run metadata consistent after
        # tolerance-safe geometric screens and censoring-label corrections.
        for name,row in stable.items():
            for key in ['optical_host_screen','continuum_gap_2nm_screen','classification','limitations','selected_event']:
                row[key]=data[name]['summary'][key]
        table('stable_pairs',list(stable.values()))
    return data


def replay_powers(data):
    rows=[]
    for name,meta in data.items():
        m=load_model(meta);p=R/'raw'/(name+'.npz')
        with np.load(p) as zz:
            keys=['t_s','R_m','U_m_s','T_g_K','T_l_rise_K','T_Au_rise_K','Q_l_to_g_W','Q_Au_to_l_W','Q_bath_W','work_by_gas_W']
            z={k:zz[k] for k in keys}
            fi=zz['field_indices'] if 'field_indices' in zz else np.arange(len(z['t_s']))
        positions=np.unique(np.r_[np.linspace(0,len(fi)-1,min(70,len(fi))).astype(int),
                          np.searchsorted(fi,np.argmax(z['T_Au_rise_K']))])
        positions=np.clip(positions,0,len(fi)-1);errors=[]
        for j in positions:
            i=fi[j];y=m.initial();y[:3]=[z['R_m'][i]/m.p.R0,z['U_m_s'][i]/m.p.c,z['T_g_K'][i]/T0]
            y[3:m.au]=z['T_l_rise_K'][:,j]/T0;y[m.au]=z['T_Au_rise_K'][i]/T0
            powers=m.powers(z['t_s'][i],y)[:4]
            errors.append([float(powers[k]-z[key][i]) for k,key in enumerate(['Q_l_to_g_W','Q_Au_to_l_W','Q_bath_W','work_by_gas_W'])])
        err=np.max(abs(np.array(errors)),axis=0)
        scales=np.array([np.max(abs(z[k])) for k in ['Q_l_to_g_W','Q_Au_to_l_W','Q_bath_W','work_by_gas_W']])
        relative=err/np.maximum(scales,1e-15)
        passed=bool(np.all((relative<1e-8)|(err<1e-17)))
        rows.append(dict(run_id=name,sampled_fields=len(positions),max_Qgas_error_W=err[0],max_QAu_error_W=err[1],
            max_Qbath_error_W=err[2],max_work_error_W=err[3],max_normalized_error=float(relative.max()),passed=passed))
    table('constitutive_replay',rows)
    if not all(r['passed'] for r in rows):raise AssertionError('A raw field differs from current physical kernel')


def paired_energy(data):
    rows=[];decomp=[]
    old_energy={r['run_id']:r for r in pd.read_csv(R/'paired_energy.csv',float_precision='round_trip').to_dict('records')} if (R/'paired_energy.csv').exists() else {}
    old_decomp={r['run_id']:r for r in pd.read_csv(R/'feedback_decomposition.csv',float_precision='round_trip').to_dict('records')} if (R/'feedback_decomposition.csv').exists() else {}
    for name,meta in data.items():
        if not name.startswith('stable_'):continue
        baseid=name[len('stable_'):].split('_F')[0]
        path=R/'raw'/(name+'.npz')
        with np.load(path) as zz:z={k:zz[k] for k in zz.files}
        if 'storage_version' in z:
            # Full scalar residual histories were computed before selecting exact
            # field snapshots. Recheck their recorded maxima; do not reconstruct
            # missing field times or silently interpolate energy.
            row=old_energy[name]
            E=max(abs(z['paired_delta_ledgers_J'][4,-1]),1e-30)
            for channel,column in [('gas','gas'),('water','water'),('Au','Au'),('total','total')]:
                value=float(np.max(abs(z['paired_'+channel+'_residual_J']))/E)
                if not np.isclose(value,row[column+'_over_absorbed'],rtol=1e-12,atol=1e-15):raise AssertionError('Stored energy residual changed')
            rows.append(row);decomp.append(old_decomp[name]);continue
        b,_=archived(baseid);m=b.model
        t=z['t_s'];yb=b.sol(t);Tb=yb[2]*T0
        dUg=m.mix.du_array(z['T_g_K'])-m.mix.du_array(Tb)
        dEl=m.grid.capacity@z['paired_delta_Tl_K']
        dEa=m.particle.count*m.particle.capacity*(z['T_Au_rise_K']-yb[m.au]*T0)
        dW,dQg,dQa,dQb,dE,dI=z['paired_delta_ledgers_J']
        rg=dUg+dW-dQg;rl=dEl-dQa+dQg+dQb;ra=dEa-dE+dQa;rt=rg+rl+ra
        for k,v in [('paired_gas_residual_J',rg),('paired_water_residual_J',rl),('paired_Au_residual_J',ra),('paired_total_residual_J',rt)]:z[k]=v
        np.savez_compressed(path,**z)
        E=max(abs(dE[-1]),1e-30)
        row=dict(run_id=name,gas_over_absorbed=np.max(abs(rg))/E,water_over_absorbed=np.max(abs(rl))/E,
             Au_over_absorbed=np.max(abs(ra))/E,total_over_absorbed=np.max(abs(rt))/E)
        rows.append(row)
        tp=meta['summary']['baseline_t_Tmax_s'] if 'baseline_t_Tmax_s' in meta['summary'] else data[baseid]['summary']['t_Tmax_s']
        j=np.searchsorted(t,tp);j=min(j,len(t)-1)
        decomp.append(dict(run_id=name,Pa_bar=m.cfg.amplitude_bar,t_reference_s=tp,
             dT_at_reference_K=z['paired_delta_T_K'][j],delta_Ugas_J=dUg[j],delta_work_J=dW[j],delta_Qgas_J=dQg[j],
             delta_absorbed_J=dE[j],peak_global_delta_Tmax_K=meta['summary'].get('global_delta_Tmax_K')))
    table('paired_energy',rows);table('feedback_decomposition',decomp)


def timing_first_cycle(data):
    df=pd.read_csv(R/'timing_sensitivity.csv',float_precision='round_trip')
    for j,row in df.iterrows():
        bid=row.baseline_id;bm=data[bid];z=np.load(R/'raw'/(bid+'.npz'))
        cand=z['peak_candidates_s'];temps=z['peak_candidates_K'];mask=cand<=1/30000
        tp=float(cand[mask][np.argmax(temps[mask])])
        src=np.load(R/'raw'/(row.run_id+'.npz'))
        df.loc[j,'first_cycle_t_Tmax_s']=tp
        df.loc[j,'first_cycle_dT_dF_K_m2_J']=np.interp(tp,src['t_s'],src['dT_dF_K_m2_J'])
        df.loc[j,'incremental_Qgas_at_first_peak_J_per_F']=np.interp(tp,src['t_s'],src['dledgers_dF_J'][1])
        df.loc[j,'first_cycle_complete']=bool(bm['summary']['end_s']>=1/30000)
        df.loc[j,'response_scope']='global_peak_derivative' if row.completed else 'formal_fixed_reference_endpoint_derivative'
        df.loc[j,'global_dTmax_dF_K_m2_J']=row.dTmax_dF_K_m2_J if row.completed else np.nan
    df.to_csv(R/'timing_sensitivity.csv',index=False,float_format='%.17g')


def refresh_tables(data):
    for filename in ['summary','sensitivity_runs']:
        old=pd.read_csv(R/(filename+'.csv'),float_precision='round_trip');rows=[]
        for run in old.run_id:
            rows.append(data[run]['summary'])
        table(filename,rows)
    df=pd.read_csv(R/'pairs.csv',float_precision='round_trip')
    for j,r in df.iterrows():
        a,b=data[r.baseline_id]['summary'],data[r.run_id]['summary']
        for key,source in [('baseline_observed_Tmax_K',a),('laser_observed_Tmax_K',b),('baseline_t_Tmax_s',a),('laser_t_Tmax_s',b)]:
            df.loc[j,key]=source['observed_Tmax_K'] if 'observed' in key else source['t_Tmax_s']
        complete=bool(a['completed'] and b['completed'] and a['requested_end_s']==b['requested_end_s'])
        df.loc[j,'completed_pair']=complete
        df.loc[j,'global_delta_Tmax_K']=b['global_Tmax_K']-a['global_Tmax_K'] if complete else np.nan
        if complete:df.loc[j,'common_prefix_delta_Tmax_K']=df.loc[j,'global_delta_Tmax_K']
    df.to_csv(R/'pairs.csv',index=False,float_format='%.17g')
    table('all_runs',[v['summary'] for v in data.values()])
    links={r['run_id']:r for r in df.to_dict('records')}
    if (R/'stable_pairs.csv').exists():
        for r in pd.read_csv(R/'stable_pairs.csv',float_precision='round_trip').to_dict('records'):
            r['baseline_id']=r['base_id'];links[r['run_id']]=r
    if (R/'independent_Radau.csv').exists():
        for r in pd.read_csv(R/'independent_Radau.csv',float_precision='round_trip').to_dict('records'):
            links[r['laser_id']]=r
    app=[]
    for name,v in data.items():
        if v['summary'].get('tangent',False):continue
        s=dict(v['summary']);link=links.get(name,{})
        for key in ['baseline_id','baseline_observed_Tmax_K','laser_observed_Tmax_K','baseline_t_Tmax_s',
                    'laser_t_Tmax_s','global_delta_Tmax_K','common_prefix_delta_Tmax_K','comparison_end_s',
                    'completed_pair','pulse_offset_s','event_shift_s','scope']:
            if key in link:s[key]=link[key]
        s['finite_effect_resolution']='cancellation_controlled' if name.startswith('stable_') else ('independently_resolved' if name=='independent_Radau_laser' else 'direct_difference_unresolved_at_weak_source_scale')
        if s['fluence_J_m2']==0:s['finite_effect_resolution']='no_laser_baseline'
        app.append(s)
    table('applicability',app)


def main():
    from .data_access import compact
    if compact():raise SystemExit('Full scientific curation requires results/raw and MR7_DATA_MODE=full. Use src.reference_verify for the compact shipped data.')
    data=all_metadata();replay_powers(data);paired_energy(data);timing_first_cycle(data);refresh_tables(data)
    inv=json.loads((R/'input_inventory.json').read_text())
    bad=[x['path'] for x in inv['files'] if hashlib.sha256((REPO/x['path']).read_bytes()).hexdigest()!=x['sha256']]
    if bad:raise AssertionError('Upstream inputs changed: '+str(bad))
    (R/'input_preservation.json').write_text(json.dumps(dict(verified_files=len(inv['files']),all_unchanged=True,initial_sha=inv['initial_sha']),indent=2)+'\n')
    print('Admitted',len(data),'dynamic archives; physical powers replayed; paired energy checked; upstream preserved.')

if __name__=='__main__':main()
