"""Bounded MR6 correction checks on saved data; no new sweep or coupled model."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.special import erfcx
from .upstream import ROOT,optics
from .coupling import Geometry,FrozenHeatModel,ThermalHistory

R=ROOT/'results'
ORIGINAL='faa2828f554517a0e2e43d9ae42cff37d682fbd0'

def table(name,rows):
    pd.DataFrame(rows).to_csv(R/(name+'.csv'),index=False,float_format='%.17g')

def array_hash(a):
    a=np.asarray(a)
    return hashlib.sha256(str(a.dtype).encode()+str(a.shape).encode()+a.tobytes(order='C')).hexdigest()

def record_original_inputs(original_root):
    """Retain baseline byte hashes and old source-array hashes before renaming."""
    records=[];sources=[]
    paths=sorted((original_root/'results').glob('*.csv'))+sorted((original_root/'results').rglob('*.npz'))
    paths += [original_root/'results'/s for s in ['parameters.json','design_extrema.json','validation_summary.json','optical_contract.json','local_flux_metrics.json']]
    for p in paths:
        rel=str(p.relative_to(original_root))
        if p.name.endswith('_source.npz'):
            with np.load(p,allow_pickle=False) as z:
                sources.append(dict(path=rel,arrays={k:array_hash(z[k]) for k in
                    ['time_s','q_wall_laser_increment_W','q_wall_baseline_W','q_wall_total_W']}))
        else:
            records.append(dict(path=rel,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    data=dict(original_MR6_commit=ORIGINAL,unchanged_files=records,renamed_source_arrays=sources,
              scope='Only source field names/metadata change; original numerical traces remain identical.')
    (R/'correction_input_manifest.json').write_text(json.dumps(data,indent=2)+'\n')

def migrate_sources():
    """Rename reference diagnostics, preserving all old numerical source values."""
    for p in sorted((R/'histories').glob('*_source.npz')):
        with np.load(p,allow_pickle=False) as z:
            if 'q_wall_reference_sum_W' in z and 'reference_bubble_radius_m' in z:continue
            tt=z['time_s'].copy();laser=z['q_wall_laser_increment_W'].copy()
            base=z['q_wall_reference_baseline_W' if 'q_wall_reference_baseline_W' in z else 'q_wall_baseline_W'].copy()
            meta=json.loads(str(z['metadata_json']))
        with np.load(p.with_name(p.name.replace('_source.npz','.npz')),allow_pickle=False) as z:
            if not np.array_equal(tt,z['time_s']):raise RuntimeError('Source/gas time mismatch')
            reference_T=z['T_K'].copy();reference_R=z['R_m'].copy()
        meta['baseline']='MR4 closure at the paired exported reference_gas_temperature_K; reference trace only'
        meta['reference_temperature_state']='paired T_K = Tb + first-order delta_T; not a supplied future gas state'
        ThermalHistory(tt,laser,base,meta,reference_T,reference_R).save(p)

def check_original_inputs():
    data=json.loads((R/'correction_input_manifest.json').read_text())
    for row in data['unchanged_files']:
        if hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()!=row['sha256']:
            raise RuntimeError('Original numeric result changed: '+row['path'])
    mapping={'time_s':'time_s','q_wall_laser_increment_W':'q_wall_laser_increment_W',
             'q_wall_baseline_W':'q_wall_reference_baseline_W','q_wall_total_W':'q_wall_reference_sum_W'}
    for row in data['renamed_source_arrays']:
        with np.load(ROOT/row['path'],allow_pickle=False) as z:
            if 'q_wall_total_W' in z:raise RuntimeError('Ambiguous total source field remains')
            for old,new in mapping.items():
                if array_hash(z[new])!=row['arrays'][old]:raise RuntimeError('Changed source values '+row['path'])
    return len(data['unchanged_files']),len(data['renamed_source_arrays'])

def audit_handoff_numbers():
    """Each rounded scientific handoff claim has its source and rounding allowance."""
    rows=[]
    def claim(name,actual,stated,tolerance,source,unit,relation='rounded'):
        rows.append(dict(claim=name,source=source,unit=unit,stated_value=stated,actual_value=float(actual),
                         absolute_tolerance=tolerance,relation=relation,
                         passed=bool(actual<=stated+tolerance if relation=='upper_bound' else abs(actual-stated)<=tolerance)))
    with np.load(R/'reference.npz') as z:
        E=z['absorbed_so_far_J'][-1];j=np.argmax(z['q_wall_laser_increment_W'])
        for args in [('reference_absorption',E*1e15,3.4337324474063743,1e-12,'fJ'),
                     ('Au_peak',max(z['gold_rise_K']),2.3066,5e-5,'K rise'),
                     ('facing_water_peak',max(z['nearest_water_rise_K']),1.6322,5e-5,'K rise'),
                     ('gas_power_peak_time',z['time_s'][j]*1e9,53.94,.005,'ns from start'),
                     ('gas_power_peak_lag',z['time_s'][j]*1e9-5,48.94,.005,'ns from center'),
                     ('gas_peak_power',z['q_wall_laser_increment_W'][j]*1e12,14.0702,5e-5,'pW'),
                     ('gas_energy_2us',z['gas_energy_J'][-1]*1e18,11.1865,5e-5,'aJ'),
                     ('gas_fraction_2us',z['gas_energy_J'][-1]/E*100,.3258,5e-5,'%')]:
            name,value,stated,tol,unit=args;claim(name,value,stated,tol,'reference.npz',unit)
    optical=float(optics.optical_cross_sections(20e-9,532e-9,1.33).absorption_m2)
    claim('sigma_abs',optical,3.4337324474063743e-15,1e-27,'MR5 optical API','m^2')
    local=json.loads((R/'local_flux_metrics.json').read_text())
    for name,value,tol in [('max_local_flux_W_m2',161.64,.005),('surface_mean_flux_W_m2',.03110,.000005),('max_over_mean',5197,.5)]:
        claim(name,local[name],value,tol,'local_flux_metrics.json','SI or dimensionless as named')
    claim('area90_percent',local['area90_fraction']*100,.0481,.00005,'local_flux_metrics.json','%')
    t=pd.read_csv(R/'timing_scan.csv',float_precision='round_trip');best=t.loc[t.efficiency_late_window.idxmax()]
    claim('late_window_fraction',best.efficiency_late_window*100,.00817,.000005,'timing_scan.csv','%')
    claim('best_timing_center',best.delta_t_s*1e9,-50.,1e-10,'timing_scan.csv','ns')
    design=pd.read_csv(R/'finite_design.csv');admit=design[design.moderate_static_candidate]
    claim('finite_design_count',len(design),216,0,'finite_design.csv','rows')
    claim('moderate_max_gas',admit.gas_energy_J.max()*1e15,3.9517,.00005,'finite_design.csv','fJ')
    gasrow=admit.loc[admit.gas_energy_J.idxmax()]
    claim('absorption_at_gas_max',gasrow.energy_absorbed_J*1e12,3.48462,.000005,'finite_design.csv','pJ')
    claim('moderate_max_absorption',admit.energy_absorbed_J.max()*1e12,3.48476,.000005,'finite_design.csv','pJ')
    claim('formal_all_grid_gas_max',design.gas_energy_J.max()*1e15,9.43494,.000005,'finite_design.csv','fJ')
    events=pd.read_csv(R/'mechanical_events.csv',float_precision='round_trip')
    for amp,tc,rc in [(1.,20.092624,4.39309),(1.15,15.664319,1.22816),(1.3,18.341064,.550809)]:
        row=events[events.Pa_bar==amp].iloc[0]
        claim(f'tc_{amp}',row.tc_s*1e6,tc,.0000005,'mechanical_events.csv','us')
        claim(f'Rc_{amp}',row.Rc_m*1e6,rc,.000005,'mechanical_events.csv','um')
    for amp,width in [(1.15,21.53),(1.3,2.278)]:
        row=events[events.Pa_bar==amp].iloc[0]
        claim(f'radius_width_{amp}',row.radius_below_2Rc_width_s*1e9,width,.005,'mechanical_events.csv','ns')
    claim('separate_MR4_thermal_width',events[events.Pa_bar==1.15].separate_MR4_coupled_thermal_width_s.iloc[0]*1e9,17.85,.005,'mechanical_events.csv','ns')
    gas=pd.read_csv(R/'gas_diagnostic.csv',float_precision='round_trip')
    for mechanism,stated,tol in [('A_preheat',1.887e-5,5e-9),('B_wall',4.680e-7,5e-11),('C_B_energy',5.602e-5,5e-9),('C_ideal_upper',.003506,.0000005),('B_after',0.,0.)]:
        row=gas[(gas.Pa_bar==1.15)&(gas.mechanism==mechanism)].iloc[0]
        claim('delta_T_at_peak_'+mechanism,row.delta_T_at_baseline_peak_K,stated,tol,'gas_diagnostic.csv','K')
    row=gas[(gas.Pa_bar==1.15)&(gas.mechanism=='B_wall')].iloc[0]
    claim('max_B_delta_T',row.max_delta_T_K,6.150e-6,5e-10,'gas_diagnostic.csv','K')
    for amp,stated in [(1.,16.),(1.15,119.),(1.3,193.)]:
        actual=gas[gas.Pa_bar==amp].max_pressure_mismatch.max()*100
        claim(f'baseline_pressure_mismatch_{amp}',actual,stated,.51,'gas_diagnostic.csv','%')
    row=gas[(gas.Pa_bar==1.3)&(gas.mechanism=='baseline')].iloc[0]
    claim('NASA_boundary_time',row.end_s*1e6,18.3408990706,5e-11,'gas_diagnostic.csv','us')
    claim('NASA_boundary_temperature',row.censoring_temperature_K,6000.,1e-5,'gas_diagnostic.csv','K')
    v=json.loads((R/'validation_summary.json').read_text())
    for key,stated,tol,factor in [('production_gas_energy_spatial_relative_max',.196,.0005,100),
        ('production_gold_spatial_relative_max',1.206,.0005,100),('production_water_spatial_relative_max',1.192,.0005,100),
        ('isolated_gold_peak_finest_relative_max',.0893,0.,100),('independent_Radau_gold_relative',1.41e-9,5e-12,1),
        ('independent_Radau_gas_relative',1.91e-12,5e-15,1),('energy_balance_relative_max',2.66e-14,5e-17,1),
        ('gaussian_normalization_relative_max',2.22e-16,5e-19,1),('baseline_static_matching_error',4.2146,.00005,100)]:
        claim(key,abs(v[key])*factor,stated,tol,'validation_summary.json','relative or % as named',
              relation='upper_bound' if key=='isolated_gold_peak_finest_relative_max' else 'rounded')
    claim('spatial_control_count',v['spatial_cases'],15,0,'validation_summary.json','rows')
    g=pd.read_csv(R/'gas_increment_convergence.csv')
    claim('increment_time_step_discrepancy',g.relative_difference.abs().max(),8.54e-7,5e-10,'gas_increment_convergence.csv','relative')
    point=pd.read_csv(R/'analytic_point_controls.csv');far=point[(point.a_m<30e-9)&(point.time_s>1e-6)].iloc[0]
    claim('far_point_2us_discrepancy',abs(far.point_relative_difference)*100,.052,.0005,'analytic_point_controls.csv','%')
    snapshots=pd.read_csv(R/'event_snapshots.csv')
    claim('formal_snapshot_count',sum(snapshots.status=='formal_frozen_motion_omitted'),45,0,'event_snapshots.csv','rows')
    feedback=pd.read_csv(R/'gas_feedback_estimate.csv',float_precision='round_trip')
    for amp,estimate,tol,bound,btol in [(1.,.2323,.00005,74.3,.05),
        (1.15,.09284,.000005,12.1,.05),(1.3,.000225,.0000005,.127,.0005)]:
        f=feedback[(feedback.Pa_bar==amp)&(feedback.mechanism=='B_wall')].iloc[0]
        claim(f'B_feedback_estimate_{amp}',f.estimated_return_fraction*100,estimate,tol,'gas_feedback_estimate.csv','%')
        claim(f'B_feedback_bound_{amp}',f.return_bound_fraction*100,bound,btol,'gas_feedback_estimate.csv','%')
    table('handoff_numeric_audit',rows)
    failed=[r['claim'] for r in rows if not r['passed']]
    if failed:raise RuntimeError('Rounded handoff claims disagree: '+', '.join(failed))
    return len(rows)

def independent_field_energy(xy,tri,temperature,props):
    """Exact degree-2 triangle quadrature of 2*pi*rho*rho_l*cp*T."""
    points=xy[tri];u=points[:,1]-points[:,0];v=points[:,2]-points[:,0]
    area=np.abs(u[:,0]*v[:,1]-u[:,1]*v[:,0])/2
    energy=np.zeros(temperature.shape[1])
    for b in [(2/3,1/6,1/6),(1/6,2/3,1/6),(1/6,1/6,2/3)]:
        rq=points[:,:,0]@b
        Tq=sum(b[j]*temperature[tri[:,j]] for j in range(3))
        energy+=(2*np.pi/3)*(area*rq)@Tq*props.rho_water_kg_m3*props.cp_water_J_kgK
    return energy

def verify_saved_fields():
    rows=[];props=optics.ThermalProperties()
    for p in [R/'reference.npz',*sorted((R/'fields').glob('*.npz'))]:
        with np.load(p,allow_pickle=False) as z:
            E=float(z['absorbed_so_far_J'][-1])
            water=independent_field_energy(z['xy_m'],z['triangles'],z['water_rise_K'],props)
            q=float(z['h_W_m2K'])*(z['wall_node_area_m2']@z['water_rise_K'])
            row=dict(file=str(p.relative_to(R)),saved_nodes=len(z['xy_m']),
                     independently_integrated_water_residual_over_E=max(abs(water-z['water_energy_J']))/E,
                     wall_power_field_relative=max(abs(q-z['q_wall_laser_increment_W']))/max(z['q_wall_laser_increment_W']),
                     balance_residual_over_E=max(abs(z['balance_residual_J']))/E,
                     min_saved_water_rise_K=float(z['water_rise_K'].min()))
            if max(row['independently_integrated_water_residual_over_E'],row['balance_residual_over_E'])>1e-10 or row['wall_power_field_relative']>1e-10:
                raise RuntimeError('Saved field balance failure: '+str(p))
            rows.append(row)
    table('final_field_verification',rows)

def gas_boundary_step(time_s,R_m,h_W_m2K,k_W_mK,alpha_m2_s):
    """Uniform prescribed gas-temperature step -> surface-water response/K.

    Infinite exterior liquid, stationary sphere, no finite Au obstruction.
    Laplace response h/[h+k/R+k*sqrt(s/alpha)]. This is an analytic diagnostic
    evaluated on saved gas histories, not a newly coupled gas/liquid solution.
    """
    t=np.asarray(time_s,dtype=float)
    if np.any(t<0) or min(R_m,k_W_mK,alpha_m2_s)<=0 or h_W_m2K<0:raise ValueError('Invalid stationary step control')
    beta=1/R_m+h_W_m2K/k_W_mK
    return h_W_m2K/(h_W_m2K+k_W_mK/R_m)*(1-erfcx(beta*np.sqrt(alpha_m2_s*t)))

def feedback_from_saved():
    rows=[];alpha=.6/(1000*4187);nodes,weights=leggauss(4)
    for amp in [1.,1.15,1.3]:
        for mechanism in ['A_preheat','B_wall','C_ideal_upper']:
            p=R/f'histories/gas_{amp:.2f}_{mechanism}.npz'
            with np.load(p) as z:
                tt=z['time_s'];delta=z['delta_T_K'];radius=z['R_m'];source_E=float(z['laser_energy_J'][-1])
            history=ThermalHistory.load(p.with_name(p.stem+'_source.npz'))
            Rc=history.metadata['frozen_radius_m'];tc=history.metadata['collapse_time_s'];h=.0264/Rc;H=4*np.pi*Rc*Rc*h
            estimates=[]
            for n in [500,1000]:
                grid=np.unique(np.r_[np.linspace(tt[0],tt[-1],n+1),np.linspace(tc-30e-9,min(tc+30e-9,tt[-1]),n//5+1)])
                grid=grid[(grid>=tt[0])&(grid<=tt[-1])]
                T=PchipInterpolator(tt,delta)(grid);dt=np.diff(grid);dT=np.diff(T)
                surface=np.zeros_like(grid)
                for i in range(1,len(grid)):
                    # End-based form is a sum of nonnegative durations, avoiding
                    # cancellation at nearly coincident saved/added time nodes.
                    age=(grid[i]-grid[1:i+1])[:,None]+dt[:i,None]*(1-nodes)/2
                    S=gas_boundary_step(age,Rc,h,.6,alpha)
                    surface[i]=np.sum(dT[:i]*(S@weights)/2)
                estimates.append((float(H*np.trapezoid(surface,grid)),float(max(surface)),
                                  float(H*np.trapezoid(np.maximum.accumulate(np.maximum(T,0)),grid))))
            returned,peak,bound=estimates[-1]
            already_counted=float(np.trapezoid(4*np.pi*.0264*radius*delta,tt))
            rows.append(dict(Pa_bar=amp,mechanism=mechanism,observed_start_s=tt[0],observed_end_s=tt[-1],
                frozen_R_m=Rc,h_W_m2K=h,source_energy_observed_J=source_E,
                already_counted_MR4_incremental_gas_loss_J=already_counted,
                estimated_omitted_liquid_return_J=returned,estimated_return_fraction=returned/source_E,
                conservative_prescribed_trace_return_bound_J=bound,return_bound_fraction=bound/source_E,
                max_extra_surface_water_K=peak,
                refinement_relative_difference=abs(returned-estimates[0][0])/max(abs(returned),1e-100),
                scope='prescribed saved deltaT; stationary spherical liquid; no finite Au; no gas update'))
    table('gas_feedback_estimate',rows)

def bounded_replay():
    """Three old geometries, same archived mesh; no new model/parameter scan."""
    rows=[]
    for p in [R/'reference.npz',R/'fields/event_1.00.npz',R/'fields/event_1.15.npz']:
        model=FrozenHeatModel.from_saved_field(p,Geometry(gap_m=100e-9) if p.name=='reference.npz' else None)
        with np.load(p) as z:
            o=model.pulse(float(z['absorbed_so_far_J'][-1]),10e-9,z['time_s'])
            row=dict(file=str(p.relative_to(R)),archived_mesh_nodes=len(model.xy),temporal_relative_tolerance=2e-7,
                admitted_maxnorm_relative_difference=2e-5)
            for key in ['gold_rise_K','nearest_water_rise_K','q_wall_laser_increment_W','gas_energy_J','water_energy_J']:
                row[key+'_maxnorm_relative']=float(max(abs(o[key]-z[key]))/max(abs(z[key])))
                if row[key+'_maxnorm_relative']>2e-5:raise RuntimeError('Saved-mesh replay differs '+key)
            rows.append(row)
    table('saved_mesh_replay',rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-original',type=Path)
    parser.add_argument('--replay',action='store_true',help='repeat only three archived-mesh thermal controls')
    args=parser.parse_args()
    if args.record_original:record_original_inputs(args.record_original)
    migrate_sources();preserved,sources=check_original_inputs()
    verify_saved_fields();feedback_from_saved();count=audit_handoff_numbers()
    if args.replay:bounded_replay()
    print(f'{preserved} original numeric files unchanged; {sources} source traces renamed without numerical change; {count} numeric claims checked')

if __name__=='__main__':main()
