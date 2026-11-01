"""Independent limits, mesh/time/domain controls and read-only input sensitivity."""
from dataclasses import replace
import json
import numpy as np
import pandas as pd
from scipy.integrate import quad
from .upstream import ROOT,optics,saved,thermo,Parameters
from .coupling import (Geometry,FrozenHeatModel,point_robin_fraction,gaussian_source,
                       particle_trajectory,ThermalHistory)
from .experiments import response,times,summary,csv,write_json,RESULT,refine_gas

def main():
    mesh=[]
    for a,d,h in [(20e-9,100e-9,4400.),(100e-9,1e-9,1e6),(100e-9,500e-9,4400.)]:
        series=[]
        for ref in [.75,1.,1.5,2.,3.]:
            m,o=response(Geometry(particle_radius_m=a,gap_m=d),h=h,ref=ref)
            series.append({**summary(m,o),'refinement':ref,'au_area_relative_error':m.au_area/(4*np.pi*a*a)-1,
                           'bubble_area_relative_error':m.bubble_area/(4*np.pi*m.geometry.radius_m**2)-1})
        for row in series:
            for col in ['max_gold_rise_K','max_water_rise_K','gas_energy_J','max_q_W','t_q_peak_s']:
                row[col+'_relative_to_finest']=row[col]/series[-1][col]-1
        mesh.extend(series)
    csv('mesh_convergence',mesh)
    checks=[];g=Geometry(gap_m=100e-9)
    m=FrozenHeatModel(g);t=times(end=200e-9)
    tight=m.pulse(1e-15,10e-9,t,rtol=2e-10)
    for method,rtol,step in [('BDF',2e-7,np.inf),('BDF',2e-9,1e-9),('Radau',2e-9,2e-9)]:
        o=m.pulse(1e-15,10e-9,t,rtol=rtol,method=method,max_step_s=step)
        checks.append(dict(method=method,rtol=rtol,max_step_s=step,
                           gold_max_relative=np.max(abs(o['gold_rise_K']-tight['gold_rise_K']))/max(tight['gold_rise_K']),
                           wall_flux_max_relative=np.max(abs(o['q_wall_laser_increment_W']-tight['q_wall_laser_increment_W']))/max(tight['q_wall_laser_increment_W']),
                           gas_energy_relative=o['gas_energy_J'][-1]/tight['gas_energy_J'][-1]-1,
                           independently_trapezoidal_gas_relative=np.trapezoid(o['q_wall_laser_increment_W'],t)/o['gas_energy_J'][-1]-1,
                           independently_trapezoidal_bath_absolute_J=np.trapezoid(o['q_bath_W'],t)-o['bath_energy_J'][-1],
                           balance_scaled=max(abs(o['balance_residual_J']))/1e-15))
    csv('time_method_convergence',checks)
    domains=[]
    for pad in [.5e-6,1e-6,2e-6]:
        for bc in ['bath','insulated']:
            mm,oo=response(replace(g,outer_padding_m=pad),outer_bc=bc,ref=1.5)
            domains.append({**summary(mm,oo),'outer_padding_m':pad,'outer_bc':bc,
                            'bath_energy_J':oo['bath_energy_J'][-1]})
    csv('domain_convergence',domains)
    isolated=[]
    for a,duration in [(20e-9,10e-9),(100e-9,100e-9)]:
        for ref in [1.,2.,3.]:
            geom=Geometry(particle_radius_m=a,gap_m=5*a,outer_padding_m=3e-6)
            model,out=response(geom,duration,ref=ref,isolated=True)
            exact=optics.RadialHeatModel(a,cells=720).pulse(1e-15,1.,duration,out['time_s'])
            isolated.append(dict(a_m=a,duration_s=duration,refinement=ref,nodes=len(model.xy),
                                 gold_peak_relative=max(out['gold_rise_K'])/max(exact['gold_rise_K'])-1,
                                 gold_history_relative=max(abs(out['gold_rise_K']-exact['gold_rise_K']))/max(exact['gold_rise_K']),
                                 water_energy_relative=out['water_energy_J'][-1]/exact['water_energy_J'][-1]-1,
                                 Au_area_relative_error=model.au_area/(4*np.pi*a*a)-1))
    csv('isolated_MR5',isolated)
    analytic=[]
    for a,d in [(20e-9,500e-9),(100e-9,1e-9)]:
        geom=Geometry(particle_radius_m=a,gap_m=d,outer_padding_m=2e-6)
        mm,oo=response(geom,1e-9,h=1e8,ref=2.)
        t=oo['time_s'];alpha=.6/(1000*4187)
        ff=point_robin_fraction(t,geom.radius_m,geom.radius_m+a+d,alpha,.6,1e8)
        for target in [2e-9,20e-9,200e-9,2e-6]:
            j=np.argmin(abs(t-target));actual=oo['gas_energy_J'][j]/1e-15
            analytic.append(dict(a_m=a,d_m=d,time_s=t[j],finite_particle_fraction=actual,
                                 instantaneous_point_fraction=ff[j],point_relative_difference=(ff[j]/actual-1) if actual else np.nan,
                                 control='far_small_particle' if d/a>=10 else 'point_failure_near_finite_sphere'))
    csv('analytic_point_controls',analytic)
    planar=[]
    for d in [20e-9,100e-9,500e-9]:
        for plane in [False,True]:
            mm,oo=response(Geometry(gap_m=d),ref=2.,planar=plane)
            planar.append({**summary(mm,oo),'planar':plane})
    csv('planar_spherical_comparison',planar)
    interpolation=[]
    for amp in [1.,1.15,1.3]:
        full=saved.SavedMechanics('KM',amp)
        event=pd.read_csv(RESULT/'mechanical_events.csv');tc=float(event[event.Pa_bar==amp].tc_s.iloc[0])
        t=np.unique(np.r_[np.linspace(tc-200e-9,tc+200e-9,2001),full.t[(full.t>tc-200e-9)&(full.t<tc+200e-9)]])
        rr,uu=full(t)
        for stride,grid in [(2,'export'),(4,'export'),(1,'accepted')]:
            coarse=saved.SavedMechanics('KM',amp,stride=stride,grid=grid)
            r,u=coarse(t)
            interpolation.append(dict(Pa_bar=amp,stride=stride,grid=grid,
                                      max_radius_difference_over_R0=max(abs(r-rr))/6e-6,
                                      max_velocity_difference_over_peak=max(abs(u-uu))/max(abs(uu)),
                                      max_relative_h_difference=max(abs(rr/r-1))))
    csv('mechanics_interpolation',interpolation)
    # Optical reference uses published table, and the separately saved MR5 optimum.
    c=optics.optical_cross_sections(20e-9,532e-9,1.33)
    prior=pd.read_csv(ROOT.parent/'05_Plasmonic_Heating/results/bounded_single_particle.csv')
    r=prior[(np.isclose(prior.radius_m,100e-9,rtol=1e-12,atol=0)) & (np.isclose(prior.duration_s,100e-9,rtol=1e-12,atol=0))].iloc[0]
    # Field names are audited rather than assuming a rounded prompt value.
    optical=dict(reference_sigma_abs_m2=float(c.absorption_m2),reference_Eabs_J=float(optics.pulse_energy(1,c.absorption_m2)),
                 published_single_particle_row=r.to_dict())
    write_json(RESULT/'optical_contract.json',optical)
    integrals=[]
    for center,sigma,interval in [(0.,2e-9,(-3e-9,1e-9)),(1e-9,.1e-9,(0.,2e-9))]:
        a,b=interval;val=quad(lambda x:float(gaussian_source(a+(b-a)*x,3e-15,center,sigma,interval))*(b-a),0,1,epsabs=1e-27,epsrel=1e-11)[0]
        integrals.append(dict(center_s=center,sigma_s=sigma,interval_start_s=a,interval_end_s=b,target_J=3e-15,integral_J=val,relative_error=val/3e-15-1))
    csv('gaussian_normalization',integrals)
    write_json(RESULT/'validation_summary.json',dict(
        spatial_cases=len(mesh),time_method_cases=len(checks),domain_cases=len(domains),
        production_gas_energy_spatial_relative_max=max(abs(r['gas_energy_J_relative_to_finest']) for r in mesh if r['refinement']==1.),
        production_gold_spatial_relative_max=max(abs(r['max_gold_rise_K_relative_to_finest']) for r in mesh if r['refinement']==1.),
        production_water_spatial_relative_max=max(abs(r['max_water_rise_K_relative_to_finest']) for r in mesh if r['refinement']==1.),
        isolated_gold_peak_finest_relative_max=max(abs(r['gold_peak_relative']) for r in isolated if r['refinement']==3.),
        independent_Radau_gold_relative=checks[-1]['gold_max_relative'],
        independent_Radau_gas_relative=checks[-1]['gas_energy_relative'],
        energy_balance_relative_max=max(r['max_balance_relative'] for r in mesh),
        gaussian_normalization_relative_max=max(abs(r['relative_error']) for r in integrals),
        mechanics_input_radius_difference_over_R0_max=max(r['max_radius_difference_over_R0'] for r in interpolation),
        baseline_static_matching_error=.0264/(.6+.0264),
        numerical_scope='empirical selected-case discrepancies, not rigorous error bounds',
        physical_validation=False))
    refine_gas()
    print('validation complete',flush=True)

if __name__=='__main__': main()
