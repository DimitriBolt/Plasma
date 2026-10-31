"""Numerical, independent-library, analytic and finite-energy controls."""
import json
from dataclasses import replace
import numpy as np
import pandas as pd
import miepython
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar
from .plasmonics import (GoldData,optical_cross_sections,mie_cross_sections,
    dipole_cross_sections,CrossSections,RadialHeatModel,ThermalProperties,
    thermal_scales,ensemble_energy,ROOT)
from .experiments import csv,save_json


def run():
    p=ThermalProperties(); data=GoldData()
    grid=np.load(ROOT/'results/optical_grid.npz')
    waves=grid['wavelength_m']; radii=grid['radius_m']
    ab=grid['sigma_abs_m2']; sc=grid['sigma_sca_m2']; ex=grid['sigma_ext_m2']
    metrics=dict(optical_points=int(ab.size),min_absorption_m2=float(ab.min()),
        max_extinction_balance_relative=float(np.max(np.abs(ex-ab-sc)/ex)))
    # All measured knots recovered; all interpolated n,k positive. Alternative
    # interpolation quantifies model uncertainty, not solver convergence.
    nk=data.index(data.wavelength_m)
    metrics['max_knot_absolute_error']=float(np.max(np.abs(nk-(data.table[:,1]+1j*data.table[:,2]))))
    dense=np.linspace(data.wavelength_m[0],data.wavelength_m[-1],5000)
    metrics['minimum_interpolated_n']=float(data.index(dense).real.min())
    metrics['minimum_interpolated_kappa']=float(data.index(dense).imag.min())
    checks=[]; interp=[]; mie_conv=[]; diplimit=[]
    for n in grid['n_medium']:
        for a in radii:
            c=optical_cross_sections(a,waves,n,data=data)
            cplus=mie_cross_sections(a,waves,n,data.index(waves),extra_orders=8)
            mie_conv.append(dict(radius_m=a,n_medium=n,max_abs_relative=float(np.max(np.abs(cplus.absorption_m2/c.absorption_m2-1)))))
            linear=optical_cross_sections(a,waves,n,interpolation='linear',data=data)
            peak=int(np.argmax(c.absorption_m2)); jlin=int(np.argmax(linear.absorption_m2))
            interp.append(dict(radius_m=a,n_medium=n,max_spectral_relative=float(np.max(np.abs(linear.absorption_m2/c.absorption_m2-1))),peak_cross_section_relative=linear.absorption_m2[jlin]/c.absorption_m2[peak]-1,peak_shift_m=waves[jlin]-waves[peak]))
            for w in np.unique(np.r_[waves[peak],400e-9,532e-9,800e-9,1000e-9]):
                c=optical_cross_sections(a,w,n,data=data)
                qe,qs,*_=miepython.efficiencies(data.index(w).conjugate(),2*a,w,n_env=n)
                independent=(qe-qs)*np.pi*a*a
                checks.append(dict(radius_m=a,wavelength_m=w,n_medium=n,abs_relative=float(c.absorption_m2/independent-1),ext_relative=float(c.extinction_m2/(qe*np.pi*a*a)-1),sca_relative=float(c.scattering_m2/(qs*np.pi*a*a)-1)))
    for a in [5e-9,2e-9,1e-9,.5e-9,.1e-9]:
        c=optical_cross_sections(a,532e-9,data=data)
        d=optical_cross_sections(a,532e-9,model='dipole',data=data)
        diplimit.append(dict(radius_m=a,size_parameter=float(c.size_parameter),abs_relative=float(d.absorption_m2/c.absorption_m2-1)))
    csv('independent_mie.csv',checks); csv('mie_convergence.csv',mie_conv)
    csv('interpolation_sensitivity.csv',interp); csv('dipole_limit.csv',diplimit)
    metrics.update(max_miepython_abs_relative=float(max(abs(r['abs_relative']) for r in checks)),max_miepython_ext_relative=float(max(abs(r['ext_relative']) for r in checks)),max_mie_series_relative=float(max(r['max_abs_relative'] for r in mie_conv)),max_interpolation_spectral_relative=float(max(r['max_spectral_relative'] for r in interp)),max_interpolation_peak_relative=float(max(abs(r['peak_cross_section_relative']) for r in interp)),max_interpolation_peak_shift_m=float(max(abs(r['peak_shift_m']) for r in interp)),dipole_limit_at_0p1nm_relative=abs(diplimit[-1]['abs_relative']))
    # Optical maximum refinement around each host's discrete maximum, keeping
    # radii in the original interval. This is a local numerical resolution check.
    refined=[]
    for ni,n in enumerate(grid['n_medium']):
        ai,wi=np.unravel_index(np.argmax(ab[ni]),ab[ni].shape)
        test_r=np.arange(max(5,radii[ai]/1e-9-1),min(100,radii[ai]/1e-9+1)+.01,.1)*1e-9
        best=None
        for a in test_r:
            fit=minimize_scalar(lambda nm:-float(optical_cross_sections(a,nm*1e-9,n,data=data).absorption_m2),bounds=(max(400,waves[wi]/1e-9-3),min(1000,waves[wi]/1e-9+3)),method='bounded',options={'xatol':1e-9})
            row=dict(n_medium=n,radius_m=a,wavelength_m=fit.x*1e-9,absorption_m2=-fit.fun,relative_gain=-fit.fun/ab[ni,ai,wi]-1)
            if best is None or row['absorption_m2']>best['absorption_m2']: best=row
        refined.append(best)
    csv('optical_refinement.csv',refined)
    thermal=[]
    for a in [20e-9,80e-9,100e-9]:
        for dt in [100e-12,1e-9,10e-9,100e-9]:
            for g in [50e6,100e6,200e6]:
                hp=replace(p,interface_W_m2K=g)
                ref=RadialHeatModel(a,hp,cells=720)
                tref=np.array([dt,dt+10e-9,dt+100e-9])
                rr=ref.pulse(1,1e-15,dt,tref)
                for cells in [90,180,360]:
                    h=RadialHeatModel(a,hp,cells=cells); r=h.pulse(1,1e-15,dt,tref)
                    peak,pt=h.peak_water(1e-15/dt,dt,dt+300e-9)
                    peakref,ptref=ref.peak_water(1e-15/dt,dt,dt+300e-9)
                    thermal.append(dict(radius_m=a,duration_s=dt,interface_W_m2K=g,cells=cells,
                        gold_peak_relative=r['gold_rise_K'][0]/rr['gold_rise_K'][0]-1,
                        liquid_100ns_relative=r['liquid_net_energy_J'][-1]/rr['liquid_net_energy_J'][-1]-1,
                        water_surface_peak_relative=peak/peakref-1,water_peak_time_s=pt,reference_water_peak_time_s=ptref,
                        max_energy_residual_relative=float(np.max(np.abs(r['balance_residual_J']))/1e-15)))
    csv('thermal_convergence.csv',thermal)
    prod=[r for r in thermal if r['cells']==180]
    conservative=[r for r in prod if r['duration_s']>=1e-9]
    metrics.update(max_heat_gold_relative=float(max(abs(r['gold_peak_relative']) for r in conservative)),max_heat_water_surface_relative=float(max(abs(r['water_surface_peak_relative']) for r in conservative)),max_heat_energy_residual_relative=float(max(abs(r['max_energy_residual_relative']) for r in thermal)),max_heat_liquid_transfer_relative=float(max(abs(r['liquid_100ns_relative']) for r in conservative)))
    # Domain, steady state, analytic spherical resistances, and ODE propagation.
    other=[]
    for a in [20e-9,80e-9,100e-9]:
        for outer in [3e-6,6e-6]:
            h=RadialHeatModel(a,cells=360,outer_radius_m=outer)
            v=h.temperatures([1.],1e-7,2.)[:,0]
            rw=1/(4*np.pi*p.k_water_W_mK)*(1/a-1/outer)
            ri=1/(4*np.pi*a*a*p.interface_W_m2K)
            other.append(dict(radius_m=a,outer_radius_m=outer,kind='steady',relative_error=float(v[0]/(1e-7*(rw+ri))-1)))
    h=RadialHeatModel(20e-9,cells=90); dt=10e-9
    rate=-h.conductance/h.capacity[:,None]
    source=np.zeros(len(h.capacity)); source[0]=1e-7/h.capacity[0]
    ode=solve_ivp(lambda t,y:rate@y+source,[0,dt],np.zeros(len(source)),method='BDF',rtol=1e-10,atol=1e-12,jac=rate)
    exact=h.temperatures([dt],1e-7,dt)[:,0]
    metrics['independent_BDF_gold_relative']=float(abs(ode.y[0,-1]/exact[0]-1))
    # Independently integrate interface flux; thermodynamic E_input-E_Au is not
    # by itself an independent numerical energy test.
    conductance=-h.conductance[0,1]
    def flux(t):
        v=h.temperatures([t],1e-7,dt)[:,0]
        return conductance*(v[0]-v[1])
    heat=quad(flux,0,dt,epsabs=1e-25,epsrel=1e-10)[0]+quad(flux,dt,110e-9,epsabs=1e-25,epsrel=1e-10)[0]
    r=h.pulse(1,1e-15,dt,[110e-9])
    metrics['interface_quadrature_relative']=float(abs(heat/r['liquid_net_energy_J'][0]-1))
    r3=RadialHeatModel(80e-9,cells=360,outer_radius_m=3e-6).pulse(1,1e-15,100e-9,[100e-9,200e-9])
    r6=RadialHeatModel(80e-9,cells=360,outer_radius_m=6e-6).pulse(1,1e-15,100e-9,[100e-9,200e-9])
    metrics['outer_domain_gold_relative']=float(np.max(np.abs(r3['gold_rise_K']/r6['gold_rise_K']-1)))
    csv('analytic_thermal_controls.csv',other)
    # Finite-depth quadrature, thin/opaque limits and exact energy partition.
    energy=[]
    for n in [0,1e10,1e17,1e22]:
        a=20e-9; c=optical_cross_sections(a,532e-9)
        e=ensemble_energy(n,a,1e-3,np.pi*(100e-6)**2,1,c)
        # An opaque layer concentrates the integrand close to z=0. Resolve
        # absorption lengths explicitly so quadrature cannot skip that layer.
        mu=float(n*c.extinction_m2)
        points=[] if mu==0 else [v/mu for v in [1.,10.,100.] if 0<v/mu<1e-3]
        integrated=quad(lambda z:n*c.absorption_m2*np.exp(-mu*z),0,1e-3,points=points,epsabs=1e-13)[0]
        energy.append(dict(number_density_m3=n,integrated_absorption_fraction=integrated,formula_absorption_fraction=float(e['absorption_fraction']),partition_error=float((e['absorbed_J']+e['scattered_J']+e['transmitted_J'])/e['incident_J']-1)))
    csv('beam_energy_checks.csv',energy)
    metrics['max_finite_beam_partition_error']=max(abs(r['partition_error']) for r in energy)
    metrics['max_depth_quadrature_absolute']=max(abs(r['integrated_absorption_fraction']-r['formula_absorption_fraction']) for r in energy)
    # Literature comparisons are theoretical references, not experimental data.
    comp=[]
    for a,w,kind,reported in [(20e-9,528e-9,'Jain2006_calculated_absorption_m2',2.93e-15),(40e-9,560e-9,'Jain2006_calculated_scattering_m2',1.23e-14),(10e-9,530e-9,'Baffou2013_estimated_CW_K',5.)]:
        c=optical_cross_sections(a,w)
        computed=float(c.absorption_m2 if 'absorption' in kind else c.scattering_m2 if 'scattering' in kind else 1e9*c.absorption_m2/(4*np.pi*.6*a))
        comp.append(dict(reference=kind,radius_m=a,wavelength_m=w,reported_value=reported,computed_value=computed,relative_difference=computed/reported-1))
    csv('literature_comparison.csv',comp)
    save_json('validation_summary.json',metrics)
    print(json.dumps(metrics,indent=2))


if __name__=='__main__': run()
