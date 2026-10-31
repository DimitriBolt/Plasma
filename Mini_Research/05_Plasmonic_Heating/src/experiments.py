"""Deterministic MR5 design. Run from this study: python -m src.experiments."""
from dataclasses import asdict, replace
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import numpy as np
import pandas as pd
from .plasmonics import (ROOT, GoldData, optical_cross_sections, CrossSections,
                        RadialHeatModel, ThermalProperties, thermal_scales,
                        steady_temperature, ensemble_energy, ideal_energy_transfer_upper_bound)

RESULTS = ROOT/'results'
BASE_SHA = '803079ef68c0d3b02514dc88e211892079beeaff'


def native(value):
    if isinstance(value, dict): return {k:native(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [native(v) for v in value]
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    return value


def save_json(name, data):
    (RESULTS/name).write_text(json.dumps(native(data),indent=2,allow_nan=False)+'\n')


def csv(name, rows):
    pd.DataFrame(rows).to_csv(RESULTS/name,index=False,float_format='%.16g')


def input_inventory():
    repo=ROOT.parents[1]
    if not (repo/'.git').exists():
        # A git archive has no object database. Its saved inventory is the
        # immutable reference; a normal clone additionally checks base objects.
        inventory=json.loads((RESULTS/'input_inventory.json').read_text())
        for row in inventory['files']:
            if hashlib.sha256((repo/row['path']).read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError(f"Archived upstream differs: {row['path']}")
        return
    names=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE_SHA,'Mini_Research/01_Bubble_Dynamics','Mini_Research/02_Compressible_Bubble_Dynamics','Mini_Research/03_In_Bubble_Thermodynamics','Mini_Research/04_Temperature_Barrier','Plasma-plan.tex','Library'],cwd=repo,text=True).splitlines()
    rows=[]
    for name in names:
        current=(repo/name).read_bytes()
        published=subprocess.check_output(['git','show',f'{BASE_SHA}:{name}'],cwd=repo)
        if current!=published: raise ValueError(f'Upstream differs: {name}')
        rows.append(dict(path=name,sha256=hashlib.sha256(current).hexdigest()))
    save_json('input_inventory.json',dict(base_sha=BASE_SHA,unchanged_inputs=len(rows),files=rows))


def run():
    RESULTS.mkdir(exist_ok=True)
    data=GoldData(); p=ThermalProperties()
    radii=np.arange(5,101,dtype=float)*1e-9
    waves=np.arange(400,1001,dtype=float)*1e-9
    hosts=[1.,1.33,1.5]
    grids=[]; resonances=[]
    absorption=np.empty((len(hosts),len(radii),len(waves)))
    scattering=np.empty_like(absorption); extinction=np.empty_like(absorption)
    dipole=np.empty_like(absorption)
    for ni,n in enumerate(hosts):
        for ai,a in enumerate(radii):
            c=optical_cross_sections(a,waves,n,data=data)
            d=optical_cross_sections(a,waves,n,model='dipole',data=data)
            for dest,values in [(absorption,c.absorption_m2),(scattering,c.scattering_m2),(extinction,c.extinction_m2),(dipole,d.absorption_m2)]: dest[ni,ai]=values
            j=int(np.argmax(c.absorption_m2)); je=int(np.argmax(c.extinction_m2))
            resonances.append(dict(radius_m=a,n_medium=n,lambda_abs_m=waves[j],lambda_ext_m=waves[je],sigma_abs_m2=c.absorption_m2[j],sigma_sca_m2=c.scattering_m2[j],sigma_ext_m2=c.extinction_m2[j],albedo=c.scattering_m2[j]/c.extinction_m2[j],size_parameter=c.size_parameter[j],sigma_abs_per_gold_volume_m1=c.absorption_m2[j]/(4*np.pi*a**3/3),dipole_relative_error_at_Mie_peak=d.absorption_m2[j]/c.absorption_m2[j]-1))
    np.savez_compressed(RESULTS/'optical_grid.npz',radius_m=radii,wavelength_m=waves,n_medium=hosts,sigma_abs_m2=absorption,sigma_sca_m2=scattering,sigma_ext_m2=extinction,sigma_abs_dipole_m2=dipole)
    csv('resonances.csv',resonances)
    # CSV cross-section function at explicitly relevant laboratory wavelengths.
    crossrows=[]
    for n in hosts:
        for a in [5e-9,10e-9,20e-9,40e-9,60e-9,80e-9,100e-9]:
            for w in [405e-9,488e-9,520e-9,532e-9,560e-9,633e-9,800e-9,1064e-9]:
                c=optical_cross_sections(a,w,n,data=data)
                crossrows.append(dict(radius_m=a,wavelength_m=w,n_medium=n,sigma_abs_m2=c.absorption_m2,sigma_sca_m2=c.scattering_m2,sigma_ext_m2=c.extinction_m2))
    csv('cross_sections.csv',crossrows)
    thermal=[]
    # Worst-interface conservative pulse cap at Au rise=30 K. Gold peak is
    # exactly pulse end for a passive diffusion network with positive step input.
    for a in radii:
        for g in [50e6,100e6,200e6]:
            hp=replace(p,interface_W_m2K=g)
            h=RadialHeatModel(a,hp)
            for dt in [100e-12,1e-9,10e-9,100e-9]:
                ts=np.array([dt,dt+1e-9,dt+10e-9,dt+100e-9])
                r=h.pulse(1,1e-15,dt,ts)
                thermal.append(dict(radius_m=a,interface_W_m2K=g,duration_s=dt,
                    gold_peak_K_per_J=(r['gold_rise_K'][0]/1e-15),
                    liquid_fraction_at_end=r['liquid_net_energy_J'][0]/1e-15,
                    liquid_fraction_10ns_after=r['liquid_net_energy_J'][2]/1e-15,
                    liquid_fraction_100ns_after=r['liquid_net_energy_J'][3]/1e-15,
                    **thermal_scales(a,hp)))
    csv('thermal_response.csv',thermal)
    tf=pd.DataFrame(thermal)
    water_res=pd.DataFrame(resonances).query('n_medium == 1.33')
    allowed=[]; formal=[]
    for ai,a in enumerate(radii):
        cs=CrossSections(absorption[1,ai],scattering[1,ai],extinction[1,ai],2*np.pi*1.33*a/waves)
        for dt in [1e-9,10e-9,100e-9]:
            tr=tf[(np.isclose(tf.radius_m/a,1,rtol=1e-12))&(tf.interface_W_m2K==50e6)&(tf.duration_s==dt)].iloc[0]
            cap=30/tr.gold_peak_K_per_J
            # Finite laser design: F<=1000 J/m2, Ipeak<=1e10 W/m2.
            f=np.minimum(np.minimum(cap/cs.absorption_m2,1000.),1e10*dt)
            energies=f*cs.absorption_m2
            tied=np.isclose(energies,np.max(energies),rtol=1e-12,atol=0)
            j=int(np.argmax(np.where(tied,cs.absorption_m2,-np.inf)))
            # a=20..100 nm is the conservative geometric subdomain, not a new
            # optical validation threshold. The full 5..100 nm results persist.
            allowed.append(dict(radius_m=a,wavelength_m=waves[j],duration_s=dt,
                fluence_J_m2=f[j],peak_intensity_W_m2=f[j]/dt,
                sigma_abs_m2=cs.absorption_m2[j],sigma_sca_m2=cs.scattering_m2[j],
                energy_absorbed_J=energies[j],gold_peak_rise_K=energies[j]*tr.gold_peak_K_per_J,
                liquid_fraction_100ns_after=tr.liquid_fraction_100ns_after,
                energy_liquid_100ns_after_J=energies[j]*tr.liquid_fraction_100ns_after,
                conservative_geometry=bool(20e-9-1e-20<=a<=100e-9+1e-20)))
    csv('bounded_single_particle.csv',allowed)
    safe=[r for r in allowed if r['conservative_geometry']]
    single=max(safe,key=lambda r:r['energy_absorbed_J'])
    # Optimize a different objective: fraction of a finite beam absorbed in a
    # layer. All water spectra searched; intensity at entrance sets worst Au rise.
    # Independent scattering/dispersion gates are explicit, not a proof against
    # rare aggregates. Thermal overlap is checked only over pulse+100 ns.
    ngrid=np.geomspace(1e14,1e19,61)
    erows=[]
    for a in radii[(radii>=20e-9-1e-20)&(radii<=100e-9+1e-20)]:
        ai=int(round(a/1e-9))-5
        tr=tf[np.isclose(tf.radius_m/a,1,rtol=1e-12)&(tf.interface_W_m2K==50e6)&(tf.duration_s==100e-9)].iloc[0]
        cap=30/tr.gold_peak_K_per_J
        ab=absorption[1,ai]; sc=scattering[1,ai]; ex=extinction[1,ai]
        f=np.minimum(np.minimum(cap/ab,1000.),1000.)
        for length in [1e-4,1e-3,1e-2]:
            for b in [25e-6,100e-6,500e-6]:
                density=ngrid[:,None]
                tau=density*ex*length
                etaa=ab/ex*(-np.expm1(-tau))
                spacing=density**(-1/3)
                valid=(density*4*np.pi*a**3/3<=1e-3)&(spacing>=10*a)&(density*sc*length<=0.1)&(length<=np.pi*1.33*b*b/waves)&(spacing>=2*a+4*np.sqrt(p.k_water_W_mK/(p.rho_water_kg_m3*p.cp_water_J_kgK)*200e-9))
                if not np.any(valid): continue
                for objective,score in [('fraction',etaa),('energy',etaa*f[None,:])]:
                    candidates=np.where(valid,score,-np.inf)
                    ix,jx=np.unravel_index(np.argmax(candidates),candidates.shape)
                    eta=float(etaa[ix,jx]); n=float(ngrid[ix]); area=np.pi*b*b
                    erows.append(dict(optimized_objective=objective,radius_m=a,wavelength_m=waves[jx],number_density_m3=n,thickness_m=length,beam_radius_m=b,duration_s=100e-9,fluence_J_m2=f[jx],incident_J=f[jx]*area,absorbed_J=eta*f[jx]*area,absorption_fraction=eta,tau_ext=tau[ix,jx],tau_sca=n*sc[jx]*length,volume_fraction=n*4*np.pi*a**3/3,spacing_m=n**(-1/3),liquid_energy_100ns_after_J=eta*f[jx]*area*tr.liquid_fraction_100ns_after,layer_mean_sensible_rise_K=eta*f[jx]/(p.rho_water_kg_m3*p.cp_water_J_kgK*length),sigma_abs_m2=ab[jx],sigma_sca_m2=sc[jx]))
    csv('bounded_ensemble.csv',erows)
    ensemble=max(erows,key=lambda r:r['absorbed_J'])
    fraction=max(erows,key=lambda r:r['absorption_fraction'])
    # Formal intensity/fluence scan: no censoring at boiling/melting or at 30 K.
    for a in [20e-9,40e-9,80e-9,100e-9]:
        wr=water_res.iloc[int(round(a/1e-9))-5]
        c=optical_cross_sections(a,wr.lambda_abs_m,data=data)
        for dt in [100e-15,1e-12,100e-12,1e-9,10e-9,100e-9]:
            h=RadialHeatModel(a)
            water_per_J,water_peak_time=h.peak_water(1e-15/dt,dt,dt+300e-9)
            water_per_J/=1e-15
            for fluence in np.geomspace(1e-3,1e3,37):
                e=float(fluence*c.absorption_m2)
                instantaneous=e/thermal_scales(a)['capacity_J_K']
                rise=float(h.temperatures([dt],e/dt,dt)[0,0])
                valid_t=dt>=10*max(p.electron_phonon_s,thermal_scales(a)['gold_diffusion_s'])
                formal.append(dict(radius_m=a,wavelength_m=wr.lambda_abs_m,duration_s=dt,fluence_J_m2=fluence,peak_intensity_W_m2=fluence/dt,absorbed_J=e,formal_lumped_gold_peak_rise_K=rise,formal_water_surface_peak_rise_K=water_per_J*e,formal_water_peak_time_s=water_peak_time,adiabatic_gold_sensible_rise_K=instantaneous,one_temperature_time_gate=valid_t,moderate_rise_gate=rise<=30,bulk_melt_warning=p.ambient_K+rise>=1337.33,ambient_boiling_diagnostic=p.ambient_K+water_per_J*e>=373.15))
    csv('formal_regimes.csv',formal)
    concentration=[]
    for a in [20e-9,40e-9,80e-9]:
        c=optical_cross_sections(a,532e-9,data=data)
        for n in np.geomspace(1e13,1e20,141):
            e=ensemble_energy(n,a,1e-3,np.pi*(100e-6)**2,1,c)
            concentration.append(dict(radius_m=a,number_density_m3=n,**e))
    csv('concentration.csv',concentration)
    cw=[]
    for a in [10e-9,20e-9,40e-9,80e-9]:
        c=optical_cross_sections(a,532e-9,data=data)
        for intensity in [1e4,1e6,1e8,1e9,1e10]:
            power=float(intensity*c.absorption_m2); v=steady_temperature(a,power)
            cw.append(dict(radius_m=a,wavelength_m=532e-9,intensity_W_m2=intensity,power_absorbed_W=power,water_surface_rise_K=float(v['water_rise_K']),gold_rise_K=float(v['gold_rise_K'])))
    csv('cw_heating.csv',cw)
    # Saved time series and profile, 20 nm particle, 10 ns, F=1 J/m2.
    reference=dict(radius_m=20e-9,wavelength_m=532e-9,fluence_J_m2=1.,duration_s=10e-9)
    c=optical_cross_sections(reference['radius_m'],reference['wavelength_m'])
    h=RadialHeatModel(reference['radius_m'])
    time=np.unique(np.r_[0,np.geomspace(1e-12,1e-6,700),reference['duration_s']])
    rr=h.pulse(reference['fluence_J_m2'],c.absorption_m2,reference['duration_s'],time)
    np.savez_compressed(RESULTS/'pulse_reference.npz',radius_water_m=h.r,**rr)
    peak,timepeak=h.peak_water(float(c.absorption_m2)/reference['duration_s'],reference['duration_s'],1e-6)
    reference.update(sigma_abs_m2=float(c.absorption_m2),gold_peak_rise_K=float(h.temperatures([reference['duration_s']],float(c.absorption_m2)/reference['duration_s'],reference['duration_s'])[0,0]),water_surface_peak_rise_K=peak,water_peak_time_s=timepeak)
    # Full-grid extrema have explicitly separate optical objectives.
    maxima={}
    for name,v in [('single_cross_section',absorption),('per_gold_volume',absorption/(4*np.pi*radii[None,:,None]**3/3))]:
        ni,ai,wi=np.unravel_index(np.argmax(v),v.shape)
        maxima[name]=dict(radius_m=radii[ai],wavelength_m=waves[wi],n_medium=hosts[ni],objective_value=v[ni,ai,wi],sigma_abs_m2=absorption[ni,ai,wi])
    upper=[]
    for label,e in [('bounded_single',single['energy_absorbed_J']),('bounded_ensemble',ensemble['absorbed_J'])]:
        for rad in [1e-6,6e-6,10e-6]:
            upper.append(dict(source=label,selected_radius_m=rad,**ideal_energy_transfer_upper_bound(e,4*np.pi*rad**3/3)))
    csv('ideal_upper_bound.csv',upper)
    save_json('summary.json',dict(optical_grid_maxima=maxima,bounded_single_particle=single,bounded_ensemble_energy=ensemble,bounded_ensemble_fraction=fraction,pulse_reference=reference))
    save_json('parameters.json',dict(thermal_properties=asdict(p),optical_radius_m=radii,optical_wavelength_m=waves,host_indices=hosts,thermal_cells=180,outer_radius_m=3e-6,conservative_radius_interval_m=[20e-9,100e-9],pulse_durations_s=[1e-9,10e-9,100e-9],max_fluence_J_m2=1000,max_intensity_W_m2=1e10,max_gold_rise_K=30,conservative_interface_W_m2K=50e6,number_density_m3=ngrid,layer_thickness_m=[1e-4,1e-3,1e-2],beam_radius_m=[25e-6,100e-6,500e-6],max_volume_fraction=1e-3,max_tau_sca=0.1,min_spacing_over_a=10,thermal_nonoverlap_time_s=200e-9))
    save_json('environment.json',dict(python=platform.python_version(),platform=platform.platform(),packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','pandas','miepython']}))
    input_inventory()
    print(json.dumps(native(dict(single=single,ensemble=ensemble,fraction=fraction,reference=reference)),indent=2))


if __name__=='__main__': run()
