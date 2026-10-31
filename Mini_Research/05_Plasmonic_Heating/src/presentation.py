"""Saved-data scientific figures (JPG only), English labels, generated tables."""
import json
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch,FancyBboxPatch
from .plasmonics import ROOT,RadialHeatModel,optical_cross_sections,steady_temperature


def run():
    out=ROOT/'figures'; out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,
        'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
        'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
    def save(fig,name):
        fig.savefig(out/(name+'.jpg'),dpi=300,bbox_inches='tight',pil_kwargs={'quality':95});plt.close(fig)
    g=np.load(ROOT/'results/optical_grid.npz'); a=g['radius_m']; w=g['wavelength_m']
    r=pd.read_csv(ROOT/'results/resonances.csv'); r=r[r.n_medium==1.33]
    s=json.loads((ROOT/'results/summary.json').read_text())
    fig,ax=plt.subplots(figsize=(7.1,1.85));ax.set(xlim=(-.2,10.2),ylim=(0,2));ax.axis('off')
    for x,label,color in [(0.1,'Laser\nlocal I(t), fluence F','#e8edf2'),(2.7,'Au sphere\nMie absorption','#fff1c4'),(5.3,'Finite interface\nconductance G','#e8edf2'),(7.9,'Water\nheat diffusion','#e2f1f6')]:
        ax.add_patch(FancyBboxPatch((x,.45),2,1.1,boxstyle='round,pad=.12',fc=color,ec='.25'))
        ax.text(x+1,1,label,ha='center',va='center')
    for x in [2.2,4.8,7.4]:ax.add_patch(FancyArrowPatch((x,1),(x+.35,1),arrowstyle='-|>',mutation_scale=12,color='.3'))
    ax.text(5,.03,'Incident energy = transmitted + first-event scattered + absorbed; absorbed = Au storage + liquid transfer',ha='center',fontsize=7.7)
    save(fig,'01_energy_path')
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.7),layout='constrained')
    for rad in [10,20,40,80,100]:
        i=rad-5; axes[0].plot(w*1e9,g['sigma_abs_m2'][1,i]/1e-14,label=f'a = {rad} nm')
        axes[1].plot(w*1e9,g['sigma_sca_m2'][1,i]/g['sigma_ext_m2'][1,i],label=f'{rad} nm')
    axes[0].set(xlabel='Vacuum wavelength (nm)',ylabel=r'$\sigma_{abs}$ ($10^{-14}$ m$^2$)',title='Absorption in water (n = 1.33)')
    axes[1].set(xlabel='Vacuum wavelength (nm)',ylabel=r'$\sigma_{sca}/\sigma_{ext}$',ylim=(0,1),title='Fraction of extinction due to scattering')
    axes[0].legend(ncol=2);axes[1].legend(ncol=2)
    save(fig,'02_optical_spectra')
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.7),layout='constrained')
    for rad in [5,10,20,40]:
        i=rad-5; axes[0].plot(w*1e9,g['sigma_abs_dipole_m2'][1,i]/g['sigma_abs_m2'][1,i]-1,label=f'a = {rad} nm')
    axes[0].set(xlabel='Vacuum wavelength (nm)',ylabel='Dipole / Mie absorption - 1',title='Radiative correction preserves passivity')
    axes[0].legend(ncol=2)
    axes[1].plot(r.radius_m*1e9,r.lambda_abs_m*1e9,label='Absorption maximum')
    axes[1].plot(r.radius_m*1e9,r.lambda_ext_m*1e9,'--',label='Extinction maximum')
    axes[1].set(xlabel='Sphere radius (nm)',ylabel='Vacuum wavelength (nm)',title='Different spectral objectives')
    axes[1].legend()
    save(fig,'03_dipole_mie')
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.7),layout='constrained')
    im=axes[0].pcolormesh(w*1e9,a*1e9,g['sigma_abs_m2'][1]/1e-14,shading='auto',cmap='cividis')
    fig.colorbar(im,ax=axes[0],label=r'$\sigma_{abs}$ ($10^{-14}$ m$^2$)')
    axes[0].plot(r.lambda_abs_m*1e9,r.radius_m*1e9,color='white',lw=1)
    im=axes[1].pcolormesh(w*1e9,a*1e9,g['sigma_sca_m2'][1]/g['sigma_ext_m2'][1],vmin=0,vmax=1,shading='auto',cmap='cividis')
    fig.colorbar(im,ax=axes[1],label='Scattering / extinction')
    for ax in axes:ax.set(xlabel='Vacuum wavelength (nm)',ylabel='Sphere radius (nm)');ax.grid(False)
    save(fig,'04_radius_wavelength_map')
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.6),layout='constrained')
    rad=20e-9;c=optical_cross_sections(rad,532e-9);radii=np.geomspace(rad,rad*30,150)
    for intensity in [1e7,1e8,1e9]:
        val=steady_temperature(rad,float(c.absorption_m2)*intensity,radii)
        line=axes[0].plot(radii/rad,val['water_rise_K'],label=f'I = {intensity:.0e} W/m²')[0]
        axes[0].plot([.6,1],[val['gold_rise_K']]*2,color=line.get_color(),ls='--')
    axes[0].set(xscale='log',xlabel='Distance from center / a',ylabel='Steady temperature rise (K)',title='a = 20 nm, 532 nm; Au dashed, water solid')
    axes[0].legend(fontsize=7)
    th=pd.read_csv(ROOT/'results/thermal_response.csv');th=th[(th.interface_W_m2K==1e8)&(th.duration_s==1e-8)]
    axes[1].loglog(th.radius_m*1e9,th.gold_diffusion_s*1e9,label='Au diffusion')
    axes[1].loglog(th.radius_m*1e9,th.water_diffusion_s*1e9,label='Water diffusion')
    axes[1].loglog(th.radius_m*1e9,th.interface_relaxation_s*1e9,label='Interface RC')
    axes[1].axhline(.0017,color='.2',ls=':',label='Electron-phonon')
    axes[1].set(xlabel='Sphere radius (nm)',ylabel='Characteristic time (ns)',title='Separate thermal time scales')
    axes[1].legend(fontsize=7)
    save(fig,'05_steady_profile_and_scales')
    rr=np.load(ROOT/'results/pulse_reference.npz'); t=rr['time_s']; good=t>0
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.6),layout='constrained')
    axes[0].semilogx(t[good]*1e9,rr['gold_rise_K'][good],label='Gold')
    axes[0].semilogx(t[good]*1e9,rr['water_rise_K'][0,good],label='First water cell')
    axes[0].axvline(10,color='.4',ls=':');axes[0].legend()
    axes[0].set(xlabel='Time from pulse start (ns)',ylabel='Temperature rise (K)',title='a = 20 nm, F = 1 J/m², 10 ns, 532 nm')
    ein=float(s['pulse_reference']['sigma_abs_m2'])
    for key,label in [('energy_input_J','Absorbed to time t'),('gold_energy_J','Gold storage'),('water_energy_J','Water storage'),('liquid_net_energy_J','Net transfer to liquid')]:axes[1].semilogx(t[good]*1e9,rr[key][good]/ein,label=label)
    axes[1].set(xlabel='Time from pulse start (ns)',ylabel='Energy / final absorbed energy',title='Finite pulse: energy is counted once',ylim=(-.02,1.05));axes[1].legend(fontsize=7)
    save(fig,'06_pulse_heating')
    cr=pd.read_csv(ROOT/'results/concentration.csv')
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.7),layout='constrained')
    for ax,rad in zip(axes,[20,80]):
        q=cr[np.isclose(cr.radius_m/(rad*1e-9),1,rtol=1e-12)]
        for key,label,style in [('absorption_fraction','Absorption','-'),('scattering_fraction','First-event scattering','--'),('transmission_fraction','Transmission',':')]:ax.semilogx(q.number_density_m3,q[key],style,label=label)
        ix=q[q.tau_sca>.1]
        if len(ix):ax.axvspan(ix.number_density_m3.iloc[0],q.number_density_m3.max(),color='.85',zorder=-2)
        ax.set(xlabel=r'Number density (m$^{-3}$)',ylabel='Fraction of incident energy',ylim=(0,1.03),title=f'a = {rad} nm, 532 nm, L = 1 mm')
        ax.legend(fontsize=7)
    save(fig,'07_ensemble_attenuation')
    f=pd.read_csv(ROOT/'results/formal_regimes.csv');f=f[np.isclose(f.radius_m/40e-9,1,rtol=1e-12)]
    durations=np.sort(f.duration_s.unique()); fluences=np.sort(f.fluence_J_m2.unique())
    pivot=f.pivot(index='duration_s',columns='fluence_J_m2',values='formal_lumped_gold_peak_rise_K')
    fig,ax=plt.subplots(figsize=(7.1,2.8),layout='constrained')
    im=ax.pcolormesh(fluences,durations*1e9,np.log10(pivot.values),shading='nearest',cmap='cividis')
    fig.colorbar(im,ax=ax,label=r'log$_{10}$(formal Au rise / K)')
    limits=f.groupby('duration_s').first()
    for rise,label,style in [(30,'Moderate-rise screen: 30 K','-'),(1337.33-293.15,'Bulk Au melt diagnostic','--')]:
        cap=rise/(limits.formal_lumped_gold_peak_rise_K/limits.fluence_J_m2)
        ax.plot(cap,durations*1e9,style,color='white',lw=1.5,label=label)
    cap=(373.15-293.15)/(limits.formal_water_surface_peak_rise_K/limits.fluence_J_m2)
    ax.plot(cap,durations*1e9,':',color='white',lw=1.6,label='Water surface reaches 373.15 K')
    ax.axhspan(durations.min()*1e9,.126,facecolor='none',hatch='///',edgecolor='white',lw=0)
    ax.set(xscale='log',yscale='log',xlim=(fluences.min(),fluences.max()),ylim=(durations.min()*1e9,durations.max()*1e9),xlabel='Fluence (J/m²)',ylabel='Pulse duration (ns)',title='a = 40 nm; formal continuation, not phase-change predictions')
    ax.legend(fontsize=7,loc='lower right',facecolor='white');ax.grid(False)
    save(fig,'08_formal_regime_map')
    # Modest tables for the report, from saved values only.
    sel=r[np.round(r.radius_m*1e9).isin([5,10,20,40,60,80,100])]
    lines=[r'\begin{tabular}{rrrrrr}',r'\toprule',r'$a$ (nm)&$\lambda_{abs}$ (nm)&$\lambda_{ext}$ (nm)&$\sigma_{abs}$ (m$^2$)&$\sigma_{sca}/\sigma_{ext}$&Dipole error\\',r'\midrule']
    for row in sel.itertuples():lines.append(f'{row.radius_m*1e9:.0f} & {row.lambda_abs_m*1e9:.0f} & {row.lambda_ext_m*1e9:.0f} & {row.sigma_abs_m2:.3e} & {row.albedo:.3f} & {row.dipole_relative_error_at_Mie_peak*100:.1f}\\% \\\\')
    lines.extend([r'\bottomrule',r'\end{tabular}'])
    (ROOT/'results/table_resonances.tex').write_text('\n'.join(lines)+'\n')
    val=json.loads((ROOT/'results/validation_summary.json').read_text())
    upper=s['bounded_single_particle']['energy_absorbed_J']/(1000*4187*4*np.pi*(6e-6)**3/3)
    macros={'SinglePJ':f"{s['bounded_single_particle']['energy_absorbed_J']*1e12:.3f}",
            'SingleLiquidPJ':f"{s['bounded_single_particle']['energy_liquid_100ns_after_J']*1e12:.3f}",
            'EnsembleMicroJ':f"{s['bounded_ensemble_energy']['absorbed_J']*1e6:.3f}",
            'SingleUpperK':f'{upper:.6f}',
            'HeatGoldPercent':f"{val['max_heat_gold_relative']*100:.4f}",
            'HeatWaterPercent':f"{val['max_heat_water_surface_relative']*100:.4f}"}
    for key,field in [('HeatLiquidRelative','max_heat_liquid_transfer_relative'),('HeatBalanceRelative','max_heat_energy_residual_relative')]:
        mantissa,exponent=f'{val[field]:.2e}'.split('e')
        macros[key]=f'${mantissa}\\times10^{{{int(exponent)}}}$'
    (ROOT/'results/report_values.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items())+'\n')
    parameters=json.loads((ROOT/'results/parameters.json').read_text())
    rows=[]
    for key,value in parameters['thermal_properties'].items():
        unit=('K' if key=='ambient_K' else 's' if key=='electron_phonon_s' else
              'kg/m^3' if key.startswith('rho') else 'J/(kg K)' if key.startswith('cp') else
              'W/(m K)' if key.startswith('k_') else 'W/(m^2 K)')
        rows.append(dict(parameter=key,value=value,unit=unit,
            provenance='Baffou-Rigneault 2011 Table I' if key not in ['ambient_K','electron_phonon_s','interface_W_m2K'] else 'design bath / literature relaxation and representative interface as documented'))
    for key,value in parameters.items():
        if key=='thermal_properties':continue
        unit=('m' if key.endswith('_m') else 's' if key.endswith('_s') else
              'J/m^2' if key.endswith('_J_m2') else 'W/m^2' if key.endswith('_W_m2') else
              'K' if key.endswith('_K') else 'W/(m^2 K)' if key.endswith('_W_m2K') else
              'm^-3' if key.endswith('_m3') else 'dimensionless')
        rows.append(dict(parameter=key,value=json.dumps(value),unit=unit,
            provenance='explicit finite design or numerical resolution; not a measured threshold'))
    with (ROOT/'results/parameters.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)


if __name__=='__main__':run()
