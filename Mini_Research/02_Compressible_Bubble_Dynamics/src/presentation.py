"""Figures and LaTeX tables read saved CSVs; no mechanical integration."""
import gzip
import json
import os
from functools import lru_cache
os.environ.setdefault('MPLCONFIGDIR','/tmp/plasma_mr2_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .model import ROOT
from .experiments import read_csv, RESULTS, CANONICAL

FIGURES=ROOT/'figures'
COLORS={'RP':'#b65425','KM':'#126a91'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.grid':True,'grid.alpha':.18,'savefig.dpi':300,
                     'pdf.fonttype':42,'font.family':'DejaVu Sans'})


def numeric(name):
    rows=read_csv(RESULTS/name)
    for r in rows:
        for k,v in r.items():
            if v in ('True','False'): r[k]=v=='True'
            else:
                try: r[k]=float(v)
                except ValueError: pass
    return rows


@lru_cache(maxsize=16)
def trajectory(model,bar):
    with gzip.open(RESULTS/'trajectories'/f'{model}_{bar:.3f}bar.csv.gz','rt') as f:
        return np.genfromtxt(f,delimiter=',',names=True)


def save(fig,name):
    fig.tight_layout()
    for ext in ('pdf','png'): fig.savefig(FIGURES/f'{name}.{ext}',bbox_inches='tight')
    plt.close(fig)


def figures():
    FIGURES.mkdir(exist_ok=True)
    rows=numeric('amplitude_summary.csv')
    def series(model,k):
        s=sorted((r for r in rows if r['model']==model),key=lambda r:r['Pa_bar'])
        return np.array([r['Pa_bar'] for r in s]),np.array([r[k] if r['completed'] else np.nan for r in s])
    for field,scale,ylabel,name in [('R_m',1e6,'Radius (um)','01_radius'),
                                    ('U_m_s',1,'Wall velocity (m/s)','02_velocity')]:
        fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True)
        for ax,bar in zip(axes.flat,CANONICAL):
            for m in ('RP','KM'):
                x=trajectory(m,bar)
                ax.plot(x['t_s']*1e6,x[field]*scale,label=m,color=COLORS[m],lw=.9,alpha=.85)
            if field=='U_m_s' and bar==1.3: ax.set_yscale('symlog',linthresh=100)
            ax.set_title(f'{bar:g} bar'); ax.set_ylabel(ylabel)
            ax.legend(fontsize=8,loc='best')
        for ax in axes[-1]: ax.set_xlabel('Time (us)')
        save(fig,name)
    first={}
    for m in ('RP','KM'):
        r=next(r for r in rows if r['model']==m and r['Pa_bar']==1.3)
        first[m]=r['first_strong_time_s']
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    for j,m in enumerate(('RP','KM')):
        x=trajectory(m,1.3); tc=first[m]
        for i,(key,scale,label) in enumerate([('R_m',1e6,'Radius (um)'),('U_m_s',1,'Velocity (m/s)')]):
            ax=axes[i,j]; mask=abs(x['t_s']-tc)<.4e-6
            ax.plot((x['t_s'][mask]-tc)*1e9,x[key][mask]*scale,color=COLORS[m],lw=1)
            ax.set(xlabel='Time from own first strong minimum (ns)',ylabel=label,
                   title=f'{m}: first cycle, center {tc*1e6:.6f} us')
            ax.axvline(0,color='k',lw=.5)
    save(fig,'03_strong_collapse')
    for k,ylabel,name,log in [
        ('Rmin_m','Minimum radius (um)','04_minimum_radius',True),
        ('max_wall_Mach','Maximum wall Mach number','05_wall_mach',True),
        ('Vmax_over_Vmin','Global maximum / minimum volume','09_volume_ratio',True)]:
        fig,ax=plt.subplots(figsize=(8,4))
        for m in ('RP','KM'):
            a,b=series(m,k); mask=a>=1.
            ax.plot(a[mask],b[mask]*(1e6 if k=='Rmin_m' else 1),'.-',label=m,color=COLORS[m],ms=4,lw=.8)
        if log: ax.set_yscale('log')
        ax.set(xlabel='Acoustic amplitude (bar)',ylabel=ylabel)
        ax.legend(); ax.set_title('Five cycles; line gaps mark incomplete RP runs')
        save(fig,name)
    cycles=numeric('cycle_comparison.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,3.6))
    for m in ('RP','KM'):
        s=sorted((r for r in cycles if r['cycle']==1 and r['Pa_bar']>=1),key=lambda r:r['Pa_bar'])
        for ax,k,label in zip(axes,('time_s','duration_s'),('First-cycle deepest minimum (us)','Paired collapse duration (us)')):
            ax.plot([r['Pa_bar'] for r in s],[r[m+'_'+k]*1e6 for r in s],'.-',color=COLORS[m],label=m,lw=.8)
            ax.set(xlabel='Acoustic amplitude (bar)',ylabel=label); ax.legend()
    save(fig,'06_collapse_timing')
    fig,axes=plt.subplots(1,3,figsize=(11,3.5))
    for ax,k,label in zip(axes,('max_wall_Mach','max_liquid_pressure_change_over_bulk_modulus','min_D_over_R'),
                         ('Wall Mach number','Pressure-change / rho c^2 (proxy)','Minimum D/R')):
        a,b=series('KM',k); mask=a>=1; ax.plot(a[mask],b[mask],'.-',color=COLORS['KM'],ms=3)
        ax.set(xlabel='Amplitude (bar)',ylabel=label)
    fig.suptitle('KM diagnostics: continuous indicators, no universal validity cutoff',fontsize=11)
    save(fig,'07_applicability')
    fig,axes=plt.subplots(1,3,figsize=(11,3.8))
    for ax,bar in zip(axes,(1.,1.15,1.3)):
        for m in ('RP','KM'):
            x=trajectory(m,bar); ax.plot(x['R_m']*1e6,x['U_m_s'],lw=.65,color=COLORS[m],label=m)
        ax.set(xlabel='Radius (um)',ylabel='Wall velocity (m/s)',title=f'{bar:g} bar')
        if bar>1.: ax.set_yscale('symlog',linthresh=100)
        ax.legend()
    save(fig,'08_phase_portraits')
    comp=numeric('comparison.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    for ax,k,label in zip(axes,('Rmin_m','max_abs_U_m_s'),('Minimum radius','Maximum absolute speed')):
        s=[r for r in comp if r['Pa_bar']>=1]
        ax.plot([r['Pa_bar'] for r in s],[100*r['relative_difference_'+k] for r in s],'.-',color=COLORS['KM'],lw=.8)
        ax.set(xlabel='Amplitude (bar)',ylabel=f'{label}: (KM-RP)/RP (%)')
    save(fig,'10_relative_differences')
    fig,axes=plt.subplots(2,3,figsize=(11,6))
    for i,m in enumerate(('RP','KM')):
        x=trajectory(m,1.3)
        row=next(r for r in rows if r['model']==m and r['Pa_bar']==1.3)
        tc=row['collapse_time_s']
        for ax,width in zip(axes[i],(500.,1.,.0005) if m=='RP' else (500.,20.,1.)):
            delta=(x['t_s']-tc)*1e9; mask=abs(delta)<=width
            unit='ps' if width<.01 else 'ns'; scale=1e3 if unit=='ps' else 1
            ax.plot(delta[mask]*scale,x['R_m'][mask]*1e6,color=COLORS[m],lw=1,marker='.',ms=1.5)
            ax.set(xlabel=f'Time from own deepest minimum ({unit})',ylabel='Radius (um)',title=f'{m}, cycle {int(row["characteristic_minimum_cycle"])}: +/- {width*scale:g} {unit}')
    save(fig,'11_multiscale_collapse')
    km=sorted((r for r in rows if r['model']=='KM' and r['Pa_bar']>=1),key=lambda r:r['Pa_bar'])
    grid=np.round(np.arange(1.,1.30001,.005),3)
    z=np.full((2,len(grid)),np.nan)
    for r in km:
        j=int(np.argmin(abs(grid-r['Pa_bar'])))
        for i,k in enumerate(('max_wall_Mach','max_liquid_pressure_change_over_bulk_modulus')):
            z[i,j]=np.log10(r[k])
    fig,axes=plt.subplots(2,1,figsize=(9,4.5),gridspec_kw={'height_ratios':[2,1]})
    im=axes[0].pcolormesh(np.r_[grid-.0025,grid[-1]+.0025],[-.5,.5,1.5],
                          np.ma.masked_invalid(z),cmap='viridis',shading='flat')
    axes[0].set_xlim(1.,1.3)
    axes[0].set_title('White band: amplitudes not sampled',fontsize=10)
    axes[0].set_yticks([0,1],['Wall Mach','Pressure / bulk-modulus proxy'])
    axes[0].set_xlabel('Acoustic amplitude (bar)');fig.colorbar(im,ax=axes[0],label='log10 diagnostic')
    axes[1].axis('off');axes[1].text(.02,.5,'Unquantified for all amplitudes: shape stability, gas pressure uniformity,\nthermal and mass transfer, finite molecular volume, continuum validity.\nA converged ODE trajectory does not resolve these uncertainties.',va='center',fontsize=10)
    save(fig,'12_diagnostic_map')
    cl=numeric('c_limit.csv')+numeric('c_limit_extended.csv'); cv=numeric('convergence.csv')
    fig,axes=plt.subplots(1,3,figsize=(11,3.7))
    for bar in (1.,1.3):
        s=[r for r in cl if r['Pa_bar']==bar]
        axes[0].loglog([r['c_factor'] for r in s],[r['relative_Rmin_error'] for r in s],'.-',label=f'{bar:g} bar')
        axes[1].loglog([r['c_factor'] for r in s],[r['max_R_difference_over_R0'] for r in s],'.-',label=f'{bar:g} bar')
    for ax,label in zip(axes[:2],('Relative minimum-radius error','Max trajectory difference / R0')):
        ax.set(xlabel='Sound-speed multiplier',ylabel=label);ax.legend()
    for bar in (1.,1.1,1.15,1.3):
        s=[r for r in cv if r['model']=='KM' and r['Pa_bar']==bar and r['method']=='DOP853' and r['rtol']!=1e-11]
        axes[2].loglog([r['rtol'] for r in s],[max(r['relative_error_Rmin_m'],1e-16) for r in s],'.-',label=f'{bar:g} bar')
    axes[2].set(xlabel='DOP853 relative tolerance',ylabel='Relative Rmin error (KM)');axes[2].legend(fontsize=8)
    save(fig,'13_validation')


def table(name,headers,rows,align=None):
    lines=['\\begin{tabular}{'+(align or 'l'*len(headers))+'}','\\toprule',
           ' & '.join(headers)+r'\\\midrule']
    lines += [' & '.join(map(str,r))+r'\\' for r in rows]
    lines += ['\\bottomrule','\\end{tabular}']
    (RESULTS/name).write_text('\n'.join(lines)+'\n')


def tables():
    rows=numeric('amplitude_summary.csv'); cv=numeric('convergence.csv'); cl=numeric('c_limit.csv')+numeric('c_limit_extended.csv')
    def fmt(x): return f'{x:.6g}' if np.isfinite(x) else '--'
    from .model import Parameters
    p=Parameters()
    provenance=[
        ('R0','$R_0$',p.R0,'m','L','Sadighi-Bonabi 2016, Table II, water'),
        ('f','$f$',p.f,'Hz','L','Sadighi-Bonabi 2016, Table II, water'),
        ('T_liquid','$T_{liquid,0}$',293.15,'K','L','Ambient state; not a computed bubble temperature'),
        ('p0','$p_0$',p.p0,'Pa','A','Standard atmosphere adopted in Mini-Research 1'),
        ('rho','$\\rho$',p.rho,'kg/m$^3$','P','IAPWS-95; accepted rounded reference'),
        ('mu','$\\mu$',p.mu,'Pa s','P','IAPWS viscosity; accepted rounded reference'),
        ('sigma','$\\sigma$',p.sigma,'N/m','P','IAPWS surface tension; accepted rounded reference'),
        ('pv','$p_v$',p.pv,'Pa','P','IAPWS saturation; fixed partial pressure'),
        ('kappa','$\\kappa$',p.kappa,'1','A','Constant diatomic-gas polytropic approximation'),
        ('c','$c$',p.c,'m/s','P','IAPWS-95; accepted rounded reference'),
        ('pg0','$p_{g0}$',p.pg0,'Pa','D','Static gas, vapor and capillary balance'),
        ('T','$T$',p.T,'s','D','Reciprocal driving frequency')]
    table('table_parameters.tex',['Parameter','Value','Unit','Status'],
          [[label,fmt(value),unit,status] for key,label,value,unit,status,source in provenance], 'lrll')
    from .experiments import save_json
    meanings={'R0':'Equilibrium bubble radius','f':'Acoustic driving frequency',
              'T_liquid':'Initial ambient liquid temperature','p0':'Static far-field absolute pressure',
              'rho':'Reference liquid density','mu':'Dynamic liquid viscosity',
              'sigma':'Surface tension at the gas-liquid interface','pv':'Fixed vapor partial pressure',
              'kappa':'Constant gas polytropic exponent','c':'Reference liquid sound speed',
              'pg0':'Initial noncondensable gas partial pressure','T':'Acoustic period'}
    save_json(RESULTS/'parameter_provenance.json',[
        dict(parameter=key,value=value,unit=unit.replace('$^3$','^3'),status=status,source=source,
             meaning=meanings[key]) for key,label,value,unit,status,source in provenance]+[
        dict(parameter='Pa',value=sorted({r['Pa_bar']*1e5 for r in rows}),unit='Pa',status='A',
             source='Prescribed six reference cases and behavior-guided finite mesh',
             meaning='Acoustic pressure amplitude; values are scenarios, not measurements')])
    main=[]
    for bar in CANONICAL:
        for m in ('RP','KM'):
            r=next(r for r in rows if r['model']==m and r['Pa_bar']==bar)
            main.append([f'{bar:g}',m,*[fmt(r[k]*s) for k,s in
                         [('Rmax_m',1e6),('Rmin_m',1e6),('max_abs_U_m_s',1),('max_wall_Mach',1)]],
                         str(int(r['local_minima_count']))])
    table('table_main.tex',['$P_a$ (bar)','Model','$R_{max}$ ($\\mu$m)','$R_{min}$ ($\\mu$m)',
                            '$\\max|U|$ (m/s)','$M_{max}$','$N_{min}$'],main,'rlrrrrr')
    timings=[]
    for bar in (1.,1.1,1.15,1.3):
        for m in ('RP','KM'):
            r=next(r for r in rows if r['model']==m and r['Pa_bar']==bar)
            timings.append([f'{bar:g}',m,fmt(r['collapse_time_s']*1e6),
                            str(int(r['characteristic_minimum_cycle'])),fmt(r['collapse_duration_s']*1e6),
                            fmt(r['first_strong_time_s']*1e6),fmt(r['Vmax_over_Vmin'])])
    table('table_timing.tex',['bar','Model','$t_c$ ($\\mu$s)','Cycle','$\\tau_c$ ($\\mu$s)',
                              'First strong ($\\mu$s)','$V_{max}/V_{min}$'],timings,'rlrrrrr')
    conv=[]
    for r in cv:
        if r['method']=='Radau' or r['rtol']==1e-10:
            conv.append([r['model'],fmt(r['Pa_bar']),r['method'],fmt(r['relative_error_Rmin_m']),
                         fmt(r['relative_error_max_abs_U_m_s']),fmt(r['collapse_time_error_s']*1e9)])
    table('table_convergence.tex',['Model','bar','Method','$\\epsilon(R_{min})$','$\\epsilon(\\max|U|)$',
                                   '$\\Delta t_c$ (ns)'],conv,'lrlrrr')
    # Every displayed numeric result is computed here from saved data.
    km=next(r for r in rows if r['model']=='KM' and r['Pa_bar']==1.3)
    rp=next(r for r in rows if r['model']=='RP' and r['Pa_bar']==1.3)
    km1=next(r for r in rows if r['model']=='KM' and r['Pa_bar']==1.)
    rp1=next(r for r in rows if r['model']=='RP' and r['Pa_bar']==1.)
    complete=[r for r in rows if r['model']=='KM' and r['completed']]
    extreme=max(complete,key=lambda r:r['max_abs_U_m_s'])
    vals={
        'ExtremeBar':extreme['Pa_bar'],'ExtremeRminUm':extreme['Rmin_m']*1e6,
        'ExtremeRmaxUm':extreme['Rmax_m']*1e6,'ExtremeSpeed':extreme['max_abs_U_m_s'],
        'ExtremeMach':extreme['max_wall_Mach'],'ExtremeVolume':extreme['Vmax_over_Vmin'],
        'ExtremeGasGPa':extreme['max_gas_pressure_Pa']/1e9,
        'ExtremeChi':extreme['max_liquid_pressure_change_over_bulk_modulus'],
        'ExtremeD':extreme['min_D_over_R'],'ExtremeTimeUs':extreme['collapse_time_s']*1e6,
        'ExtremeDurationUs':extreme['collapse_duration_s']*1e6,
        'RadiusFactor':km['Rmin_m']/rp['Rmin_m'],
        'SpeedReduction':100*(1-km['max_abs_U_m_s']/rp['max_abs_U_m_s']),
        'OneRadiusDifference':100*(km1['Rmin_m']/rp1['Rmin_m']-1),
        'AmplitudeCount':len(complete),'IncompleteRP':sum(not r['completed'] for r in rows if r['model']=='RP'),
        'KMConvergenceR':max(r['relative_error_Rmin_m'] for r in cv if r['model']=='KM' and r['rtol']<=1e-10),
        'KMConvergenceU':max(r['relative_error_max_abs_U_m_s'] for r in cv if r['model']=='KM' and r['rtol']<=1e-10),
        'RPLimitError':next(r['relative_Rmin_error'] for r in cl if r['Pa_bar']==1.3 and r['c_factor']==1e7),
        'RPLimitTrajectory':cl[-1]['max_R_difference_over_R0'],
        'RPLimitSpeedTrajectory':cl[-1]['max_U_difference_over_reference_max'],
        'PgZero':101325+2*.07274/6e-6-2339,
    }
    (RESULTS/'report_values.tex').write_text('\n'.join(
        '\\newcommand{\\'+k+'}{'+fmt(float(v))+'}' for k,v in vals.items())+'\n')
    summary={k:v for k,v in vals.items()}
    (RESULTS/'report_values.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    tables(); figures()
