"""English figures and report numbers generated only from saved results."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch
from .bubble_thermo import ROOT,MR2,Parameters,adiabatic_temperature

OUT=ROOT/'results';FIG=ROOT/'figures'
COLORS=['#1764a0','#b65313','#207344']
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.alpha':.22,'figure.facecolor':'white','savefig.facecolor':'white',
    'font.family':'DejaVu Sans','lines.linewidth':1.25})


def load():
    return pd.read_csv(OUT/'adiabatic_summary.csv'),pd.read_csv(OUT/'thermal_summary.csv')


def save(fig,name):
    fig.savefig(FIG/(name+'.jpg'),dpi=300,bbox_inches='tight',pil_kwargs={'quality':95})
    plt.close(fig)


def thermal_series(row):
    return np.load(OUT/(row.run_id+'.npz'))


def get(b,amp,eta=1,gas='air_like',mode='B_fully_coupled'):
    return b[(b.Pa_bar==amp)&(b.eta==eta)&(b.gas==gas)&(b.thermodynamics==mode)].iloc[0]


def schematic():
    fig,ax=plt.subplots(figsize=(8,2.35));ax.axis('off');ax.set_xlim(-3.9,3.9);ax.set_ylim(-1.3,1.3)
    ax.add_patch(Circle((0,0),.88,facecolor='#f4f4f4',edgecolor='black',lw=1.3))
    ax.text(0,.12,'Fixed noncondensable gas\n$U_g=nC_VT$\n$p_gV=nR_uT$',ha='center',va='center',fontsize=7)
    ax.text(0,-.52,r'$V=4\pi R^3/3$',ha='center',fontsize=8)
    for a,z in [((-3.,0),(-.91,0)),((.91,0),(3.,0))]:
        ax.add_patch(FancyArrowPatch(a,z,arrowstyle='<->',mutation_scale=14,color='black'))
    ax.text(-2.55,.3,'Compression work\n'+r'$-p_g\,dV/dt$',ha='center',va='bottom')
    ax.text(2.55,.3,'Signed heat transfer\n'+r'$4\pi\eta k_gR(T-T_\ell)$',ha='center',va='bottom')
    ax.text(2.55,-.55,r'Liquid reservoir: $T_\ell$ fixed',ha='center')
    ax.text(-2.55,-.55,'External heat source: 0 W',ha='center')
    ax.text(0,1.12,'Spherical bubble in water; vapor pressure enters mechanics separately',ha='center')
    save(fig,'01_energy_balance')


def histories(a,b):
    fig,axes=plt.subplots(3,2,figsize=(8,6.0),sharex=True)
    for i,amp in enumerate([1.,1.15,1.3]):
        ow=thermal_series(get(b,amp,mode='B_one_way_post_processing'));co=thermal_series(get(b,amp))
        axes[i,0].plot(ow['t_s']*1e6,ow['R_m']*1e6,color=COLORS[0],label='Saved KM radius')
        axes[i,0].plot(co['t_s']*1e6,co['R_m']*1e6,color=COLORS[2],ls='--',label='Coupled radius, eta=1')
        axes[i,1].plot(ow['t_s']*1e6,adiabatic_temperature(ow['R_m']),color=COLORS[0],label='A: adiabatic')
        axes[i,1].plot(ow['t_s']*1e6,ow['T_K'],color=COLORS[1],label='B: one-way, eta=1')
        axes[i,1].plot(co['t_s']*1e6,co['T_K'],color=COLORS[2],ls='--',label='B: coupled, eta=1')
        axes[i,0].set_ylabel(f'{amp:g} bar\nRadius (um)');axes[i,1].set_ylabel('Gas temperature (K)')
        axes[i,1].set_yscale('log')
    for ax in axes[-1]:ax.set_xlabel('Time (us)')
    for ax in axes[0]:ax.legend(fontsize=7,loc='upper right')
    fig.tight_layout();save(fig,'02_radius_temperature')


def zoom(a,b):
    fig,axes=plt.subplots(1,3,figsize=(8.1,3.15))
    for ax,amp in zip(axes,[1.,1.15,1.3]):
        owrow=get(b,amp,mode='B_one_way_post_processing');corow=get(b,amp)
        ow=thermal_series(owrow);co=thermal_series(corow)
        ar=a[(a.mechanics=='KM')&(a.Pa_bar==amp)].iloc[0]
        for z,T,tc,label,color,style in [(ow,adiabatic_temperature(ow['R_m']),ar.t_Tmax_s,'A',COLORS[0],'-'),
                    (ow,ow['T_K'],owrow.t_Tmax_s,'One-way',COLORS[1],'-'),
                    (co,co['T_K'],corow.t_Tmax_s,'Coupled',COLORS[2],'--')]:
            dt=(z['t_s']-tc)*1e9;mask=abs(dt)<(150 if amp==1 else 15)
            ax.plot(dt[mask],T[mask],label=f'{label}: {tc*1e6:.4f} us',color=color,ls=style)
        ax.set_title(f'{amp:g} bar');ax.set_xlabel('Time relative to each T peak (ns)')
        ax.legend(fontsize=6.5,loc='upper center',bbox_to_anchor=(.5,-.28),frameon=False)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    axes[0].set_ylabel('Gas temperature (K)');fig.tight_layout();save(fig,'03_peak_zoom')


def sensitivity(a,b):
    fig,axes=plt.subplots(1,2,figsize=(8,3.3))
    for model,color in [('KM',COLORS[0]),('RP','#666666')]:
        s=a[(a.mechanics==model)].sort_values('Pa_bar')
        y=s.global_Tmax_K.to_numpy()
        axes[0].plot(s.Pa_bar,y,color=color,label=model+(' diagnostic' if model=='RP' else ' adiabatic'))
        if model=='RP':
            bad=s[~s.completed];axes[0].scatter(bad.Pa_bar,bad.observed_Tmax_K,marker='x',color='#aa3333',s=18,label='RP partial: observed only')
    for amp,marker in [(1.,'o'),(1.15,'s'),(1.3,'^')]:
        for mode,color,ls in [('B_one_way_post_processing',COLORS[1],'-'),('B_fully_coupled',COLORS[2],'--')]:
            s=b[(b.Pa_bar==amp)&(b.gas=='air_like')&(b.thermodynamics==mode)&(b.eta>0)].sort_values('eta')
            axes[1].plot(s.eta,s.observed_Tmax_K,marker=marker,color=color,ls=ls,label=f'{amp:g} bar, '+('one-way' if 'one_way' in mode else 'coupled'))
    axes[0].set_yscale('log');axes[0].set_xlabel('Acoustic amplitude (bar)');axes[0].set_ylabel('Peak gas temperature (K)')
    axes[0].axvspan(1.2,1.3,alpha=.08,color='black');axes[0].legend(fontsize=6.5)
    axes[1].set_xscale('log');axes[1].set_yscale('log');axes[1].set_xlabel('Heat-transfer efficiency eta');axes[1].set_ylabel('Peak gas temperature (K)')
    axes[1].legend(fontsize=6.5,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.02),frameon=False)
    fig.tight_layout();save(fig,'04_amplitude_heat_sensitivity')
    s=pd.read_csv(OUT/'adiabatic_sensitivity.csv');fig,axes=plt.subplots(1,2,figsize=(8,2.8))
    for amp,color in zip([1.,1.15,1.3],COLORS):
        z=s[(s.Pa_bar==amp)&(s.T0_K==293.15)].sort_values('gamma')
        axes[0].plot(z.gamma,z.observed_Tmax_K,'o-',color=color,label=f'{amp:g} bar')
    for gamma,color in zip([1.2,1.4,1.55,5/3],['#666666',*COLORS]):
        z=s[(s.Pa_bar==1.3)&np.isclose(s.gamma,gamma)].sort_values('T0_K')
        axes[1].plot(z.T0_K,z.observed_Tmax_K,'o-',color=color,label=f'gamma={gamma:.3g}')
    for ax in axes:ax.set_yscale('log');ax.set_ylabel('Adiabatic peak temperature (K)');ax.legend(fontsize=7)
    axes[0].set_xlabel('Constant gamma');axes[1].set_xlabel('Initial gas temperature (K)')
    axes[0].set_title('Fixed accepted KM trajectories');axes[1].set_title('1.3 bar; fixed mechanics and water')
    fig.tight_layout();save(fig,'05_gamma_initial_temperature')


def comparison(a,b):
    events=pd.read_csv(OUT/'collapse_temperatures.csv');fig,axes=plt.subplots(1,2,figsize=(8,2.8))
    for model,color,marker in [('RP','#666666','s'),('KM',COLORS[0],'o')]:
        s=events[(events.mechanics==model)&events.first_strong].sort_values('Pa_bar')
        axes[0].plot(s.Pa_bar,s.T_K,marker+'-',color=color,label=model+' first strong')
        s=a[(a.mechanics==model)&a.Pa_bar.isin([1.,1.15,1.3])].sort_values('Pa_bar')
        axes[0].plot(s.Pa_bar,s.observed_Tmax_K,marker+'--',color=color,label=model+' global')
    axes[0].set_yscale('log');axes[0].set_xlabel('Acoustic amplitude (bar)');axes[0].set_ylabel('Adiabatic temperature (K)');axes[0].legend(fontsize=7)
    axes[0].set_title('Event selections kept separate')
    for amp,color in zip([1.,1.15,1.3],COLORS):
        row=get(b,amp,mode='B_one_way_post_processing');z=thermal_series(row)
        pgold=Parameters().pg0*(Parameters().R0/z['R_m'])**4.2
        axes[1].plot(z['t_s']*1e6,z['pg_Pa']/pgold,color=color,label=f'{amp:g} bar')
    axes[1].axhline(1,color='black',ls=':',lw=.8);axes[1].set_xlabel('Time (us)');axes[1].set_ylabel('One-way pressure / accepted pressure');axes[1].set_yscale('log');axes[1].legend(fontsize=7)
    axes[1].set_yticks([1,2,4,8],labels=['1','2','4','8'])
    fig.tight_layout();save(fig,'06_event_comparison_pressure')


def validation_figure():
    s=pd.read_csv(OUT/'solver_convergence.csv');q=pd.read_csv(OUT/'input_grid_convergence.csv')
    fig,axes=plt.subplots(1,2,figsize=(8,2.65))
    for label,color in [('coarse','#777777'),('production',COLORS[0]),('independent_Radau',COLORS[1])]:
        z=s[(s.eta==1)&(s.gas=='air_like')&(s.configuration==label)].sort_values('Pa_bar')
        axes[0].plot(z.Pa_bar,np.maximum(z.relative_Tmax_error,1e-14),'o-',label=label.replace('_',' '),color=color)
    for label,color in [('export_stride4','#777777'),('export_stride2',COLORS[0]),('accepted_stride1',COLORS[1])]:
        z=q[(q.thermodynamics=='B_one_way_post_processing')&(q.grid==label)].sort_values('Pa_bar')
        axes[1].plot(z.Pa_bar,np.maximum(z.relative_Tmax_error,1e-14),'o-',label=label.replace('_',' '),color=color)
    for ax in axes:ax.set_yscale('log');ax.set_xlabel('Acoustic amplitude (bar)');ax.set_ylabel('Relative peak-temperature difference');ax.legend(fontsize=7)
    axes[0].set_title('Coupled: tighter DOP853 reference');axes[1].set_title('One-way: full saved export reference')
    fig.tight_layout();save(fig,'07_numerical_checks')


def tex_table(file,headers,rows,align=None):
    text='\\begin{tabular}{'+(align or 'r'*len(headers))+'}\n\\toprule\n'+' & '.join(headers)+' \\\\\n\\midrule\n'
    text+='\n'.join(' & '.join(str(v) for v in row)+' \\\\' for row in rows)
    text+='\n\\bottomrule\n\\end{tabular}\n';(OUT/file).write_text(text)


def tables(a,b):
    rows=[];combined=[]
    for amp in [1.,1.15,1.3]:
        ar=a[(a.mechanics=='KM')&(a.Pa_bar==amp)].iloc[0];ow=get(b,amp,mode='B_one_way_post_processing');co=get(b,amp)
        rows.append([f'{amp:.2f}',f'{ar.observed_Tmax_K:.2f}',f'{ow.observed_Tmax_K:.2f}',f'{co.observed_Tmax_K:.2f}'])
        combined.append(dict(Pa_bar=amp,A_K=ar.observed_Tmax_K,B_one_way_K=ow.observed_Tmax_K,
                            B_fully_coupled_K=co.observed_Tmax_K,category=ar.upstream_next_stage_use))
    tex_table('table_main.tex',['$P_a$ (bar)','A: adiabatic (K)','B: one-way (K)','B: coupled (K)'],rows)
    pd.DataFrame(combined).to_csv(OUT/'main_temperature_comparison.csv',index=False)
    rows=[]
    for amp in [1.,1.15,1.3]:
        for mode,label in [('B_one_way_post_processing','one-way'),('B_fully_coupled','coupled')]:
            r=get(b,amp,mode=mode)
            rows.append([f'{amp:.2f}',label,f'{r.t_Tmax_s*1e6:.6f}',f'{r.Tmax_minus_Rmin_time_s*1e9:.3g}',
                f'{r.duration_above_1000K_s*1e9:.5g}',f'{r.duration_above_10000K_s*1e9:.5g}',
                f'{r.max_wall_Mach:.3g}'])
    tex_table('table_peaks.tex',['bar','Model','$t_T$ ($\\mu$s)','$t_T-t_R$ (ns)','$>10^3$ K (ns)','$>10^4$ K (ns)','$M_{\\max}$'],rows,align='rlrrrrr')
    rows=[]
    for amp in [1.,1.15,1.3]:
        for gas in ['air_like','argon_like']:
            r=get(b,amp,gas=gas);w=get(b,amp,gas=gas,mode='B_one_way_post_processing')
            rows.append([f'{amp:.2f}',gas.replace('_like',''),f'{w.observed_Tmax_K:.2f}',f'{r.observed_Tmax_K:.2f}',
                f'{r.pg_at_Tmax_Pa/1e9:.4g}',f'{r.density_at_Tmax_kg_m3:.4g}'])
    tex_table('table_gases.tex',['bar','Ideal gas','One-way (K)','Coupled (K)','$p_g(t_T)$ (GPa)','$\\rho_g(t_T)$ (kg/m$^3$)'],rows,align='rlrrrr')
    rows=[];events=pd.read_csv(OUT/'collapse_temperatures.csv')
    for amp in [1.,1.15,1.3]:
        for model in ['RP','KM']:
            r=a[(a.Pa_bar==amp)&(a.mechanics==model)].iloc[0]
            first=events[(events.Pa_bar==amp)&(events.mechanics==model)&events.first_strong]
            rows.append([f'{amp:.2f}',model,'--' if first.empty else f'{first.iloc[0].T_K:.2f}',
                f'{r.observed_Tmax_K:.2f}',str(r.peak_cycle),f'{r.max_wall_Mach:.4g}'])
    tex_table('table_rp_km.tex',['bar','Model','First strong (K)','Global (K)','Global cycle','$M_{\\max}$'],rows,align='rlrrrr')
    c=pd.read_csv(OUT/'solver_convergence.csv');q=pd.read_csv(OUT/'input_grid_convergence.csv')
    rows=[]
    for amp in [1.,1.15,1.3]:
        for typ,df,mask in [('coupled',c,(c.configuration=='production')&(c.eta==1)&(c.gas=='air_like')),
                           ('one-way',q,(q.grid=='full_export_production')&(q.eta==1)&(q.gas=='air_like')&(q.thermodynamics=='B_one_way_post_processing'))]:
            r=df[mask&(df.Pa_bar==amp)].iloc[0]
            rows.append([f'{amp:.2f}',typ,f'{r.relative_Tmax_error:.2e}',f'{r.abs_peak_time_error_s:.2e}',f'{r.abs_duration_1000K_error_s:.2e}'])
    tex_table('table_convergence.tex',['bar','Model','Relative $T_{\\max}$','$|\\delta t_T|$ (s)','$|\\delta t_{>1000}|$ (s)'],rows,align='rlrrr')
    macro={}
    def put(k,v,fmt='.3g'):macro[k]=format(v,fmt)
    extreme=a[(a.mechanics=='KM')&(a.Pa_bar==1.3)].iloc[0]
    for k,v in [('AirA',extreme.observed_Tmax_K),('AirOW',get(b,1.3,mode='B_one_way_post_processing').observed_Tmax_K),
                ('AirCO',get(b,1.3).observed_Tmax_K),('ArgOW',get(b,1.3,gas='argon_like',mode='B_one_way_post_processing').observed_Tmax_K),
                ('ArgCO',get(b,1.3,gas='argon_like').observed_Tmax_K)]:put(k,v,'.2f')
    put('Nmol',Gas_n(),'.8g')
    put('TauInitialUs',get(b,1.3).tau_thermal_at_R0_s*1e6,'.3f')
    put('ExtremeMach',extreme.max_wall_Mach,'.3f')
    put('ExtremeChi',extreme.pressure_bulk_modulus_proxy,'.3f')
    put('MaxResidual',b.energy_residual_scaled.max(),'.2e')
    put('MaxCoupledError',c[c.configuration=='production'].relative_Tmax_error.max(),'.2e')
    put('MaxIndependentError',c[c.configuration=='independent_Radau'].relative_Tmax_error.max(),'.2e')
    fine=q[(q.grid=='export_stride2')&(q.thermodynamics=='B_one_way_post_processing')]
    put('MaxInputError',fine.relative_Tmax_error.max(),'.2e')
    put('MaxInputDurationError',fine.relative_duration_1000K_error.max(),'.2e')
    put('MaxRecoveryError',pd.read_csv(OUT/'adiabatic_recovery.csv').relative_Tmax_vs_saved_mechanics.max(),'.2e')
    put('AdiabaticConsistency',a.refined_vs_interpolant_relative_error.max(),'.2e')
    macro.update(AmplitudeCount=str(a.Pa_bar.nunique()),PartialCount=str(sum(~a.completed)),ThermalCount=str(len(b)))
    # Scientific notation is TeX-safe and all scientific values come from saved CSVs.
    def texnum(s):
        if 'e' in s:
            x,e=s.split('e');return x+'\\times10^{'+str(int(e))+'}'
        return s
    (OUT/'report_values.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+texnum(v)+'}' for k,v in macro.items())+'\n')


def Gas_n():
    return json.loads((OUT/'parameter_provenance.json').read_text())['n_mol']


def main():
    FIG.mkdir(exist_ok=True);a,b=load();schematic();histories(a,b);zoom(a,b);sensitivity(a,b);comparison(a,b)
    validation_figure();tables(a,b)


if __name__=='__main__':main()
