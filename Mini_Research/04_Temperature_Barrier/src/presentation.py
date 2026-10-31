"""Tables and print-friendly JPG figures from saved results only."""
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .upstream import ROOT,MR3,Parameters,RU
from .thermodynamics import Mixture,species_cp
from .run_experiments import run_key

OUT=ROOT/'results';FIG=ROOT/'figures'
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.22,'figure.facecolor':'white',
                     'axes.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
COLORS={'MR3':'#4d4d4d','M1_variable':'#1565a7','M2_sat':'#b34e18'}

def latex_table(path,headers,rows,align=None):
    align=align or ('l'*len(headers))
    text='\\begin{tabular}{'+align+'}\\toprule\n'+' & '.join(headers)+' \\\\\n\\midrule\n'
    text+=''.join(' & '.join(map(str,r))+' \\\\\n' for r in rows)
    text+='\\bottomrule\\end{tabular}\n';(OUT/path).write_text(text)

def save(fig,name):
    fig.savefig(FIG/(name+'.jpg'),dpi=300,bbox_inches='tight');plt.close(fig)

def main():
    df=pd.read_csv(OUT/'summary.csv',float_precision='round_trip');old=pd.read_csv(MR3/'results/thermal_summary.csv',float_precision='round_trip')
    def get(amp,gas,mode,c,eta=1):
        return df[df.run_id==run_key(amp,gas,mode,c,eta)].iloc[0]
    def oldget(amp,gas,mode):
        return old[(old.Pa_bar==amp)&(old.gas==gas+'_like')&(old.eta==1)&
                   (old.thermodynamics==('B_fully_coupled' if mode=='coupled' else 'B_one_way_post_processing'))].iloc[0]
    rows=[];table=[]
    for amp in (1.,1.15,1.3):
        for gas in ('air','argon'):
            for mode in ('one_way','coupled'):
                a=oldget(amp,gas,mode);b=get(amp,gas,mode,'M1_variable');c=get(amp,gas,mode,'M2_sat')
                rows.append(dict(Pa_bar=amp,gas=gas,mode=mode,MR3_Tmax_K=a.observed_Tmax_K,
                    M1_observed_Tmax_K=b.observed_Tmax_K,M1_completed=b.completed,
                    MR4_observed_Tmax_K=c.observed_Tmax_K,MR4_completed=c.completed,
                    MR4_global_Tmax_K=c.global_Tmax_K,MR4_run_id=c.run_id,
                    delta_T_MR4_vs_MR3=c.observed_Tmax_K/a.observed_Tmax_K-1 if c.completed else math.nan))
                mark=lambda r:f'{r.observed_Tmax_K:.1f}'+('' if r.completed else '$^\\dagger$')
                table.append([f'{amp:.2f}','Air' if gas=='air' else 'Ar','OW' if mode=='one_way' else 'CO',
                              f'{a.observed_Tmax_K:.1f}',mark(b),mark(c),
                              f'{100*(c.observed_Tmax_K/a.observed_Tmax_K-1):+.1f}' if c.completed else '--'])
    pd.DataFrame(rows).to_csv(OUT/'main_comparison.csv',index=False)
    latex_table('table_main.tex',[r'$P_a$ (bar)','Gas','Mode','MR3 (K)','M1 (K)','M2 (K)',r'$100\delta_T$'],table,'rllrrrr')
    # Explicitly paired contrasts, same gas, mode, amplitude and conductivity.
    sensitiv=[]
    comparisons=[('variable_Cv_legacy','M1_variable','M1_const'),('vapor_pressure_convention','M2_dry','M1_variable'),
                 ('saturated_vapor','M2_sat','M2_dry'),('variable_Cv_wet','M2_sat','M2_sat_const'),
                 ('trapped_vapor_5','M2_v5','M2_dry'),('trapped_vapor_10','M2_v10','M2_dry')]
    for mechanism,new,ref in comparisons:
        for r in df[(df.config==new)&(df.eta==1)].to_dict('records'):
            b=get(r['Pa_bar'],r['gas'],r['mode'],ref);valid=r['completed'] and b.completed
            sensitiv.append(dict(mechanism=mechanism,new_run_id=r['run_id'],reference_run_id=b.run_id,
                new_observed_Tmax_K=r['observed_Tmax_K'],reference_observed_Tmax_K=b.observed_Tmax_K,
                both_five_cycles=valid,delta_T=r['observed_Tmax_K']/b.observed_Tmax_K-1 if valid else math.nan,
                observed_max_ratio_diagnostic_only=r['observed_Tmax_K']/b.observed_Tmax_K-1))
    for r in df[df.eta!=1].to_dict('records'):
        b=get(r['Pa_bar'],r['gas'],r['mode'],'M2_sat');valid=r['completed'] and b.completed
        sensitiv.append(dict(mechanism='heat_transfer',new_run_id=r['run_id'],reference_run_id=b.run_id,
            new_observed_Tmax_K=r['observed_Tmax_K'],reference_observed_Tmax_K=b.observed_Tmax_K,
            both_five_cycles=valid,delta_T=r['observed_Tmax_K']/b.observed_Tmax_K-1 if valid else math.nan,
            observed_max_ratio_diagnostic_only=r['observed_Tmax_K']/b.observed_Tmax_K-1))
    sd=pd.DataFrame(sensitiv);sd.to_csv(OUT/'sensitivity.csv',index=False)
    energy=df[(df.config=='M2_sat')&(df.eta==1)].copy();energy.to_csv(OUT/'energy_budgets.csv',index=False)
    erows=[]
    for r in energy.to_dict('records'):
        erows.append([f"{r['Pa_bar']:.2f}",r['gas'],'OW' if r['mode']=='one_way' else 'CO',
          *[f"{r[k]*1e9:.4f}" for k in ('delta_U_peak_J','work_by_gas_peak_J','signed_heat_loss_peak_J')],
          f"{r['energy_residual_scaled']:.1e}"])
    latex_table('table_energy.tex',[r'$P_a$','Gas','Mode',r'$\Delta U$ (nJ)',r'$W_{\rm by}$ (nJ)',r'$Q_{\rm loss}$ (nJ)',r'$\epsilon_E$'],erows,'rllrrrr')
    boundary=df[~df.completed][['run_id','end_s','boundary_temperature_K','solver_message','observed_Tmax_K']]
    boundary.to_csv(OUT/'boundary_events.csv',index=False)
    # Capacity is per mole of the entire mixture, not per mole of dry gas.
    fig,axs=plt.subplots(1,2,figsize=(9.5,2.6),layout='constrained');T=np.geomspace(200,6000,400)
    for gas,label,color in [('air','Dry N2/O2','#1565a7'),('argon','Ground-state Ar','#4d4d4d')]:
        mix=Mixture(gas=gas,vapor_factor=0);wet=Mixture(gas=gas)
        axs[0].plot(T,[mix.cv(float(t))/RU for t in T],label=label,color=color)
        axs[1].plot(T,[1+RU/wet.cv(float(t)) for t in T],label=label+' + saturated inventory',color=color)
    axs[0].plot(T,[(species_cp('H2O',float(t))-RU)/RU for t in T],label='H2O',color='#b34e18')
    axs[0].set(ylabel=r'Molar $c_V/R_u$',xlabel='Temperature (K)',xscale='log')
    axs[1].set(ylabel=r'Diagnostic $1+R_u/c_V$',xlabel='Temperature (K)',xscale='log')
    for ax in axs: ax.legend(fontsize=7)
    save(fig,'01_heat_capacity')
    fig,axs=plt.subplots(2,2,figsize=(9.5,4.5),layout='constrained')
    for i,gas in enumerate(('air','argon')):
        for c,label in [('M0','MR3 limit'),('M1_variable','Variable dry heat capacity'),('M2_sat','Fixed saturated inventory')]:
            key=run_key(1.3,gas,'coupled',c,1);r=df[df.run_id==key].iloc[0]
            with np.load(OUT/'trajectories'/f'{key}.npz') as z:
                axs[i,0].plot(z['t_s']*1e6,z['T_K'],label=label,lw=.8)
                axs[i,1].plot((z['t_s']-r.t_Tmax_s)*1e9,z['T_K'],label=label,lw=1)
                if not r.completed:
                    for ax,x in [(axs[i,0],r.end_s*1e6),(axs[i,1],(r.end_s-r.t_Tmax_s)*1e9)]:
                        ax.scatter(x,r.boundary_temperature_K,s=20,marker='x',color='black',zorder=5)
        axs[i,0].set(title=gas.title()+', coupled, 1.3 bar',xlabel='Time (microseconds)',ylabel='Temperature (K)',yscale='log')
        axs[i,1].set(xlabel=r'Time relative to each $T_{\max}$ (ns)',ylabel='Temperature (K)',xlim=(-15,15))
    axs[0,0].legend(fontsize=7);save(fig,'02_temperature_histories')
    fig,axs=plt.subplots(1,2,figsize=(9.5,2.8),layout='constrained')
    for ax,gas in zip(axs,('air','argon')):
        for mode,style in [('one_way','o--'),('coupled','s-')]:
            rr=[get(1.3,gas,mode,c) for c in ['M2_dry','M2_sat','M2_v5','M2_v10']]
            ax.plot([r.vapor_mole_fraction*100 for r in rr],[r.observed_Tmax_K if r.completed else np.nan for r in rr],style,label='One-way' if mode=='one_way' else 'Coupled',lw=1)
            for r in rr:
                if not r.completed: ax.scatter(r.vapor_mole_fraction*100,r.observed_Tmax_K,s=65,marker='^',facecolor='white',edgecolor='black',zorder=4)
        ax.set(title=gas.title()+', 1.3 bar',xlabel='Fixed H2O mole fraction (%)',ylabel='Observed peak temperature (K)');ax.legend(fontsize=7)
        if gas=='argon': ax.set_yscale('log')
    save(fig,'03_vapor_sensitivity')
    fig,axs=plt.subplots(1,2,figsize=(9.5,2.8),layout='constrained')
    for ax,(amp,gas) in zip(axs,[(1.15,'air'),(1.3,'argon')]):
        r=get(amp,gas,'coupled','M2_sat')
        with np.load(OUT/'trajectories'/f'{r.run_id}.npz') as z:
            for field,label in [('delta_U_J',r'$\Delta U$'),('work_by_gas_J',r'$-W_{\rm by}$'),('signed_heat_loss_J',r'$-Q_{\rm loss}$')]:
                sign=1 if field=='delta_U_J' else -1
                ax.plot(z['t_s']*1e6,sign*z[field]*1e9,label=label,lw=.9)
        ax.set(title=f'{gas.title()}, {amp:g} bar, coupled',xlabel='Time (microseconds)',ylabel='Energy (nJ)');ax.legend(fontsize=7)
    save(fig,'04_energy_budgets')
    fig,axs=plt.subplots(2,2,figsize=(9.5,4.2),layout='constrained');x=np.arange(3)
    for i,gas in enumerate(('air','argon')):
        for j,mode in enumerate(('one_way','coupled')):
            ax=axs[i,j]
            for k,(c,label) in enumerate([('M0','MR3'),('M1_variable','M1'),('M2_sat','M2')]):
                rr=[get(a,gas,mode,c) for a in (1.,1.15,1.3)]
                bars=ax.bar(x+(k-1)*.25,[r.observed_Tmax_K for r in rr],width=.23,label=label,color=list(COLORS.values())[k])
                for bar,r in zip(bars,rr):
                    if not r.completed: bar.set_hatch('///');bar.set_edgecolor('black')
            ax.set(xticks=x,xticklabels=['1.0','1.15','1.3'],xlabel='Acoustic amplitude (bar)',ylabel='Observed peak (K)',yscale='log',title=gas.title()+', '+mode.replace('_',' '))
    axs[0,0].legend(fontsize=7,ncol=3);save(fig,'05_model_comparison')
    # All report numbers are emitted from saved tables.
    conv=pd.read_csv(OUT/'convergence.csv',float_precision='round_trip');reg=pd.read_csv(OUT/'mr3_regression.csv',float_precision='round_trip');qq=pd.read_csv(OUT/'independent_quadrature.csv',float_precision='round_trip')
    macros={'RunCount':len(df),'CompletedCount':int(df.completed.sum()),'CensoredCount':int((~df.completed).sum()),
        'MaximumResidual':f'{df.energy_residual_scaled.max():.2e}','MaximumRegression':f'{abs(reg.relative_difference).max():.2e}',
        'MaximumStrictDifference':f"{abs(conv[conv.check=='strict'].peak_relative_difference).max():.2e}",
        'MaximumRadauDifference':f"{abs(conv[conv.check=='Radau'].peak_relative_difference).max():.2e}",
        'MaximumInputDifference':f"{abs(conv[conv.check.str.startswith('input_')].peak_relative_difference).max():.2e}",
        'MaximumQuadratureResidual':f'{qq[qq.order==16].scaled_residual.max():.2e}',
        'SaturatedFraction':f'{Mixture().fractions["H2O"]*100:.4f}'}
    for gas,prefix in [('air','Air'),('argon','Argon')]:
        for mode,suffix in [('one_way','OW'),('coupled','CO')]:
            for c,tag in [('M1_variable','Dry'),('M2_sat','Wet')]: macros[prefix+tag+suffix]=f'{get(1.3,gas,mode,c).observed_Tmax_K:.1f}'
    budgets=pd.read_csv(OUT/'chemistry_inventory_budgets.csv',float_precision='round_trip')
    def chem(gas,process):
        return budgets[(budgets.run_id==run_key(1.3,gas,'coupled','M2_sat',1))&(budgets.process==process)].iloc[0]
    for key,gas,process in [('ChemWaterFirst','air','H2O -> OH + H'),('ChemWaterFull','air','H2O -> 2H + O'),('ChemAirNtwo','air','N2 -> 2N'),('ChemArgonIon','argon','Ar -> Ar+ + e-')]:
        r=chem(gas,process);macros[key+'NJ']=f'{r.energy_for_complete_conversion_J*1e9:.3f}'
        macros[key+'ToU']=f'{r.ratio_to_peak_delta_U:.2f}';macros[key+'ToWork']=f'{r.ratio_to_collapse_compression:.2f}'
    (OUT/'report_values.tex').write_text(''.join('\\newcommand{\\'+k+'}{'+str(v)+'}\n' for k,v in macros.items()))
    # Appendix tables keep all diagnostics machine-readable in CSV, selected columns here.
    vr=[]
    for r in conv.to_dict('records'):
        if r['check']=='Radau' or (r['check']=='strict' and 'M2_sat' in r['run_id']):
            short=r['run_id'].replace('M2_sat_','').replace('bar','').replace('_eta1','').replace('one_way','OW').replace('coupled','CO').replace('_',' / ')
            vr.append([short,r['check'],f"{r['peak_relative_difference']:.2e}",f"{r['peak_time_difference_s']*1e9:.3g}",'full' if r['completed'] else 'censored'])
    latex_table('table_convergence.tex',['Case','Check',r'$\Delta T/T$',r'$\Delta t_T$ (ns)','Scope'],vr,'llrrl')
    cr=pd.read_csv(OUT/'reaction_energy_scales.csv',float_precision='round_trip');chemrows=[]
    for r in cr.to_dict('records'):
        chemrows.append([r['process'].replace('->',r'$\to$').replace('e-',r'$e^-$').replace('+',r'$+$'),
            f"{r['energy_eV']:.3f}",f"{r['delta_U_298_J_mol']/1000:.2f}",f"{r['equivalent_E_over_k_K']:.0f}"])
    latex_table('table_chemistry.tex',['Process','Energy (eV)',r'Cost (kJ/mol)',r'$E/R_u$ (K)'],chemrows,'lrrr')
    params=dict(water=vars(Parameters()),T0_K=293.15,pressure_convention='p_infty = p0 - Pa*sin(2*pi*f*t)',
                gas_composition='79/21 N2/O2 or ground-state Ar; optional fixed H2O',models={c:dict(v) for c,v in __import__('src.run_experiments',fromlist=['CONFIGS']).CONFIGS.items()})
    (OUT/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')

if __name__=='__main__': main()
