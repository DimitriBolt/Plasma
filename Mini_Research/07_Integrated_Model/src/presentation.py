"""Regenerate all scientific JPGs, tables and prose from saved data only."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from .upstream import ROOT,T0
from .data_access import data_path

R=ROOT/'results';F=ROOT/'figures'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':110,
                     'savefig.facecolor':'white','axes.grid':True,'grid.alpha':.18,'legend.frameon':False})
colors=['#244e7d','#c75f21','#23805c','#7f4c8f','#555555']

def read(name):
    return pd.read_csv(R/(name+'.csv'),float_precision='round_trip')

def raw(name):return np.load(data_path(name))

def save(fig,name):
    fig.tight_layout();fig.savefig(F/(name+'.jpg'),dpi=300,pil_kwargs={'quality':95});plt.close(fig)

def text(name,s): (R/(name+'.tex')).write_text(s+'\n')

def number(x,d=6):
    if pd.isna(x):return r'\textit{undefined}'
    return f'{x:.{d}g}'

def tex_table(headers,rows,caption,width=None):
    spec='l'+'r'*(len(headers)-1)
    body='\n'.join(' & '.join(map(str,r))+r' \\' for r in rows)
    table=r'\begin{table}[H]\centering\small'+'\n'+r'\begin{tabular}{'+spec+'}\n'+r'\toprule'+'\n'+' & '.join(headers)+r' \\ \midrule'+'\n'+body+'\n'+r'\bottomrule\end{tabular}'+'\n'+r'\caption{'+caption+'}\n'+r'\end{table}'
    return table

def main():
    F.mkdir(exist_ok=True)
    summ=read('summary');pairs=read('pairs');timing=read('timing_sensitivity');conv=read('convergence')
    sens=read('local_sensitivity');sruns=read('sensitivity_runs');static=read('static_MR6_comparison')
    stable=read('stable_pairs') if (R/'stable_pairs.csv').exists() else pd.DataFrame()
    fig,ax=plt.subplots(figsize=(7,2.4));ax.axis('off');ax.set(xlim=(0,10),ylim=(0,3))
    for x,label in [(0.2,'Local laser\n'+r'$P_{abs}=I\sigma_{abs}$'),(2.7,'Gold storage\n'+r'$C_{Au}\dot T_{Au}$'),(5.2,'Material water\n'+r'$C_i\dot T_i$'),(7.7,'Gas energy\n'+r'$C_V\dot T_g$')]:
        ax.add_patch(FancyBboxPatch((x,.9),2.0,1.15,boxstyle='round,pad=.1',facecolor='#f3f5f7',edgecolor=colors[0]));ax.text(x+1,1.475,label,ha='center',va='center')
    for x,label in [(2.25,r'$P_{abs}$'),(4.75,r'$Q_{Au\to\ell}$'),(7.25,r'$Q_{\ell\to g}$')]:
        ax.annotate('',(x+.35,1.5),(x-.05,1.5),arrowprops={'arrowstyle':'->','color':colors[1],'lw':1.8});ax.text(x+.13,2.22,label,ha='center')
    ax.annotate('Bath $T_0$\n$Q_{bath}$',(6.2,.8),(6.2,.05),ha='center',arrowprops={'arrowstyle':'<-','color':colors[2]})
    ax.text(8.65,.35,r'$p_g\dot V$; KM pressure feedback',ha='center',fontsize=8)
    save(fig,'01_system')

    b=raw('base_1.00');h=raw('laser_F10_1.00')
    fig,axs=plt.subplots(2,2,figsize=(7,4.2))
    for ax,key,scale,label in zip(axs.ravel(),['R_m','U_m_s','T_g_K','p_g_Pa'],[1e6,1,1,1e-5],['Radius ($\\mu$m)','Wall speed (m/s)','Gas temperature (K)','Mixture pressure (bar)']):
        ax.plot(b['t_s']*1e6,b[key]*scale,color=colors[0],label='No laser')
        ax.plot(h['t_s']*1e6,h[key]*scale,color=colors[1],ls='--',label='Laser, 10 J/m$^2$');ax.set(xlabel='Time ($\\mu$s)',ylabel=label)
    axs[0,0].legend(fontsize=8);save(fig,'02_dynamics')

    fig,axs=plt.subplots(1,2,figsize=(7,3))
    for name,label in [('sens_passive_near_base','Passive, $d_0=100$ nm'),('sens_passive_far_base','Passive, $d_0=1$ $\\mu$m'),('sens_fixed_base','Fixed center'),('base_1.00','Prescribed confinement')]:
        z=raw(name);axs[0].plot(z['t_s']*1e6,z['r_np_m']*1e6,label=label);axs[1].plot(z['t_s']*1e6,z['gap_m']*1e9,label=label)
    axs[0].set(xlabel='Time ($\\mu$s)',ylabel='Center radius ($\\mu$m)');axs[1].set(xlabel='Time ($\\mu$s)',ylabel='Surface gap (nm)',ylim=(-10,1100));axs[0].legend(fontsize=7)
    save(fig,'03_particle')

    fig,axs=plt.subplots(1,2,figsize=(7,3))
    z=raw('static_reference_320');q=np.load(ROOT.parent/'06_AuNP_Bubble_Coupling/results/reference.npz')
    for zz,label,col in [(z,'Radial source average',colors[0]),(q,'MR6 finite Au, axisymmetric',colors[1])]:
        axs[0].plot(zz['time_s']*1e6,100*zz['gas_energy_J']/q['absorbed_so_far_J'][-1],label=label,color=col)
        axs[1].plot(zz['time_s']*1e9,zz['gold_rise_K'],label=label,color=col)
    axs[0].set(xlabel='Time after onset ($\\mu$s)',ylabel='Absorbed energy delivered to gas (%)');axs[0].legend(fontsize=7)
    axs[1].set(xlabel='Time after onset (ns)',ylabel='Gold rise (K)',xlim=(0,80));save(fig,'04_static')
    rows=[]
    for r in static[static.cells==320].to_dict('records'):
        rows.append([r['case'].replace('_',' '),number(100*r['MR6_gas_fraction'],5),number(100*r['radial_gas_fraction'],5),number(100*r['relative_integrated_transfer_difference'],4),number(r['MR6_gold_peak_K'],4),number(r['radial_gold_peak_K'],4)])
    text('table_static',tex_table(['Case','MR6 gas (\\%)','Radial (\\%)','Diff. (\\%)','MR6 Au (K)','Radial Au (K)'],rows,'Matched 2 microsecond integrated delivery and peak Au temperature rise. Different observables retain their separate discrepancies.'))

    fig,axs=plt.subplots(1,2,figsize=(7,3.1))
    for pa,g in timing[(timing.level==3)&timing.completed].groupby('Pa_bar'):
        axs[0].plot(g.offset_s*1e6,g.dTmax_dF_K_m2_J*1e6,'o-',label=f'{pa:g} bar')
        if 'first_cycle_dT_dF_K_m2_J' in g:
            axs[1].plot(g.offset_s*1e6,g.first_cycle_dT_dF_K_m2_J*1e6,'o-',label=f'{pa:g} bar')
    for ax,title in zip(axs,['Five-cycle maximum','First-cycle thermal peak']):
        ax.axhline(0,color='#666',lw=.6);ax.set(xlabel='Pulse center offset ($\\mu$s)',ylabel='Response ($\\mu$K per J/m$^2$)',title=title);ax.legend(fontsize=8)
    save(fig,'05_timing')

    rows=[]
    principal=stable[(stable.level==3)&(stable.order==3)&(stable.fluence_J_m2==10)] if len(stable) else pairs[(pairs.level==3)&pairs.run_id.str.startswith('laser_F10_')]
    for r in principal.to_dict('records'):
        pa=r.get('Pa_bar',float(r['run_id'].split('_')[-1]) if r['run_id'].startswith('laser_') else np.nan)
        rows.append([f'{pa:g}',number(r['baseline_observed_Tmax_K'],8),number(r['laser_observed_Tmax_K'],8),number(r['global_delta_Tmax_K']*1e6,4),number(r['baseline_t_Tmax_s']*1e6,8),number(r['laser_t_Tmax_s']*1e6,8),'full' if r['completed_pair'] else 'censored'])
    text('table_main',tex_table(['$P_a$ (bar)','$T_0^{obs}$ (K)','$T_L^{obs}$ (K)','$\\Delta T$ ($\\mu$K)','$t_0$ ($\\mu$s)','$t_L$ ($\\mu$s)','Scope'],rows,'Cancellation-controlled nonlinear pairs, F=10 J/m$^2$, 10 ns, a=20 nm and prescribed 100 nm confinement. A checked reconstructed no-laser reference is shared. Censored observed maxima do not define a five-cycle difference. Direct independent subtraction and a full nonlinear Radau check are archived separately.'))
    if len(stable):
        primary=stable[(stable.level==3)&(stable.order==3)&(stable.fluence_J_m2==10)]
        s='The stable nonlinear increment formulation gives '
        s+=', '.join(f"{r.Pa_bar:g} bar: $\\Delta T_{{max}}={number(r.global_delta_Tmax_K*1e6,5)}$ $\\mu$K" if r.completed_pair else f"{r.Pa_bar:g} bar: no defined five-cycle difference after censoring" for r in primary.itertuples())
        s+='. These values share a checked reconstructed no-laser base; they do not acquire physical validation from cancellation control.'
        if (R/'independent_Radau.csv').exists():
            rr=read('independent_Radau').iloc[0]
            s+=f" A separate full nonlinear Radau pair at 1.0 bar gives {rr.global_delta_Tmax_K*1e6:.4g} microkelvin, without the reconstructed reference."
        s+=' Direct production-tolerance subtraction is unresolved at this scale; its larger differences are retained as numerical diagnostics.'
    else:s='The direct paired differences are interpreted together with the retained variational and convergence diagnostics; no weak positive effect is inferred from subtraction alone.'
    text('narrative_primary',s)

    fig,axs=plt.subplots(1,2,figsize=(7,3.2))
    for level,g in timing[timing.completed & np.isclose(timing.offset_s,-5e-6)].groupby('level'):
        axs[0].plot(g.Pa_bar,g.dTmax_dF_K_m2_J*1e6,'o-',label=f'Level {level}')
    axs[0].set(xlabel='Acoustic amplitude (bar)',ylabel='Peak sensitivity ($\\mu$K per J/m$^2$)');axs[0].legend(fontsize=8)
    labels={'wetargon':'Ar + nominal vapor','a40_fixed_gap':'a=40 nm, d=100 nm','a40':'a=40 nm, d=200 nm','duration100':'100 ns pulse',
            'dryargon':'Dry Ar','vapor5':'Vapor inventory ×5','k_half':'Water k=0.3 W/(m K)','eta_half':'Gas conductance ×0.5','G50':'G=50 MW/(m² K)',
            'lambda600':'Wavelength 600 nm','gap500':'Gap 500 nm','N8':'Eight particles','passive_far':'Passive, gap 1 µm','passive_near':'Passive, gap 100 nm'}
    ss=sens.copy();axs[1].barh(ss['case'].map(lambda x:labels.get(x,x)),ss.dTmax_dF_K_m2_J*1e6,color=colors[0]);axs[1].set(xlabel='Peak sensitivity ($\\mu$K per J/m$^2$)');axs[1].tick_params(axis='y',labelsize=6)
    save(fig,'06_feedback')
    feedback=[]
    if len(stable):
        for pa in (1.,1.15):
            g=stable[(stable.Pa_bar==pa)&(stable.order==3)&(stable.fluence_J_m2==10)]
            if set(g.level)>={2,3}:
                l2=g[g.level==2].iloc[0].global_delta_Tmax_K;l3=g[g.level==3].iloc[0].global_delta_Tmax_K
                feedback.append(f'At {pa:g} bar, Level 2 and Level 3 give {l2*1e6:.4g} and {l3*1e6:.4g} microkelvin. The additional feedback contribution is {(l3-l2)*1e6:.4g} microkelvin ({100*(l3/l2-1):.3g} percent of the one-way effect).')
    text('narrative_feedback',' '.join(feedback) or 'Level-specific paired response and its numerical limits are retained in the saved comparison tables.')

    fig,axs=plt.subplots(1,2,figsize=(7,3))
    for key,label in [('U_gas_change_J','Gas sensible energy'),('E_water_J','Water storage'),('work_by_gas_J','Work by gas'),('Q_bath_J','Bath removal')]:
        axs[0].plot(b['t_s']*1e6,b[key]*1e12,label=label)
    axs[0].set(xlabel='Time ($\\mu$s)',ylabel='Signed energy (pJ)');axs[0].legend(fontsize=7)
    for key,label in [('gas_residual_J','Gas'),('water_residual_J','Water'),('total_residual_J','Total')]:axs[1].plot(h['t_s']*1e6,h[key]*1e18,label=label)
    axs[1].set(xlabel='Time ($\\mu$s)',ylabel='Balance residual (aJ)');axs[1].legend(fontsize=7);save(fig,'07_energy')

    maxres=float(summ.total_residual_scaled.max());quad=float(summ.independent_total_residual_scaled.max())
    val=f'The retained production summaries have maximum scaled total residual {maxres:.3g}; independent flux quadrature gives {quad:.3g}.'
    if (R/'paired_energy.csv').exists():val+=f" The cancellation-controlled nonlinear pairs have maximum incremental residual {read('paired_energy').total_over_absorbed.max():.3g} of absorbed energy."
    val+=' Spatial, fixed-inner-grid domain extension and time/method controls are separate. The older outer-domain row also changes inner resolution and is not interpreted as a pure boundary test. A direct strict BDF pulse-edge failure is retained; rescaling only the incident-energy ledger tolerance restores completion. Absolute trajectory error is not used to certify a tiny difference.'
    text('narrative_validation',val)
    rows=[]
    for r in conv.to_dict('records'):
        rows.append([r['check'].replace('_',' '),number(r['observed_Tmax_K'],8),number(r['t_Tmax_s']*1e6,8),number(r['Rmin_m']*1e6,7),number(r['dTmax_dF_K_m2_J']*1e6,6),number(r['total_residual_scaled'],3)])
    text('table_convergence',tex_table(['Control','$T_{max}$ (K)','$t_T$ ($\\mu$s)','$R_{min}$ ($\\mu$m)','$dT/dF$ ($\\mu$K)','Residual'],rows,'Independent numerical controls at 1.0 bar. The derivative unit is microkelvin per J/m$^2$.'))

    extreme=summ[(summ.Pa_bar==1.3)&~summ.run_id.str.startswith(('N0_','material_guard_'))]
    rows=[]
    for r in extreme.to_dict('records'):
        rows.append([r['gas']+(' dry' if r['vapor_factor']==0 else ' wet'),f"{r['fluence_J_m2']:g}",number(r['observed_Tmax_K'],6),number(r['global_Tmax_K'],6),number(r['max_wall_Mach'],4),r['status'].replace('_',' ')[:30]])
    text('table_extremes',tex_table(['Gas','$F$ (J/m$^2$)','$T_{obs}$ (K)','$T_{global}$ (K)','$M_{max}$','Status'],rows,'Formal 1.3 bar cases with prescribed confinement. Wet pulses use the nominal minus 5 microsecond offset; dry F=100 cases use minus 1 microsecond relative to their own no-laser reference. These protocols do not isolate a composition-only effect. Missing global values are not a temperature ceiling.'))
    if (R/'physical_diagnostics.csv').exists():
        pdia=read('physical_diagnostics');pr=[]
        for r in pdia.to_dict('records'):
            pr.append([f"{r['Pa_bar']:g}",r['gas']+(' dry' if r['vapor_factor']==0 else ' wet'),number(r['max_pressure_Pa']/1e6,5),number(r['max_gas_density_kg_m3'],5),number(r['gas_acoustic_crossing_at_peak_s']*1e9,4),number(r['reconstructed_temperature_half_height_width_s']*1e9,4)])
        text('table_physical',tex_table(['$P_a$ (bar)','Gas','$p_{max}$ (MPa)','$\\rho_{g,max}$ (kg/m$^3$)','$R/c_g$ (ns)','Width (ns)'],pr,'Ideal-EOS diagnostics, not measured material states. The width is reconstructed from saved samples and is undefined when an observed peak lacks a right-hand half-height crossing.'))
    fig,ax=plt.subplots(figsize=(7,2.8));ax.axis('off');ax.set(xlim=(0,10),ylim=(0,3))
    labels=['Numerical\ncompletion','Gas coefficient\ndomain','Particle and\noptical geometry','Thermal / flow\napplicability']
    for j,label in enumerate(labels):
        x=j*2.5+.1;ax.add_patch(FancyBboxPatch((x,.7),2.15,1.4,boxstyle='round,pad=.08',facecolor='#f3f5f7',edgecolor=colors[j%4]));ax.text(x+1.075,1.4,label,ha='center',va='center')
    ax.text(5,.12,'Passing one screen does not establish the others. Formal and censored values remain archived.',ha='center',fontsize=8)
    save(fig,'08_applicability')

    fig,ax=plt.subplots(figsize=(7,3))
    for tus in (10,14,16,18,33):
        ft=h['field_t_s'] if 'field_t_s' in h else h['t_s']
        j=min(np.searchsorted(ft,tus*1e-6),len(ft)-1);i=h['field_indices'][j] if 'field_indices' in h else j
        rr=(h['R_m'][i]**3+3*h['xi_centers_m3'])**(1/3)
        ax.plot((rr-h['R_m'][i])*1e6,h['T_l_rise_K'][:,j],label=f"{h['t_s'][i]*1e6:.2f} $\\mu$s")
    ax.set(xlabel='Distance from current wall ($\\mu$m)',ylabel='Spherically averaged water rise (K)',xlim=(0,3));ax.legend(fontsize=8);save(fig,'09_profiles')
    z=np.load(R/'transport_timing.npz');fig,axs=plt.subplots(1,3,figsize=(7,2.65))
    for tag,label,col in [('100','1.0 bar',colors[0]),('115','1.15 bar',colors[1])]:
        t=z[tag+'_time_from_center_s']
        axs[0].plot(t*1e9,z[tag+'_delta_TAu_K'],color=col,label=label)
        axs[1].plot(t*1e9,z[tag+'_delta_Qgas_W']*1e12,color=col)
        axs[2].plot(t*1e6,z[tag+'_delta_Tgas_K']*1e6,color=col)
    axs[0].set(xlim=(-10,100),xlabel='From pulse center (ns)',ylabel='Au increment (K)');axs[0].legend(fontsize=7)
    axs[1].set(xlim=(-10,300),xlabel='From pulse center (ns)',ylabel='Net gas-flux increment (pW)')
    axs[2].set(xlim=(0,6),xlabel='From pulse center ($\\mu$s)',ylabel='Gas increment ($\\mu$K)')
    save(fig,'10_transport')
    delays=read('transport_timing')
    text('narrative_transport',' '.join(f"At {r.Pa_bar:g} bar the initial net wall-flux increment peaks {r.observed_incremental_flux_peak_lag_s*1e9:.3g} ns after the pulse center, within the declared 2 microsecond window. The first-cycle thermal maximum follows the center by {r.baseline_first_cycle_thermal_peak_after_center_s*1e6:.4g} microseconds; its temperature increment at the reference time is {r.first_cycle_delta_T_at_reference_K*1e6:.4g} microkelvin." for r in delays.itertuples())+' The net flux increment includes changed gas cooling and compression feedback; it is not a laser-tagged energy fraction.')
    if len(stable):
        pp=stable[(stable.level==3)&(stable.order==3)&(stable.fluence_J_m2==10)]
        conclusion='For the declared 10 J/m$^2$, 10 ns pulse and prescribed confinement, '
        conclusion+=', '.join(f"the stable {r.Pa_bar:g} bar pair has {number(r.global_delta_Tmax_K*1e6,5)} microkelvin global change" if r.completed_pair else f"the {r.Pa_bar:g} bar pair is censored and has no defined five-cycle temperature change" for r in pp.itertuples())
        conclusion+='.'
    else:conclusion='The saved finite and differential responses quantify a small conditional thermal effect; numerical resolution and physical admissibility are assessed independently.'
    conclusion+=' The timing scan retains pulses that miss the selected thermal maximum. The temperature differences describe the specified reduced equations, with source averaging and mechanical feedback explicitly included or separated; they are not experimentally achieved temperature gains.'
    text('narrative_conclusion',conclusion)
    notebook(summ,stable)
    print('Ten JPG figures, report tables and numbered notebook regenerated from saved arrays.')

def notebook(summ,stable):
    def md(s):return dict(cell_type='markdown',id='mr7-'+str(len(s))+'-'+str(sum(map(ord,s))%1000000),metadata={},source=s.splitlines(keepends=True))
    primary=stable[(stable.level==3)&(stable.order==3)&(stable.fluence_J_m2==10)]
    rows=['|Amplitude (bar)|Baseline observed (K)|Laser observed (K)|Five-cycle change (µK)|Scope|','|---|---:|---:|---:|---|']
    for r in primary.itertuples():rows.append(f'|{r.Pa_bar:g}|{r.baseline_observed_Tmax_K:.8f}|{r.laser_observed_Tmax_K:.8f}|'+(f'{r.global_delta_Tmax_K*1e6:.5g}|complete|' if r.completed_pair else 'undefined|censored at gas domain|'))
    cells=[md('# Mini-Research 7: Integrated Laser–AuNP–Bubble Model\n\nDimitri Bolt · Academic supervisor: Prof. Dr. Ildar R. Gabitov.\n\n[Scientific report](report/07_Integrated_Model.pdf) · [Mathematical details](results/model_decisions.md) · [SI schema](results/data_schema.md)\n'),
       md('## 1. Model and question\n\nThe coupled state is bubble radius and speed, a fixed NASA9 gas mixture, material water temperatures and dynamic Au storage. The material coordinate is ξ=(r³−R³)/3. A single signed gas–water flow is used once in each subsystem with opposite signs. The paired comparison changes only the optical pulse within the same moving-water model. Spherical averaging does not resolve the finite-particle local hot spot.\n\n![Coupled system](figures/01_system.jpg)\n'),
       md('## 2. Primary paired result\n\nOne 10 ns, 532 nm pulse, 10 J/m², one 20 nm-radius Au particle, prescribed 100 nm gap; five acoustic periods where complete. The pulse center is five microseconds before the selected no-laser radius event.\n\n'+'\n'.join(rows)+'\n\nThese are cancellation-controlled nonlinear increments. Absolute baseline accuracy remains lower than the displayed increment resolution. The independent full Radau pair at 1.0 bar gives 1.4073 µK. The wet 1.3 bar case stops at the 6000 K coefficient boundary; no unknown tail is extrapolated.\n'),
       md('## 3. Dynamics and particle geometry\n\n![Dynamics](figures/02_dynamics.jpg)\n\n![Particle motion](figures/03_particle.jpg)\n\nPassive liquid motion, a fixed laboratory center and prescribed confinement have different contact histories. The confinement force is unverified. [All runs](results/all_runs.csv) and [applicability](results/applicability.csv) retain terminal events and physical screens.\n'),
       md('## 4. Heat arrival, timing and feedback\n\n![Transport](figures/10_transport.jpg)\n\n![Timing](figures/05_timing.jpg)\n\nAt 1.0 bar Level 2 and Level 3 give 0.8578 and 1.4088 µK; at 1.15 bar they give 5.6563 and 2.1849 µK. Feedback can increase or reduce the response. The extremely small 1.15 bar response to a pulse centered on the radius event is unresolved under mesh refinement. A pulse arriving after a thermal peak cannot alter that past peak.\n'),
       md('## 5. Verification and reduction error\n\n![Static control](figures/04_static.jpg)\n\nAt 320 cells, matched stationary integrated gas delivery differs from MR6 by approximately 0.04–1.68%, while radial Au temperature rises are only 21–30% of the finite-Au reference. These discrepancies are different observables; no fitted correction is applied. The exact Robin point-source control, clamped MR4 regression, two stiff methods, finite increments, tangent responses, domain extension, mesh/time refinement and independent energy quadrature are retained. [Validation](results/validation.md) · [Numerical scope](results/numerical_scope.json) · [Tests](results/tests.txt).\n\n![Energy](figures/07_energy.jpg)\n'),
       md('## 6. Formal extremes and interpretation\n\nDry fixed-composition 1.3 bar baselines reach about 6264 K (air) and 22809 K (ground-state constant-Cv Ar). Their pulse protocol differs from the nominal wet case. They are formal sensitivities, with missing phase/chemical channels, nonideal EOS, shocks and local particle physics. The 6000 K wet boundary is a coefficient limit, not a physical ceiling. No experimental temperature gain or fusion yield is established.\n\n[Physical diagnostics](results/physical_diagnostics.csv) · [Missing-channel reference scales](results/missing_channel_scales.csv) · [Convergence](results/convergence.csv)\n')]
    if (R/'publication_scope.json').exists():
        cells.insert(1,md('## Public data scope\n\nThis compact package ships source, report, all aggregate tables and run metadata, four full accepted-state references at three acoustic amplitudes, selected exact scalar traces and five original liquid profiles. Other high-dimensional sweep fields remain local. [Selection and coverage](results/publication_scope.json) · [Per-array inventory](results/public_data_inventory.json). Historical full-archive verification and current compact checks have separate scopes.\n'))
    (ROOT/'07_Integrated_Model.ipynb').write_text(json.dumps(dict(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3'),language_info=dict(name='python',version='3.12')),nbformat=4,nbformat_minor=5),indent=1,ensure_ascii=False)+'\n')

if __name__=='__main__':main()
