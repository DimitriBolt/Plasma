"""All figures/tables read saved numerical data; JPG only, English labels."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Circle,Rectangle,FancyArrowPatch
from .upstream import ROOT

R=ROOT/'results'; F=ROOT/'figures'
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,
                     'figure.dpi':120,'savefig.dpi':300,'axes.spines.top':False,'axes.spines.right':False,
                     'font.family':'DejaVu Sans'})
def save(fig,name):
    fig.savefig(F/(name+'.jpg'),dpi=300,bbox_inches='tight',facecolor='white',pil_kwargs={'quality':95});plt.close(fig)
def read(name):return pd.read_csv(R/(name+'.csv'),float_precision='round_trip')

def main():
    ref=np.load(R/'reference.npz');t=ref['time_s'];tc=0.;d=100e-9;a=20e-9;rb=6e-6
    fig,ax=plt.subplots(figsize=(7,2.55));ax.set_aspect('equal');ax.set_xlim(-.9,1.5);ax.set_ylim(-.46,.64)
    ax.add_patch(Circle((0,-1.2),1.2,color='#d8ebfa',ec='#2171b5'))
    ax.add_patch(Circle((0,.24),.075,color='#e6b333',ec='#8b6400'))
    ax.annotate('AuNP: lumped $T_{Au}$',(0,.29),(.35,.48),arrowprops={'arrowstyle':'->'})
    ax.annotate('Bubble: $R$, uniform reference $T_g$',(.45,-.15),(.6,-.34),arrowprops={'arrowstyle':'->'})
    ax.annotate('Surface gap $d=r_{np}-R-a$',(0,.1),(.35,.13),arrowprops={'arrowstyle':'->'})
    ax.annotate('Local laser $I_{local}\sigma_{abs}$',(-.07,.27),(-.8,.49),arrowprops={'arrowstyle':'->','color':'#c44e52'})
    ax.annotate('Resolved liquid diffusion',(-.13,.13),(-.85,-.05),arrowprops={'arrowstyle':'->'})
    ax.arrow(-.22,-.03,0,.22,width=.008,head_width=.06,color='#2171b5',length_includes_head=True)
    ax.text(-.62,.22,'$\\mathbf{n}_g$ into liquid');ax.axis('off');ax.text(1.45,.57,'Schematic; not to scale',ha='right',fontsize=8)
    save(fig,'01_geometry')
    q=ref['q_wall_laser_increment_W'];idx=int(np.argmax(q));xy=ref['xy_m'];tri=ref['triangles']
    field=ref['water_rise_K'][:,idx];fig,(ax,bx)=plt.subplots(1,2,figsize=(7,2.65),gridspec_kw={'width_ratios':[1.1,1]})
    triang=mtri.Triangulation(xy[:,0]*1e9,xy[:,1]*1e9,tri)
    image=ax.tripcolor(triang,field,shading='gouraud',cmap='inferno',vmin=0,vmax=max(field))
    ax.set_xlim(0,400);ax.set_ylim(-30,400);ax.set_aspect('equal');ax.set_xlabel('Cylindrical radius (nm)');ax.set_ylabel('Height above apex (nm)')
    ax.set_title(f'Water rise at flux peak, {t[idx]*1e9:.1f} ns');fig.colorbar(image,ax=ax,label='Laser temperature rise (K)',pad=.02)
    wall=ref['wall_node_area_m2']>0;xx=xy[wall,0];zz=xy[wall,1];theta=np.arctan2(xx,zz+rb)
    flux=field[wall]*float(ref['h_W_m2K']);order=np.argsort(theta);mean=q[idx]/(4*np.pi*rb**2)
    plotted=order[flux[order]>0]
    bx.semilogy(theta[plotted]*180/np.pi,flux[plotted],color='#2171b5',label='Resolved positive local flux')
    bx.axhline(mean,color='#c44e52',ls='--',label='Whole-surface mean')
    bx.set_xlabel('Polar angle from particle (degree)');bx.set_ylabel('Laser increment of flux (W/m$^2$)');bx.set_xlim(0,40)
    bx.set_ylim(max(flux)*1e-5,max(flux)*2);bx.legend();fig.tight_layout();save(fig,'02_local_temperature_flux')
    # Preserve the signed field and record local concentration metrics.
    weights=ref['wall_node_area_m2'][wall];positive=np.maximum(flux,0);sort=np.argsort(positive)[::-1]
    cumulative=np.cumsum(positive[sort]*weights[sort]);j=np.searchsorted(cumulative,.9*cumulative[-1]);area90=weights[sort[:j+1]].sum()
    local=dict(time_s=float(t[idx]),max_local_flux_W_m2=float(max(flux)),surface_mean_flux_W_m2=float(mean),
               max_over_mean=float(max(flux)/mean),area90_m2=float(area90),area90_fraction=float(area90/(4*np.pi*rb**2)),
               signed_integral_W=float(np.sum(flux*weights)),negative_tail_flux_integral_W=float(np.sum(np.minimum(flux,0)*weights)))
    (R/'local_flux_metrics.json').write_text(json.dumps(local,indent=2)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(7,2.2))
    axes[0].plot(t*1e9,ref['gold_rise_K'],label='Au');axes[0].plot(t*1e9,ref['nearest_water_rise_K'],label='Facing Au surface water')
    axes[0].set_ylabel('Laser temperature rise (K)');axes[0].legend();axes[0].set_xlim(0,100)
    axes[1].plot(t*1e9,q*1e12);axes[1].set_ylabel('Liquid-to-gas laser power (pW)');axes[1].set_xlim(0,250)
    axes[2].plot(t*1e9,ref['gas_energy_J']*1e18);axes[2].set_ylabel('Energy reaching gas (aJ)');axes[2].set_xlim(0,2000)
    for ax in axes:ax.set_xlabel('Time from laser start (ns)');ax.grid(alpha=.2)
    fig.tight_layout();save(fig,'03_time_response')
    gap=read('gap_scan');sens=read('transport_sensitivity');fig,(ax,bx)=plt.subplots(1,2,figsize=(7,2.5))
    ax.loglog(gap.d_m*1e9,gap.gas_fraction*100,'o-');ax.axvspan(.8,2,color='#c44e52',alpha=.15)
    ax.axvline(100,ls=':',color='gray');ax.set_xlabel('Surface gap $d$ (nm)');ax.set_ylabel(r'Gas delivery by 2 $\mu$s (%)')
    ax.set_title('20 nm Au, 10 ns, $h=4400$ W/(m$^2$ K)')
    hh=sens[(sens.varied_parameter=='h')&(sens.value_SI>0)];bx.semilogx(hh.value_SI,hh.gas_fraction*100,'o-',color='#b57128')
    bx.set_xlabel('Effective gas transfer $h$ (W/m$^2$/K)');bx.set_ylabel('Transfer by 2 $\mu$s (%)');bx.set_title('Gap = 100 nm; high $h$ is a mathematical limit')
    for ax in [ax,bx]:ax.grid(alpha=.2)
    fig.tight_layout();save(fig,'04_gap_transport')
    timing=read('timing_scan');matrix=read('gap_timing_map');pivot=matrix.pivot(index='d_m',columns='delta_t_s',values='late_window_fraction')
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7,2.5));ax.plot(timing.delta_t_s*1e9,timing.efficiency_by_tc*100,'o-',label='Delivered before $t_c$')
    ax.plot(timing.delta_t_s*1e9,timing.efficiency_late_window*100,'s-',label='$t_c \pm 10$ ns window');ax.set_xlim(-1000,150)
    ax.set_xlabel('$\Delta t=t_{laser,center}-t_c$ (ns)');ax.set_ylabel('Fraction of local absorption (%)');ax.legend()
    image=bx.imshow(pivot.to_numpy()*100,aspect='auto',origin='lower',cmap='viridis')
    bx.set_xticks(np.arange(len(pivot.columns))[::2],[f'{x*1e9:.0f}' for x in pivot.columns[::2]])
    bx.set_yticks(np.arange(len(pivot.index)),[f'{x*1e9:.0f}' for x in pivot.index]);bx.set_xlabel('Laser center relative to $t_c$ (ns)');bx.set_ylabel('Surface gap (nm)')
    fig.colorbar(image,ax=bx,label='Late-window fraction (%)',pad=.02);fig.tight_layout();save(fig,'05_timing')
    kin=read('kinematics');fig,axes=plt.subplots(1,2,figsize=(7,2.5))
    for ax,amp in zip(axes,[1.15,1.3]):
        for mode,label in [('fixed','Fixed center'),('passive','Passive incompressible flow'),('capture','Prescribed capture')]:
            z=kin[(kin.Pa_bar==amp)&(kin['mode']==mode)]
            ax.plot(z.t_minus_tc_s*1e9,z.surface_gap_m*1e9,label=label)
        ax.axhline(0,color='black',lw=.6);ax.set_xlim(-400,100);ax.set_yscale('symlog',linthresh=100)
        ax.set_xlabel('Time relative to radius minimum (ns)');ax.set_ylabel('Surface gap (nm)');ax.set_title(f'{amp:g} bar: kinematic diagnostic')
    axes[0].legend();fig.tight_layout();save(fig,'06_kinematics')
    gas=read('gas_diagnostic');amp=1.15;tab=gas[gas.Pa_bar==amp];base=tab[tab.mechanism=='baseline'].iloc[0]
    with np.load(R/'histories'/f'gas_{amp:.2f}_baseline.npz') as z:bt=z['time_s'].copy();bT=z['T_K'].copy()
    from scipy.interpolate import PchipInterpolator
    bf=PchipInterpolator(bt,bT,extrapolate=False)
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7,2.5))
    for mechanism,label in [('A_preheat','A: early heating'),('B_wall','B: diffusive wall'),('B_after','B: after collapse'),('C_B_energy','C: same B energy'),('C_ideal_upper','C: ideal 100%')]:
        with np.load(R/'histories'/f'gas_{amp:.2f}_{mechanism}.npz') as z:
            x=z['time_s'];y=z['delta_T_K'];xx=z['time_minus_tc_s']*1e9
        ax.plot(xx,y,label=label)
    ax.set_xlim(-150,250);ax.set_yscale('symlog',linthresh=1e-7);ax.set_ylim(0.,.007)
    ax.set_xlabel('Time relative to $t_c$ (ns)');ax.set_ylabel('Gas temperature increment (K)');ax.legend(fontsize=6.5)
    labels=['A','B','C: B energy','C: 100%'];mechs=['A_preheat','B_wall','C_B_energy','C_ideal_upper']
    vals=[tab[tab.mechanism==mech].direct_source_energy_J.iloc[0]*1e15 for mech in mechs]
    bx.bar(np.arange(4),vals,color=['#70a896','#5089ad','#d2aa62','#c77373']);bx.set_xticks(np.arange(4),labels,rotation=15)
    bx.set_yscale('log');bx.set_ylabel('Source energy (fJ)');bx.set_title('1.15 bar\nFormal frozen diagnostics')
    fig.tight_layout();save(fig,'07_ABC_diagnostic')
    fig,(ax,bx)=plt.subplots(1,2,figsize=(7,2.5))
    E=ref['absorbed_so_far_J'][-1];select=[np.argmin(abs(t-v)) for v in [10e-9,100e-9,2e-6]]
    bottom=np.zeros(3)
    for key,label,color in [('gold_energy_J','Au storage','#d7ad43'),('water_energy_J','Water storage','#69a3c3'),('gas_energy_J','Gas received','#5a956e'),('bath_energy_J','Outer bath','#ad7aaa')]:
        value=ref[key][select]/E*100;ax.bar(np.arange(3),value,bottom=bottom,label=label,color=color);bottom+=value
    ax.set_xticks(np.arange(3),['10 ns','100 ns','2 $\mu$s']);ax.set_ylabel('Fraction of absorbed energy (%)');ax.legend(fontsize=7,loc='lower right')
    bx.semilogx(t[1:]*1e9,ref['gas_energy_J'][1:]/E*100,color='#5a956e')
    bx.set_xlabel('Time from laser start (ns)');bx.set_ylabel('Gas received (%)');bx.set_ylim(bottom=0)
    fig.tight_layout();save(fig,'08_energy_partition')
    # English numbered notebook, saved-result navigation only.
    cells=[]
    def md(s):cells.append(dict(cell_type='markdown',metadata={},source=[s]))
    def code(s):cells.append(dict(cell_type='code',metadata={},source=[s],execution_count=None,outputs=[]))
    md('# Mini-Research 6: Frozen-Geometry AuNP-Bubble Thermal Coupling\nDimitri Bolt. Academic supervisor: Prof. Dr. Ildar R. Gabitov.\n\nThe completed scope is stationary thermal coupling. Motion/advection and self-consistent mechanics belong to MR7. This notebook reads saved results. [Scientific report](report/06_AuNP_Bubble_Coupling.pdf); [model details](results/model_details.md).')
    md('## 1. Model and energy channel\nThe axisymmetric liquid equation resolves a finite Au sphere and a spherical gas boundary in frozen geometry. The signed MR4 gas loss is supplemented only by the laser-tagged increment. All dynamic event applications are conditional diagnostics.\n\n'
       f'Reference: a=20 nm, surface gap=100 nm, R=6 micrometres, 532 nm, local F=1 J/m², 10 ns. Absorbed energy: **{E*1e15:.5f} fJ**. Gas energy by 2 microseconds: **{ref["gas_energy_J"][-1]*1e18:.4f} aJ**, or **{ref["gas_energy_J"][-1]/E*100:.4f}%**.\n\n![Reference transient](figures/03_time_response.jpg)')
    code("from pathlib import Path\nimport json, numpy as np, pandas as pd\nfrom IPython.display import display, Image\nr=Path('results')\ndisplay(pd.read_csv(r/'gap_scan.csv'))")
    md('## 2. Spatial localization and timing\nThe whole-surface mean does not represent the local hot patch. Timing shifts below use the calculated diffusion response.\n\n![Local field and wall flux](figures/02_local_temperature_flux.jpg)\n\n![Gap and timing](figures/05_timing.jpg)')
    code("display(Image(filename='figures/02_local_temperature_flux.jpg'))\ndisplay(Image(filename='figures/05_timing.jpg'))")
    md('## 3. Particle kinematics and geometry admission\nThe passive invariant is r_np³-R³=constant; the exterior surface gap increases during collapse. Fixed-center and prescribed capture paths are distinct. All late strong-event frozen snapshots fail the flow/displacement screens.\n\n![Kinematic scenarios](figures/06_kinematics.jpg)')
    code("display(pd.read_csv(r/'mechanical_events.csv'))\ndisplay(pd.read_csv(r/'event_snapshots.csv'))\ndisplay(Image(filename='figures/06_kinematics.jpg'))")
    md('## 4. Finite design objectives and formal extremes\nThe 216-point design has separate absorbed-energy, gas-energy and fraction objectives. The moderate static maximum gas delivery is about 3.95 fJ by 2 microseconds from about 3.48 pJ absorption. Wavelength plateaus and unresolved small gap differences do not identify a physical optimum. Formal short-pulse/near-contact states remain with their failure flags.\n\n![Gap and transfer coefficient](figures/04_gap_transport.jpg)')
    code("display(pd.read_csv(r/'finite_design.csv'))\ndisplay(pd.read_csv(r/'formal_regimes.csv'))\nprint((r/'design_extrema.json').read_text())")
    md('## 5. One-way gas diagnostic and coefficient censoring\nA is a **transient early-heating diagnostic**, without an equilibrium update of radius, vapor pressure or composition. A common nonlinear gas baseline and a directly integrated first-order laser increment avoid subtractive numerical noise. At 1.15 bar the B increment at the baseline peak is about 0.468 microkelvin; ideal 100% C delivery gives about 3.506 millikelvin. These are formal frozen diagnostics with substantial baseline pressure mismatch. The 1.3 bar histories stop at the actual NASA coefficient boundary; 6000 K is not a physical ceiling.\n\n![A/B/C diagnostic](figures/07_ABC_diagnostic.jpg)')
    code("display(pd.read_csv(r/'gas_diagnostic.csv'))\ndisplay(Image(filename='figures/07_ABC_diagnostic.jpg'))")
    md('## 6. Numerical checks and interpretation\nProduction-to-finest gas-energy discrepancy is at most 0.196% in the selected controls. The finest isolated Au peaks agree with MR5 within 0.0893%. Numerical convergence does not validate capture, interfacial optics, fixed transport, or a frozen strong collapse. The history interface cannot be applied to changed mechanics without evolving the liquid state.\n\n![Separate energy channels](figures/08_energy_partition.jpg)')
    code("print((r/'validation_summary.json').read_text())\ndisplay(pd.read_csv(r/'isolated_MR5.csv'))\ndisplay(pd.read_csv(r/'analytic_point_controls.csv'))")
    nb=dict(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.12'}},nbformat=4,nbformat_minor=5)
    (ROOT/'06_AuNP_Bubble_Coupling.ipynb').write_text(json.dumps(nb,indent=1)+'\n')
    # Numerical macros are derived exclusively from saved outputs.
    v=json.loads((R/'validation_summary.json').read_text());ext=json.loads((R/'design_extrema.json').read_text());tim=read('timing_scan')
    best=tim.loc[tim.efficiency_late_window.idxmax()];gasmax=ext['moderate_max_gas_energy'];bound=ext['moderate_max_absorbed']
    values=dict(RefEabsFJ=f'{E*1e15:.5f}',RefGasAJ=f'{ref["gas_energy_J"][-1]*1e18:.4f}',
                RefFraction=f'{ref["gas_energy_J"][-1]/E*100:.4f}',RefGold=f'{max(ref["gold_rise_K"]):.4f}',
                RefWater=f'{max(ref["nearest_water_rise_K"]):.4f}',RefFluxTime=f'{t[idx]*1e9:.2f}',
                RefFluxPW=f'{max(q)*1e12:.4f}',LocalRatio=f'{local["max_over_mean"]:.0f}',
                AreaPercent=f'{local["area90_fraction"]*100:.4f}',MeshGasPercent=f'{v["production_gas_energy_spatial_relative_max"]*100:.3f}',
                MeshGoldPercent=f'{v["production_gold_spatial_relative_max"]*100:.3f}',IsolatedPercent=f'{v["isolated_gold_peak_finest_relative_max"]*100:.3f}',
                ModerateGasFJ=f'{gasmax["gas_energy_J"]*1e15:.4f}',ModerateEabsPJ=f'{bound["energy_absorbed_J"]*1e12:.4f}',
                BestOffsetNS=f'{best.delta_t_s*1e9:.0f}',BestLatePercent=f'{best.efficiency_late_window*100:.5f}')
    (R/'report_values.tex').write_text(''.join('\\newcommand{\\'+key+'}{'+val+'}\n' for key,val in values.items()))
    print('eight JPG figures, notebook and report values written')

if __name__=='__main__': main()
