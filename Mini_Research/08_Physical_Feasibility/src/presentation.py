"""Saved-data tables and print-quality JPG figures. Does not integrate models."""
import json,shutil,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from .upstream import ROOT,REPO

R=ROOT/'results';F=ROOT/'figures';FINAL=REPO/'Final_Report'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})

def save(fig,name):
    fig.savefig(F/(name+'.jpg'),dpi=300,bbox_inches='tight',pil_kwargs={'quality':95,'subsampling':0})
    plt.close(fig)
    shutil.copy2(F/(name+'.jpg'),FINAL/'figures'/(name+'.jpg'))

def table(filename,headers,rows,widths=None):
    col='l'+'r'*(len(headers)-1) if widths is None else widths
    s='\\begin{tabular}{'+col+'}\\toprule\n'+' & '.join(headers)+' \\\\ \\midrule\n'
    s+='\n'.join(' & '.join(str(x) for x in row)+' \\\\' for row in rows)+'\n\\bottomrule\\end{tabular}\n'
    (R/filename).write_text(s);shutil.copy2(R/filename,FINAL/'tables'/filename)

def main():
    F.mkdir(exist_ok=True)
    for d in ['figures','tables']:(FINAL/d).mkdir(exist_ok=True)
    m=pd.read_csv(R/'mechanics.csv');s=m[m.control=='strict']
    fig,ax=plt.subplots(figsize=(10,4.3));ax.set(xlim=(0,10),ylim=(0,4.5));ax.axis('off')
    def box(x,y,w,h,text,color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.07',fc=color,ec='#273744',lw=1))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=10)
    def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,lw=1.5,color='#273744'))
    box(.15,3.05,2.4,.85,'Acoustic pressure\nR(t), U(t), p(t)','#e9edf0')
    box(3.3,3.05,2.4,.85,'Compression work\n−p dV/dt','#e9edf0');arrow((2.65,3.5),(3.2,3.5))
    box(.15,.8,1.8,.85,'Laser pulse\nI(t), F, λ','#fff1d6')
    box(2.35,.8,1.8,.85,'Au absorption\nand storage','#fff1d6');arrow((2,.1+1.1),(2.25,1.2))
    box(4.55,.8,1.8,.85,'Moving liquid\nheat transport','#fff1d6');arrow((4.2,1.2),(4.45,1.2))
    box(7.2,1.6,2.4,1.1,'Gas energy balance\nCv(T) dT/dt\n= −p dV/dt + Qℓ→g','#deebf5')
    arrow((5.8,3.45),(7.4,2.8));arrow((6.45,1.2),(7.4,1.5))
    ax.text(7.8,3.45,'Pressure feedback\nchanges mechanics',ha='center',fontsize=9)
    arrow((8.3,2.8),(5.9,3.85))
    ax.text(4.9,.18,'Finite particle placement, signed flux, composition and EOS limit interpretation.',ha='center',fontsize=9)
    save(fig,'01_framework')
    fig,axs=plt.subplots(2,3,figsize=(10,5.5))
    for j,a in enumerate([1.,1.15,1.3]):
        for model,color,style in [('KM','#1f4e79','-'),('Gilmore','#b35115','--')]:
            z=np.load(R/f'mechanics_{model}_{a:.2f}.npz');r=s[(s.Pa_bar==a)&(s.model==model)].iloc[0]
            axs[0,j].plot(z['t_s']*1e6,z['R_m']*1e6,style,c=color,label=model)
            mask=abs(z['t_s']-r.t_Rmin_s)<(1e-6 if a==1 else 80e-9)
            axs[1,j].plot((z['t_s'][mask]-r.t_Rmin_s)*1e9,z['R_m'][mask]*1e6,style,c=color)
        axs[0,j].set(title=f'{a:.2f} bar',xlabel='Absolute time (µs)',ylabel='Radius (µm)')
        axs[1,j].set(xlabel='Time from own minimum (ns)',ylabel='Radius (µm)')
        for ax in axs[:,j]:ax.grid(alpha=.2)
    axs[0,0].legend();fig.tight_layout();save(fig,'02_mechanics')
    scales=json.loads((R/'scales.json').read_text())
    names=['Au internal diffusion','Electron–phonon reference','Au/water interface','Laser pulse','Dry air hot width, formal','Wet 1.15 bar hot width','100 nm water diffusion','120 nm center distance diffusion','500 nm water diffusion']
    times=[scales['tau_Au_internal_s'],1.7e-12,scales['tau_Au_interface_s'],1e-8,3.231635e-9,17.861676e-9,scales['tau_diffusion_gap100nm_s'],scales['tau_diffusion_center120nm_s'],scales['tau_diffusion_gap500nm_s']]
    fig,ax=plt.subplots(figsize=(9,4.4));ax.barh(names,np.array(times)*1e9,color=['#516779']*4+['#a36b36']*2+['#447784']*3)
    ax.set_xscale('log');ax.set(xlabel='Characteristic time (ns); widths are inherited MR7 values',xlim=(.0005,10000));ax.invert_yaxis();ax.grid(axis='x',alpha=.25);fig.tight_layout();save(fig,'03_timescales')
    ideal=pd.read_csv(R/'ideal_caloric_allowance.csv');i=ideal[(ideal.gas=='air')&(ideal.vapor_factor==1)&(ideal.fluence_J_m2==10)].sort_values('baseline_T_K')
    fig,ax=plt.subplots(figsize=(7.6,3.7));x=np.arange(2);width=.32
    ax.bar(x-width/2,[1.4088162743e-6,2.1848791221e-6],width,label='MR7 coupled global peak change',color='#1f4e79')
    ax.bar(x+width/2,i.ideal_all_absorbed_at_fixed_volume_delta_T_K,width,label='MR8 ideal fixed-volume caloric allowance',color='#b35115')
    ax.set_yscale('log');ax.set_xticks(x,['1.00 bar: conditional','1.15 bar: formal']);ax.set_ylabel('Gas temperature increment (K)');ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2);fig.tight_layout();save(fig,'04_caloric_comparison')
    rows=[]
    for a in [1.,1.15,1.3]:
        for model in ['KM','Gilmore']:
            q=s[(s.Pa_bar==a)&(s.model==model)].iloc[0]
            rows.append([f'{a:.2f}',model,f'{q.Rmin_m*1e6:.6f}',f'{q.t_Rmin_s*1e6:.6f}',f'{q.inward_speed_m_s:.3f}',f'{q.max_M_ambient:.4f}',f'{q.chi:.4f}'])
    table('table_mechanics.tex',[r'$P_a$ (bar)','Model',r'$R_{\min}$ ($\mu$m)',r'$t_R$ ($\mu$s)',r'$|U_-|$ (m/s)',r'$M_0$',r'$\chi$'],rows)
    ea=pd.read_csv(R/'energy_audit.csv')
    rows=[]
    for name,label in [('stable_base_1.00_F10_L3_q3','1.00'),('stable_base_1.15_F10_L3_q3','1.15'),('stable_extra_high_1.30_BDF_base_F10_L3_q3','1.30*')]:
        q=ea[ea.run_id==name].iloc[0]
        rows.append([label]+[f'{q[k]*1e12:.6g}' for k in ['E_absorbed_J','E_Au_J','E_liquid_J','delta_U_gas_J','W_by_gas_J','Q_l_to_g_J']])
    table('table_energy.tex',[r'$P_a$ (bar)',r'$E_{\rm abs}$',r'$E_{\rm Au}$',r'$E_\ell$',r'$\Delta U_g$',r'$W_{\rm by}$',r'$Q_{\ell g}$'],rows)
    rows=[[f'{r.baseline_T_K:.2f}',r.gas+(' wet' if r.vapor_factor else ' dry'),f'{r.ideal_all_absorbed_at_fixed_volume_delta_T_K:.6g}'] for _,r in ideal[ideal.fluence_J_m2==10].iterrows()]
    table('table_ideal.tex',['Reference $T$ (K)','Frozen gas',r'Ideal $\Delta T$ (K)'],rows)
    rows=[]
    d=pd.read_csv(R/'particle_kinematics.csv')
    for _,r in d.iterrows():rows.append([r.run_id.split('_')[-1],f'{r.max_required_slip_m_s:.5g}',f'{r.max_particle_Re:.5g}'])
    table('table_retention.tex',[r'$P_a$ (bar)','Required slip (m/s)',r'Particle Reynolds number'],rows)
    c=pd.read_csv(R/'convergence.csv');p=pd.read_csv(R/'paired_energy_audit.csv')
    values=dict(max_Rmin_discrepancy=float(c.relative_Rmin.max()),max_speed_discrepancy=float(c.relative_speed.max()),
        max_absolute_event_difference_s=float(abs(c.event_time_difference_s).max()),
        max_thermal_scaled_residual=float(ea.total_scaled_residual.max()),
        max_paired_gas_residual_over_absorption=float(p.independent_gas_increment_residual_over_absorption.max()),
        max_stored_paired_residual_over_absorption=float(p.stored_paired_total_residual_over_absorption.max()))
    (R/'verification_values.json').write_text(json.dumps(values,indent=2)+'\n')
    for key in ['max_Rmin_discrepancy','max_speed_discrepancy','max_absolute_event_difference_s','max_thermal_scaled_residual','max_paired_gas_residual_over_absorption']:
        pass
    srcs=[]
    for stage,filename,target in [('05_Plasmonic_Heating','02_optical_spectra.jpg','05_optical_spectra.jpg'),('06_AuNP_Bubble_Coupling','02_local_temperature_flux.jpg','06_local_temperature_flux.jpg'),('07_Integrated_Model','10_transport.jpg','07_integrated_transport.jpg')]:
        src=ROOT.parent/stage/'figures'/filename;dst=FINAL/'figures'/target;shutil.copy2(src,dst)
        srcs.append(dict(source=str(src.relative_to(REPO)),destination=str(dst.relative_to(REPO)),sha256=hashlib.sha256(src.read_bytes()).hexdigest(),modification='Unmodified copy'))
    (R/'figure_provenance.json').write_text(json.dumps(srcs,indent=2)+'\n')
    print(values)

if __name__=='__main__':main()
