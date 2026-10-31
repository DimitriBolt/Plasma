"""Energy scales only: no composition evolution and no corrected temperature."""
import math
import pandas as pd
from .upstream import ROOT,RU
from .run_experiments import CONFIGS,write_json
from .thermodynamics import Mixture

F=96485.33212331002  # NA * exact elementary charge, J mol^-1 eV^-1
TF=298.15
# NIST WebBook CODATA review values, gas, 298.15 K; OH: Chase (1998).
HF={'H':217998.,'O':249180.,'N':472680.,'OH':38990.,'H2O':-241826.,'N2':0.,'O2':0.}
IE={'Ar':15.759610,'H':13.598433,'O':13.61805,'N':14.5341}

def reactions():
    rows=[]
    for name,dH,dn,parent in [('N2 -> 2N',2*HF['N'],1,'N2'),('O2 -> 2O',2*HF['O'],1,'O2'),
       ('H2O -> OH + H',HF['OH']+HF['H']-HF['H2O'],1,'H2O'),
       ('H2O -> 2H + O',2*HF['H']+HF['O']-HF['H2O'],2,'H2O')]:
        rows.append(dict(process=name,parent=parent,type='dissociation_298K_energy_scale',
           delta_H_298_J_mol=dH,delta_U_298_J_mol=dH-dn*RU*TF,energy_eV=dH/F,
           equivalent_E_over_k_K=dH/RU,energy_scale_reference_K=TF))
    for s,e in IE.items():
        rows.append(dict(process=s+' -> '+s+'+ + e-',parent=s,type='ground_state_first_ionization',
            delta_H_298_J_mol=math.nan,delta_U_298_J_mol=e*F,energy_eV=e,
            equivalent_E_over_k_K=e*F/RU,energy_scale_reference_K=0.))
    return pd.DataFrame(rows)

def main():
    out=ROOT/'results';rr=reactions();rr.to_csv(out/'reaction_energy_scales.csv',index=False)
    df=pd.read_csv(out/'summary.csv',float_precision='round_trip');rows=[]
    for m in df[(df.Pa_bar==1.3)&(df.config.isin(['M1_variable','M2_sat','M2_v10']))&(df.eta==1)].to_dict('records'):
        mix=Mixture(gas=m['gas'],**CONFIGS[m['config']]);nn={s:x*mix.nt for s,x in mix.fractions.items()}
        nn['Ar']=nn.get('Ar_ground',0.)
        # Hypothetical atomic inventories after complete dissociation; not actual populations.
        nn['H']=2*nn.get('H2O',0.);nn['O']=2*nn.get('O2',0.)+nn.get('H2O',0.);nn['N']=2*nn.get('N2',0.)
        for r in rr.to_dict('records'):
            n=nn.get(r['parent'],0.)
            if not n: continue
            E=n*r['delta_U_298_J_mol']
            rows.append(dict(run_id=m['run_id'],process=r['process'],hypothetical_reacting_mol=n,
                energy_for_complete_conversion_J=E,ratio_to_gross_compression_at_peak=E/m['gross_compression_peak_J'],
                ratio_to_collapse_compression=E/m['collapse_compression_J'],ratio_to_peak_delta_U=E/m['delta_U_peak_J'],
                frozen_observed_Tmax_K=m['observed_Tmax_K'],completed=m['completed'],
                no_reaction_temperature_recomputed=True))
    pd.DataFrame(rows).to_csv(out/'chemistry_inventory_budgets.csv',index=False)
    boltz=[]
    for r in rr.to_dict('records'):
        for T in (1000.,3000.,6000.,10000.,20000.,30000.):
            boltz.append(dict(process=r['process'],T_K=T,energy_scale_exp_factor=math.exp(-r['energy_eV']*F/(RU*T)),
                             interpretation='Energy-scale diagnostic only; neither equilibrium fraction nor rate'))
    pd.DataFrame(boltz).to_csv(out/'chemical_temperature_scales.csv',index=False)
    write_json(out/'chemical_parameter_sources.json',dict(formation_enthalpy_J_mol=HF,temperature_K=TF,
        formation_sources={s:'https://webbook.nist.gov/cgi/cbook.cgi?ID='+id+'&Units=SI&Mask=1' for s,id in
          [('H','C12385136'),('O','C17778802'),('N','C17778880'),('OH','C3352576'),('H2O','C7732185')]},
        first_ionization_eV=IE,ionization_sources={s:'https://physics.nist.gov/PhysRefData/Handbook/Tables/'+name+'table1.htm' for s,name in
          [('Ar','argon'),('H','hydrogen'),('O','oxygen'),('N','nitrogen')]},
        description='Dissociation costs are standard reaction enthalpy/internal energy at 298.15 K, not D0 or activation barriers. Ionization costs are isolated-atom ground-state thresholds. Not simultaneous additive channels.'))

if __name__=='__main__': main()
