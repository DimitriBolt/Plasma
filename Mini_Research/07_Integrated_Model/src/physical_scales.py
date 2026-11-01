"""Missing-channel reference demands; not a reacting or phase-change temperature model."""
import json
import numpy as np
import pandas as pd
from .upstream import ROOT,thermo,RU
from .scientific_audit import load_model

R=ROOT/'results'

def main():
    ref=json.loads((ROOT.parent/'04_Temperature_Barrier/results/chemical_parameter_sources.json').read_text())
    hf=ref['formation_enthalpy_J_mol'];Tr=298.15
    phase=dict(reference_temperature_K=Tr,Hf_H2O_liquid_J_mol=-285830.,Hf_H2O_gas_J_mol=hf['H2O'],
       liquid_source='https://janaf.nist.gov/tables/H-063.html',
       gas_source='https://webbook.nist.gov/cgi/cbook.cgi?Name=water&cTG=on',
       interpretation='Standard reference demands, not high-pressure latent heat or a reaction extent. No temperature recomputed.')
    (R/'missing_channel_parameters.json').write_text(json.dumps(phase,indent=2)+'\n')
    rows=[];screen=[]
    names=['base_1.00','base_1.15','base_1.30','dry_air_base','dry_argon_base']
    for name in names:
        p=R/'raw'/(name+'.json')
        if not p.exists():continue
        meta=json.loads(p.read_text());m=load_model(meta);s=meta['summary'];z=np.load(p.with_suffix('.npz'))
        n=m.mix.nt;fr=m.mix.fractions;nw=n*fr.get('H2O',0.)
        demands={'H2O_vaporization_magnitude_298K':nw*(hf['H2O']+285830-RU*Tr),
            'H2O_to_H_OH_298K':nw*(hf['H']+hf['OH']-hf['H2O']-RU*Tr),
            'H2O_to_2H_O_298K':nw*(2*hf['H']+hf['O']-hf['H2O']-2*RU*Tr),
            'N2_to_2N_298K':n*fr.get('N2',0)*(2*hf['N']-hf['N2']-RU*Tr),
            'O2_to_2O_298K':n*fr.get('O2',0)*(2*hf['O']-hf['O2']-RU*Tr),
            'Ar_first_ionization_ground_threshold':n*fr.get('Ar_ground',0)*ref['first_ionization_eV']['Ar']*96485.33212331001}
        Epeak=float(np.max(z['U_gas_change_J']))
        for channel,E in demands.items():
            if E:
                rows.append(dict(run_id=name,process=channel,reference_energy_J=E,
                    ratio_to_frozen_peak_sensible_energy=E/Epeak,
                    no_reaction_or_phase_temperature_recomputed=True,
                    overlapping_processes_not_additive=True))
        t=z['t_s'];T=z['T_g_K'];j=int(np.argmax(T));peak=T[j];half=thermo.T0+(peak-thermo.T0)/2
        left=np.flatnonzero(T[:j]<half);right=np.flatnonzero(T[j+1:]<half)+j+1
        width=np.nan
        if len(left) and len(right):
            l=left[-1];r=right[0]
            tl=np.interp(half,T[l:l+2],t[l:l+2]);tr=np.interp(half,T[r-1:r+1][::-1],t[r-1:r+1][::-1]);width=tr-tl
        screen.append(dict(run_id=name,Pa_bar=m.cfg.amplitude_bar,gas=m.cfg.gas,vapor_factor=m.cfg.vapor_factor,
            observed_Tmax_K=s['observed_Tmax_K'],global_Tmax_K=s['global_Tmax_K'],status=s['status'],
            max_wall_Mach=s['max_wall_Mach'],max_pressure_Pa=s['max_pressure_Pa'],
            max_gas_density_kg_m3=s['max_gas_density_kg_m3'],
            pressure_over_liquid_bulk_modulus=s['pressure_over_liquid_bulk_modulus'],
            gas_acoustic_crossing_at_peak_s=s['gas_acoustic_crossing_at_peak_s'],
            reconstructed_temperature_half_height_width_s=width,
            sound_crossing_over_thermal_width=s['gas_acoustic_crossing_at_peak_s']/width if width>0 else np.nan))
    pd.DataFrame(rows).to_csv(R/'missing_channel_scales.csv',index=False,float_format='%.17g')
    pd.DataFrame(screen).to_csv(R/'physical_diagnostics.csv',index=False,float_format='%.17g')
    print('Saved reference energy demands and independent EOS/time-scale diagnostics.')

if __name__=='__main__':main()
