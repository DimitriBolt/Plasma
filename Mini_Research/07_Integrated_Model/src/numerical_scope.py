"""Quantitative resolution statements derived from retained controls."""
from dataclasses import asdict
import json
import numpy as np
import pandas as pd
from .upstream import ROOT,optics,Parameters,BASE
from .experiments import DEFAULT
from .particle import Particle
R=ROOT/'results'

def main():
    conv=pd.read_csv(R/'convergence.csv',float_precision='round_trip')
    prod=conv[conv.check=='production'].iloc[0];fine=conv[conv.check=='fine'].iloc[0]
    native=pd.read_csv(R/'independent_Radau.csv',float_precision='round_trip').iloc[0]
    stable=pd.read_csv(R/'stable_pairs.csv',float_precision='round_trip')
    primary=stable[(stable.Pa_bar==1)&(stable.level==3)&(stable.order==3)&(stable.fluence_J_m2==10)].iloc[0]
    hm=pd.read_csv(R/'high_mesh_sensitivity.csv',float_precision='round_trip')
    timing=pd.read_csv(R/'timing_sensitivity.csv',float_precision='round_trip')
    high=pd.read_csv(R/'high_Radau_sensitivity.csv',float_precision='round_trip').iloc[0]
    rows=[]
    for r in hm.to_dict('records'):
        c=timing[(timing.Pa_bar==1.15)&(timing.level==3)&np.isclose(timing.offset_s,r['offset_s'],atol=1e-15)].iloc[0]
        diff=abs(r['dTmax_dF_K_m2_J']/c.dTmax_dF_K_m2_J-1)
        rows.append(dict(Pa_bar=1.15,offset_s=r['offset_s'],cells48_slope=c.dTmax_dF_K_m2_J,
            cells96_slope=r['dTmax_dF_K_m2_J'],relative_mesh_change=diff,
            numerical_resolution='unresolved_late_tail' if diff>.1 else 'empirically_resolved_reduced_model'))
    pd.DataFrame(rows).to_csv(R/'strong_collapse_resolution.csv',index=False,float_format='%.17g')
    early=timing[(timing.Pa_bar==1.15)&(timing.level==3)&np.isclose(timing.offset_s,-5e-6,atol=1e-15)].iloc[0]
    obj=dict(one_bar_mesh48_to96_relative_change=abs(fine.dTmax_dF_K_m2_J/prod.dTmax_dF_K_m2_J-1),
        one_bar_independent_finite_Radau_relative_difference=abs(native.global_delta_Tmax_K/primary.global_delta_Tmax_K-1),
        one_bar_full_baseline_temperature_spread_K=float(conv.observed_Tmax_K.max()-conv.observed_Tmax_K.min()),
        bar115_Radau_directional_relative_difference=abs(high.dTmax_dF_K_m2_J/early.dTmax_dF_K_m2_J-1),
        late_centered_pulse='unresolved: large 48-to-96-cell variation, values retained',
        physical_validation='not established; local Au field, force, chemistry, EOS and high-Mach limits remain')
    (R/'numerical_scope.json').write_text(json.dumps(obj,indent=2)+'\n')
    params=dict(model_id='moving-material-radial-KM-NASA9-Au-v1',accepted_input_sha=BASE,
        production_config=asdict(DEFAULT),particle=Particle().metadata(),water_mechanics=asdict(Parameters()),
        thermal_properties=asdict(optics.ThermalProperties()),
        acoustic_amplitudes_bar=[1,1.15,1.3],pulse_F_J_m2=[0,10,100],pulse_duration_s=1e-8,
        nominal_pulse_offset_s=-5e-6,formal_dry_case_offset_s=-1e-6,
        timing_offsets_s=[-5e-6,-1e-6,-1e-7,0,1e-7,1e-6],
        local_counts=[0,1,8],local_count_is_number_density=False,
        gas_model='closed fixed ideal mixture; no additional constant vapor pressure',
        gas_heat_flow='single current-state signed Robin flow; no additive MR4 or MR6 baseline',
        chemistry_and_phase_changes='not evolved; separate reference scales only',
        wet_domain_K=[200,6000],dry_air_domain_K=[200,20000],dry_ground_Ar_domain='positive T, constant Cv; formal at high T',
        MR6_source_contract='frozen-reference-baseline-plus-laser-v2; no frozen trace used as dynamic input')
    (R/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    print(json.dumps(obj,indent=2))

if __name__=='__main__':main()
