"""Separate outer-domain-shape context from angular source error."""
import numpy as np
import pandas as pd
from .upstream import ROOT
from .benchmarks import stationary

def main():
    old=pd.read_csv(ROOT/'results/static_MR6_comparison.csv',float_precision='round_trip')
    rows=[]
    for r in old[old.cells==320].to_dict('records'):
        padding=1e-6+2*r['a_m']+r['d_m']
        out=stationary(r['R_m'],r['a_m'],r['d_m'],r['h_W_m2K'],cells=320,padding=padding,energy=r['energy_J'])
        np.savez_compressed(ROOT/'results/raw'/('static_'+r['case']+'_matched_clearance.npz'),**out)
        rows.append(dict(case=r['case'],outer_radial_padding_m=padding,MR6_outer_padding_m=1e-6,
          nearest_Au_to_radial_bath_m=1e-6,MR6_nearest_Au_to_bath_m=1e-6,
          bath_shape_note='nearest Au-bath clearance matched; sphere versus cylinder shapes differ',
          radial_gas_fraction=out['gas_energy_J'][-1]/r['energy_J'],MR6_gas_fraction=r['MR6_gas_fraction'],
          relative_difference=out['gas_energy_J'][-1]/(r['energy_J']*r['MR6_gas_fraction'])-1,
          radial_gold_peak_K=out['gold_rise_K'].max(),MR6_gold_peak_K=r['MR6_gold_peak_K']))
    pd.DataFrame(rows).to_csv(ROOT/'results/static_domain_shape.csv',index=False,float_format='%.17g')
    print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':main()
