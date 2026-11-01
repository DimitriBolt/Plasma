"""Exact equation-limit regression, with independent first-cycle and global extrema."""
from dataclasses import replace
import numpy as np
import pandas as pd
from .upstream import ROOT
from .model import Model,Config
from .particle import Particle
from .diagnostics import extrema

def main():
    source=ROOT.parent/'04_Temperature_Barrier/results';df=pd.read_csv(source/'summary.csv',float_precision='round_trip')
    rows=[]
    for pa in (1.,1.15,1.3):
        m=Model(Config(amplitude_bar=pa,cells=8,clamp_water=True,formal_material=True,method='Radau',rtol=2e-11,atol=1e-13),Particle(count=0))
        s=m.solve();e=extrema(s);key=f'M2_sat_air_{pa:.2f}bar_coupled_eta1'
        r=df[df.run_id==key].iloc[0];z=np.load(source/'trajectories'/(key+'.npz'))
        end=min(s.end,1/m.p.f);first=extrema(s,end)
        full_first=bool(min(s.end,r.end_s)>=1/m.p.f)
        mask=z['t_s']<=end;ids=np.flatnonzero(mask);j=ids[np.argmax(z['T_K'][mask])]
        rows.append(dict(Pa_bar=pa,MR7_Tmax_K=e['observed_Tmax_K'],MR4_Tmax_K=r.observed_Tmax_K,
          peak_relative_difference=e['observed_Tmax_K']/r.observed_Tmax_K-1,
          MR7_t_Tmax_s=e['t_Tmax_s'],MR4_t_Tmax_s=r.t_Tmax_s,
          global_peak_time_difference_s=e['t_Tmax_s']-r.t_Tmax_s,
          first_interval_end_s=end,first_MR7_Tmax_K=first['observed_Tmax_K'],first_MR4_Tmax_K=z['T_K'][j],
          first_cycle_comparison_defined=full_first,
          first_peak_relative_difference=first['observed_Tmax_K']/z['T_K'][j]-1 if full_first else np.nan,
          first_peak_time_difference_s=first['t_Tmax_s']-z['t_s'][j] if full_first else np.nan,
          MR7_completed=s.completed,MR4_completed=bool(r.completed),status=s.status,
          MR7_end_s=s.end,MR4_end_s=r.end_s,observed_end_time_difference_s=s.end-r.end_s))
        pd.DataFrame(rows).to_csv(ROOT/'results/MR4_regression_full.csv',index=False,float_format='%.17g')
        print(rows[-1],flush=True)

if __name__=='__main__':main()
