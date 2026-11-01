"""Delay and compression diagnostics from saved nonlinear pairs, without ODE replay."""
import numpy as np
import pandas as pd
from .upstream import ROOT,T0
from .stable_pairs import archived
R=ROOT/'results'

def main():
    data={};rows=[]
    for pa,tag in [(1.,'100'),(1.15,'115')]:
        name=f'base_{pa:.2f}';base,_=archived(name);m=base.model
        z=np.load(R/'raw'/(f'stable_{name}_F10_L3_q3.npz'))
        meta=pd.read_csv(R/'stable_pairs.csv',float_precision='round_trip')
        r=meta[meta.run_id==f'stable_{name}_F10_L3_q3'].iloc[0];center=r.pulse_center_s
        ids=np.flatnonzero((z['t_s']>=center-10e-9)&(z['t_s']<=center+8e-6))
        t=z['t_s'][ids];yb=base.sol(t)
        q0=np.array([m.powers(tt,yb[:,j])[0] for j,tt in enumerate(t)])
        dq=z['Q_l_to_g_W'][ids]-q0;dta=z['T_Au_rise_K'][ids]-yb[m.au]*T0
        for key,val in [('time_from_center_s',t-center),('delta_Qgas_W',dq),
                        ('delta_Tgas_K',z['paired_delta_T_K'][ids]),('delta_TAu_K',dta)]:data[tag+'_'+key]=val
        pre=np.flatnonzero((t>=center)&(t<=center+2e-6));j=pre[np.argmax(dq[pre])]
        raw=np.load(R/'raw'/(name+'.npz'));tt=raw['peak_candidates_s'];tv=raw['peak_candidates_K'];mask=tt<=1/30000
        tp=float(tt[mask][np.argmax(tv[mask])])
        rows.append(dict(Pa_bar=pa,pulse_center_s=center,observed_incremental_flux_peak_lag_s=t[j]-center,
             incremental_flux_peak_W=dq[j],flux_peak_window_s=2e-6,
             baseline_first_cycle_thermal_peak_after_center_s=tp-center,
             first_cycle_delta_T_at_reference_K=np.interp(tp,z['t_s'],z['paired_delta_T_K']),
             meaning='net wall-flux increment includes gas feedback; not a pure laser-tagged delivery fraction'))
    np.savez_compressed(R/'transport_timing.npz',**data)
    pd.DataFrame(rows).to_csv(R/'transport_timing.csv',index=False,float_format='%.17g')
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
