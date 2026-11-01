"""Matched stationary MR6 comparison and analytic radial controls, no fitting."""
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from .transport import MaterialGrid
from .upstream import ROOT,optics

def stationary(R=6e-6,a=20e-9,d=100e-9,h=4400.,G=1e8,k=.6,
               cells=160,padding=25e-6,energy=3.4337324474063743e-15,duration=10e-9,end=2e-6):
    grid=MaterialGrid(cells,padding,k,initial_radius_m=R)
    w=grid.surface_weights(R,R+a+d,a)
    ki,kk,kb=grid.conductances(R);H=4*np.pi*R**2*h;H=H*ki/(H+ki)
    C=4*np.pi*a**3/3*19320*129;Kau=4*np.pi*a*a*G
    # Assemble the physical linear operator by its sparse pairwise exchanges.
    n=cells+1;L=np.zeros((n,n));cap=np.r_[grid.capacity,C]
    for j,K in enumerate(kk):
        L[j,j]-=K;L[j+1,j+1]-=K;L[j,j+1]+=K;L[j+1,j]+=K
    L[0,0]-=H;L[cells-1,cells-1]-=kb
    # QAu=Kau*(TAu-w.T Tl); water receives w QAu, preserving reciprocity.
    L[:cells,:cells]-=Kau*np.outer(w,w)
    L[:cells,-1]+=Kau*w;L[-1,:cells]+=Kau*w;L[-1,-1]-=Kau
    A=L/cap[:,None]
    full=np.zeros((n+2,n+2));full[:n,:n]=A
    full[-2,0]=H;full[-1,cells-1]=kb
    src=np.zeros(n+2);src[cells]=energy/duration/C
    t=np.unique(np.r_[0.,np.geomspace(1e-13,end,800),duration])
    s1=solve_ivp(lambda t,y:full@y+src,(0,duration),np.zeros(n+2),method='BDF',jac=full,
                 rtol=1e-10,atol=np.r_[np.full(n,1e-12),1e-29,1e-29],dense_output=True)
    s2=solve_ivp(lambda t,y:full@y,(duration,end),s1.y[:,-1],method='BDF',jac=full,
                 rtol=1e-10,atol=np.r_[np.full(n,1e-12),1e-29,1e-29],dense_output=True)
    y=np.empty((n+2,len(t)));mask=t<=duration
    y[:,mask]=s1.sol(t[mask]);y[:,~mask]=s2.sol(t[~mask])
    Eabs=energy*np.minimum(t/duration,1)
    stored=cap@y[:n]
    return dict(time_s=t,water_rise_K=y[:cells],gold_rise_K=y[cells],
        q_gas_W=H*y[0],gas_energy_J=y[-2],bath_energy_J=y[-1],
        stored_energy_J=stored,balance_residual_J=stored+y[-2]+y[-1]-Eabs,
        xi_faces_m3=grid.faces,capacity_J_K=cap)

def main():
    refroot=ROOT.parent/'06_AuNP_Bubble_Coupling/results'
    settings=[('reference',refroot/'reference.npz',6e-6,20e-9,100e-9,4400.),
              ('far_particle',refroot/'fields/far_particle.npz',6e-6,20e-9,500e-9,4400.),
              ('near_contact',refroot/'fields/near_contact.npz',6e-6,100e-9,1e-9,1e6)]
    rows=[]
    for name,path,R,a,d,h in settings:
        z=np.load(path);E=float(z['absorbed_so_far_J'][-1])
        for cells in (80,160,320):
            v=stationary(R,a,d,h,cells=cells,energy=E)
            np.savez_compressed(ROOT/'results/raw'/f'static_{name}_{cells}.npz',**v)
            rows.append(dict(case=name,cells=cells,R_m=R,a_m=a,d_m=d,h_W_m2K=h,G_W_m2K=1e8,
              k_water_W_mK=.6,fluence_J_m2=1.,duration_s=10e-9,energy_J=E,
              radial_gas_fraction=v['gas_energy_J'][-1]/E,
              MR6_gas_fraction=z['gas_energy_J'][-1]/E,
              relative_integrated_transfer_difference=v['gas_energy_J'][-1]/z['gas_energy_J'][-1]-1,
              radial_gold_peak_K=v['gold_rise_K'].max(),MR6_gold_peak_K=z['gold_rise_K'].max(),
              radial_wall_power_max_W=v['q_gas_W'].max(),MR6_wall_power_max_W=z['q_wall_laser_increment_W'].max(),
              balance_over_absorbed=np.max(abs(v['balance_residual_J']))/E))
    pd.DataFrame(rows).to_csv(ROOT/'results/static_MR6_comparison.csv',index=False,float_format='%.17g')
    print(pd.DataFrame(rows)[['case','cells','relative_integrated_transfer_difference','radial_gold_peak_K','MR6_gold_peak_K']].to_string(index=False),flush=True)

if __name__=='__main__':main()
