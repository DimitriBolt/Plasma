"""Transient point-shell control against the exact infinite-liquid Robin sphere."""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.special import erfc,erfcx
from .transport import MaterialGrid
from .upstream import ROOT

def exact_fraction(t,R,r,alpha,k,h):
    t=np.asarray(t);ans=np.zeros_like(t);mask=t>0
    beta=1/R+h/k;x=(r-R)/np.sqrt(4*alpha*t[mask])
    ans[mask]=R/r*(h/k)/beta*(erfc(x)-np.exp(-x*x)*erfcx(x+beta*np.sqrt(alpha*t[mask])))
    return ans

def control(cells):
    R=6e-6;r=R+520e-9;k=.6;h=4400.;E=1e-15
    g=MaterialGrid(cells,25e-6,k,initial_radius_m=R);rc,_=g.radii(R)
    w=np.zeros(cells);j=np.searchsorted(rc,r);f=(r-rc[j-1])/(rc[j]-rc[j-1]);w[j-1]=1-f;w[j]=f
    ki,kk,kb=g.conductances(R);H=4*np.pi*R*R*h;H=H*ki/(H+ki)
    A=np.zeros((cells+2,cells+2))
    L=np.zeros((cells,cells))
    for j,K in enumerate(kk):
        L[j,j]-=K;L[j+1,j+1]-=K;L[j,j+1]+=K;L[j+1,j]+=K
    L[0,0]-=H;L[-1,-1]-=kb
    A[:cells,:cells]=L/g.capacity[:,None];A[-2,0]=H;A[-1,cells-1]=kb
    y=np.zeros(cells+2);y[:cells]=w*E/g.capacity
    times=np.array([0.,.2e-6,.5e-6,1e-6,2e-6])
    sol=solve_ivp(lambda t,y:A@y,(0,times[-1]),y,method='BDF',jac=A,rtol=1e-10,
                 atol=np.r_[np.full(cells,1e-14),1e-29,1e-29],t_eval=times)
    ex=exact_fraction(times,R,r,k/(1000*4187),k,h)
    rows=[dict(cells=cells,t_s=t,radial_fraction=val,exact_fraction=e,
               relative_error=val/e-1 if e else np.nan,
               balance_over_initial=(g.capacity@sol.y[:cells,i]+sol.y[-2,i]+sol.y[-1,i]-E)/E)
          for i,(t,val,e) in enumerate(zip(times,sol.y[-2]/E,ex))]
    return rows

def main():
    rows=[]
    for n in (80,160,320):rows+=control(n)
    pd.DataFrame(rows).to_csv(ROOT/'results/analytic_point_control.csv',index=False,float_format='%.17g')
    print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':main()
