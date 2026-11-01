"""SI Gilmore-Tait mechanics with the unchanged MR2 polytropic pressure.

No temperature is calculated or transferred from any other trajectory.
The Tait exponent is the water reference in Denner (2021), Sec. 4;
B is derived from the accepted ambient density and sound speed, not fitted.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .upstream import mechanics as mr2

@dataclass(frozen=True)
class Tait:
    exponent: float = 7.15
    sound_multiplier: float = 1.

    @property
    def B(self):
        p = mr2.Parameters()
        return p.rho*(p.c*self.sound_multiplier)**2/self.exponent-p.p0

    def density(self, pressure):
        p = mr2.Parameters()
        if np.any(np.asarray(pressure)+self.B <= 0):
            raise ValueError('Outside the positive-density Tait branch')
        return p.rho*np.exp(np.log1p((np.asarray(pressure)-p.p0)/(p.p0+self.B))/self.exponent)

    def quantities(self, pl, pinf):
        rl, ri = self.density(pl), self.density(pinf)
        C = np.sqrt(self.exponent*(pl+self.B)/rl)
        # Stable enthalpy difference even close to equilibrium.
        q = (self.exponent-1)/self.exponent
        H = self.exponent/(self.exponent-1)*(pinf+self.B)/ri * np.expm1(q*np.log1p((pl-pinf)/(pinf+self.B)))
        return rl, ri, C, H

def terms(t, R, U, amplitude, eos=Tait()):
    p = mr2.Parameters()
    pg = mr2.gas_pressure(R,p)
    pinf = p.p0-amplitude*np.sin(2*np.pi*p.f*t)
    dpinf = -amplitude*2*np.pi*p.f*np.cos(2*np.pi*p.f*t)
    pl = pg+p.pv-2*p.sigma/R-4*p.mu*U/R
    bw = -3*p.kappa*pg*U/R+2*p.sigma*U/R**2+4*p.mu*U**2/R**2
    rl,ri,C,H = eos.quantities(pl,pinf)
    return pl,pinf,bw,dpinf,rl,ri,C,H

def gilmore_rhs(t, y, amplitude, eos=Tait()):
    R,U = y
    if R<=0: raise ValueError('Nonpositive trial radius')
    pl,pi,bw,dpi,rl,ri,C,H = terms(t,R,U,amplitude,eos)
    D = (1-U/C)*(R+4*mr2.Parameters().mu/(C*rl))
    N = (1+U/C)*H+(1-U/C)*R/C*(bw/rl-dpi/ri)-1.5*(1-U/(3*C))*U**2
    return np.array([U,N/D])

def implicit_residual(t,y,amplitude,eos=Tait()):
    R,U = y
    A = gilmore_rhs(t,y,amplitude,eos)[1]
    pl,pi,bw,dpi,rl,ri,C,H = terms(t,R,U,amplitude,eos)
    hd = (bw-4*mr2.Parameters().mu*A/R)/rl-dpi/ri
    left = (1-U/C)*R*A+1.5*(1-U/(3*C))*U**2
    right = (1+U/C)*H+(1-U/C)*R/C*hd
    return float(left-right), max(1.,abs(left),abs(right),abs(R*A),U**2)

def solve(amplitude_bar, model, method='DOP853', rtol=1e-9, steps=800):
    """One complete acoustic cycle with event and accepted-step dense controls."""
    p = mr2.Parameters(); amplitude = amplitude_bar*1e5
    fun = (lambda t,y: mr2.km_rhs(t,y,amplitude,p)) if model=='KM' else (lambda t,y:gilmore_rhs(t,y,amplitude))
    def minimum(t,y):return y[1]
    minimum.direction=1
    sol = solve_ivp(fun,(0,1/p.f),[p.R0,0.],method=method,rtol=rtol,
        atol=[rtol*1e-6,rtol*1e2],max_step=1/p.f/steps,dense_output=True,events=minimum)
    if not sol.success:raise RuntimeError(sol.message)
    # Locate extrema on accepted half intervals. Root tolerance is numerical only.
    grid = np.unique(np.r_[sol.t,(sol.t[:-1]+sol.t[1:])/2])
    def roots(component):
        vals = np.array([component(t) for t in grid]); ans=[]
        for a,b,fa,fb in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]):
            if fa*fb<0:ans.append(brentq(component,a,b,xtol=1e-22))
        return np.array(ans)
    rad = roots(lambda t:sol.sol(t)[1]); times=np.r_[0.,rad,sol.t[-1]]
    tspeed = roots(lambda t:fun(t,sol.sol(t))[1]); allspeed=np.r_[0.,tspeed,sol.t[-1]]
    rr = sol.sol(times)[0];j=np.argmin(rr)
    u = sol.sol(allspeed)[1]; allgrid=np.unique(np.r_[grid,times,allspeed])
    R,U=sol.sol(allgrid);pl=mr2.gas_pressure(R,p)+p.pv-2*p.sigma/R-4*p.mu*U/R
    row=dict(model=model,Pa_bar=amplitude_bar,method=method,rtol=rtol,steps_per_cycle=steps,
        completed=True,interval_s=sol.t[-1],Rmin_m=rr[j],t_Rmin_s=times[j],Rmax_m=np.max(rr),
        inward_speed_m_s=-np.min(u),max_speed_m_s=np.max(abs(u)),max_M_ambient=np.max(abs(u))/p.c,
        max_wall_pressure_Pa=np.max(pl),chi=np.max(abs(pl-p.p0))/(p.rho*p.c**2),
        minimum_farfield_pressure_Pa=p.p0-amplitude,temperature_recomputed=False)
    if model=='Gilmore':
        C=np.sqrt(Tait().exponent*(pl+Tait().B)/Tait().density(pl))
        row.update(max_M_local=float(np.max(abs(U)/C)),max_density_ratio=float(np.max(Tait().density(pl))/p.rho))
    return row,sol,allgrid
