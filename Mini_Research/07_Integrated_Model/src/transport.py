"""Material volumes, exact radial resistances and normalized finite-surface source."""
import numpy as np
from .upstream import Parameters, optics

PI = np.pi

def selected_max(x, y):
    return np.where(np.real(x) >= np.real(y), x, y)

def selected_min(x, y):
    return np.where(np.real(x) <= np.real(y), x, y)

class MaterialGrid:
    def __init__(self, cells=64, padding_m=25e-6, k=.6, rho_cp=1000*4187,
                 initial_radius_m=Parameters().R0,faces_m3=None):
        if cells < 4 or min(padding_m,k,rho_cp,initial_radius_m)<=0:
            raise ValueError('Positive grid and properties required')
        self.k, self.rho_cp, self.R0 = k, rho_cp, initial_radius_m
        # Fixed material coordinates, initially logarithmic physical wall gaps.
        gaps=np.expm1(np.linspace(0,np.log1p(padding_m/5e-9),cells+1))*5e-9
        self.faces=((self.R0+gaps)**3-self.R0**3)/3 if faces_m3 is None else np.asarray(faces_m3,dtype=float)
        if len(self.faces)!=cells+1 or self.faces[0]!=0 or not np.all(np.diff(self.faces)>0):
            raise ValueError('Material faces must start at zero and increase strictly')
        self.centers=(self.faces[1:]+self.faces[:-1])/2
        self.capacity=4*PI*rho_cp*np.diff(self.faces)
        self.cells=cells

    def radii(self,R):
        # Principal branch is analytic around every admitted positive real state.
        return (R**3+3*self.centers)**(1/3), (R**3+3*self.faces)**(1/3)

    def conductances(self,R):
        r,rf=self.radii(R)
        inner=4*PI*self.k/(1/R-1/r[0])
        interior=4*PI*self.k/(1/r[:-1]-1/r[1:])
        outer=4*PI*self.k/(1/r[-1]-1/rf[-1])
        return inner,interior,outer

    def surface_weights(self,R,rnp,a):
        r,rf=self.radii(R)
        # Integrate radial Au-surface density against linear water reconstruction.
        # Deposit with its adjoint weights for exact thermal conservation.
        low=selected_max(r[:-1],rnp-a);high=selected_min(r[1:],rnp+a)
        width=high-low;dr=r[1:]-r[:-1];denom=2*rnp*a*dr
        left=width*((r[1:]-high)*(high+low)/2+width*(high+2*low)/6)/denom
        right=width*((low-r[:-1])*(high+low)/2+width*(2*high+low)/6)/denom
        mask=np.real(width)>0
        w=np.zeros(self.cells,dtype=np.result_type(R,rnp))
        w[:-1]+=np.where(mask,left,0.);w[1:]+=np.where(mask,right,0.)
        lo=selected_max(rf[0],rnp-a);hi=selected_min(r[0],rnp+a)
        w[0]+=np.where(np.real(hi)>np.real(lo),(hi-lo)*(hi+lo)/(4*rnp*a),0.)
        lo=selected_max(r[-1],rnp-a);hi=selected_min(rf[-1],rnp+a)
        w[-1]+=np.where(np.real(hi)>np.real(lo),(hi-lo)*(hi+lo)/(4*rnp*a),0.)
        # No renormalization of an escaped source: the geometry guard stops it.
        return w

    def water_rate(self,R,theta,q_gas,q_au,weights):
        _,ki,kb=self.conductances(R)
        flux=ki*(theta[:-1]-theta[1:])
        bath=kb*theta[-1]
        power=np.zeros(self.cells,dtype=np.result_type(theta,R))
        power[:-1]-=flux; power[1:]+=flux
        power[0]-=q_gas; power[-1]-=bath
        power+=weights*q_au
        return power/self.capacity,bath

def center_radius(R,a,d,mode,R0=Parameters().R0):
    r0=R0+a+d
    if mode=='passive': return (R**3+r0**3-R0**3)**(1/3)
    if mode=='fixed': return r0+0*R
    if mode=='capture': return R+a+d
    raise ValueError('Unknown motion')
