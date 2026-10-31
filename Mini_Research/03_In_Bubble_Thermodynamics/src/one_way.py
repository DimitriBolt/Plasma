"""One-way post-processing of the saved, newly reconstructed mechanics."""
import numpy as np
from .bubble_thermo import Gas,Parameters,adiabatic_temperature,gas_pressure,temperature_rate
from .heat_transfer import heat_loss
from .heat_source import heat_source
from .numerics import integrate,Settings,InvalidState


def solve_one_way(path,gas=Gas(),eta=1.,p=Parameters(),settings=Settings()):
    # theta=T/Tad(R): exact cancellation of reversible compression, with U=dR/dt.
    def fun(t,y):
        R,U=path(t);Tad=adiabatic_temperature(R,gas,p);T=y[0]*Tad
        if T<=0 or R<=0 or not np.all(np.isfinite(y)):
            raise InvalidState('Nonpositive one-way state')
        loss=heat_loss(R,T,gas,eta);source=heat_source(t,R,T)
        work=gas_pressure(R,T,gas,p)*4*np.pi*R**2*U/gas.E0(p)
        return np.array([(source-loss)/(gas.capacity(p)*Tad),work,
                         loss/gas.E0(p),source/gas.E0(p)])
    sol=integrate(fun,[1,0,0,0],path.end,1/p.f,[1e-9,1e-8,1e-8,1e-8],settings)
    sol.path=path;sol.kind='one_way_post_processing';sol.fun=fun
    sol.gas,sol.eta,sol.amplitude,sol.p=gas,eta,path.amplitude,p
    return sol


def temperature(sol,t):
    return sol.sol(t)[0]*adiabatic_temperature(sol.path(t)[0],sol.gas,sol.p)


def temperature_derivative(sol,t):
    R,U=sol.path(t);T=temperature(sol,t)
    return temperature_rate(R,U,T,heat_source(t,R,T)-heat_loss(R,T,sol.gas,sol.eta),sol.gas,sol.p)
