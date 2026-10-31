"""Fully coupled Keller--Miksis + fixed-inventory gas energy balance."""
import numpy as np
from .bubble_thermo import Parameters,Gas,volume,gas_pressure,temperature_rate,pressure_rate
from .heat_transfer import heat_loss
from .heat_source import heat_source
from .numerics import integrate,Settings,InvalidState


def acceleration(t,R,U,T,Pa,gas,eta,p):
    pg=gas_pressure(R,T,gas,p)
    net=heat_source(t,R,T)-heat_loss(R,T,gas,eta)
    dpg=pressure_rate(R,U,T,net,gas,p)
    dp=pg+p.pv-(p.p0-Pa*np.sin(2*np.pi*p.f*t))-2*p.sigma/R-4*p.mu*U/R
    B=dpg+2*p.sigma*U/R**2+Pa*2*np.pi*p.f*np.cos(2*np.pi*p.f*t)+4*p.mu*U**2/R**2
    D=R*(1-U/p.c)+4*p.mu/(p.rho*p.c)
    if D==0:
        raise InvalidState('KM denominator is zero')
    return ((1+U/p.c)*dp/p.rho+R*B/(p.rho*p.c)-1.5*(1-U/(3*p.c))*U**2)/D


def rhs(t,y,Pa,gas=Gas(),eta=1.,p=Parameters()):
    R,U,T=y[:3]
    if R<=0 or T<=0 or not np.all(np.isfinite(y)):
        raise InvalidState('Nonpositive R or T, or nonfinite trial state')
    loss=heat_loss(R,T,gas,eta);source=heat_source(t,R,T)
    A=acceleration(t,R,U,T,Pa,gas,eta,p)
    dT=temperature_rate(R,U,T,source-loss,gas,p)
    work=gas_pressure(R,T,gas,p)*4*np.pi*R**2*U/gas.E0(p)
    # Integrated variables are dimensionless W/E0, signed Qloss/E0, Qsource/E0.
    return np.array([U,A,dT,work,loss/gas.E0(p),source/gas.E0(p)])


def solve_coupled(amplitude,gas=Gas(),eta=1.,p=Parameters(),settings=Settings(),cycles=5.):
    fun=lambda t,y:rhs(t,y,amplitude*1e5,gas,eta,p)
    guards=((lambda y:y[0]-1e-12,'technical_radius_guard_1pm'),
            (lambda y:1-y[1]/p.c+4*p.mu/(p.rho*p.c*y[0])-1e-6,'technical_denominator_guard'))
    sol=integrate(fun,[p.R0,0,gas.T0,0,0,0],cycles/p.f,1/p.f,
                  [1e-16,1e-8,1e-6,1e-8,1e-8,1e-8],settings,guards)
    sol.fun=fun
    sol.kind='fully_coupled'
    sol.gas,sol.eta,sol.amplitude,sol.p=gas,eta,amplitude,p
    return sol
