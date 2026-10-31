"""KM plus closed-mixture energy; one-way uses the identical energy RHS."""
from dataclasses import dataclass
from types import SimpleNamespace
import math
import numpy as np
from scipy.integrate import DOP853, Radau, OdeSolution
from scipy.optimize import brentq
from .upstream import Parameters
from .thermodynamics import Mixture,T0

class InvalidState(ValueError): pass

class PositiveDOP853(DOP853):
    def _step_impl(self):
        for _ in range(40):
            try: return super()._step_impl()
            except InvalidState: self.h_abs*=.25
        return False,'Trial stages left positive domain'

@dataclass(frozen=True)
class Settings:
    method: str = 'DOP853'
    rtol: float = 2e-9
    steps_per_period: int = 600

def acceleration(t,R,U,T,amplitude,mix,eta,p=Parameters()):
    pg=mix.pressure(R,T)
    dT=mix.temperature_rate(R,U,T,eta,trial=True)
    dpg=mix.pressure_rate(R,U,T,dT)
    dp=pg+mix.p_constant-(p.p0-amplitude*1e5*math.sin(2*math.pi*p.f*t))-2*p.sigma/R-4*p.mu*U/R
    B=dpg+2*p.sigma*U/R**2+amplitude*1e5*2*math.pi*p.f*math.cos(2*math.pi*p.f*t)+4*p.mu*U**2/R**2
    D=R*(1-U/p.c)+4*p.mu/(p.rho*p.c)
    if D==0: raise InvalidState('Zero KM denominator')
    return ((1+U/p.c)*dp/p.rho+R*B/(p.rho*p.c)-1.5*(1-U/(3*p.c))*U**2)/D

def integrate(fun,y0,end,settings,guards,atol):
    p=Parameters();cls=PositiveDOP853 if settings.method=='DOP853' else Radau
    options={}
    if cls is Radau:
        scales=np.maximum(np.abs(y0),1.);scales[0]=abs(y0[0])
        def jac(t,y):
            hs=np.sqrt(np.finfo(float).eps)*np.maximum(np.abs(y),scales)
            columns=[]
            for j,h in enumerate(hs):
                step=np.zeros(len(y));step[j]=h
                columns.append((fun(t,y+step)-fun(t,y-step))/(2*h))
            return np.array(columns).T
        options['jac']=jac
    solver=cls(fun,0.,np.array(y0,dtype=float),end,rtol=settings.rtol,
               atol=atol,max_step=1/p.f/settings.steps_per_period,**options)
    ts=[0.];ys=[np.array(y0)];polys=[];status=0;message='Reached five-cycle end'
    try:
        while solver.status=='running':
            msg=solver.step()
            if solver.status=='failed': status=-1;message=str(msg);break
            dense=solver.dense_output();stops=[]
            for name,guard in guards:
                if guard(solver.t,dense(solver.t))<=0:
                    root=brentq(lambda t:guard(t,dense(t)),ts[-1],solver.t,xtol=1e-24,rtol=4*np.finfo(float).eps)
                    stops.append((root,name))
            last,reason=min(stops) if stops else (solver.t,None)
            ts.append(last);ys.append(dense(last));polys.append(dense)
            if reason: status=1;message=reason;break
    except (ValueError,ArithmeticError) as exc:
        status=-2;message=str(exc)+'; last accepted state retained'
    if not polys: raise RuntimeError(message)
    return SimpleNamespace(t=np.array(ts),y=np.array(ys).T,sol=OdeSolution(ts,polys),
                           status=status,message=message,nfev=solver.nfev,settings=settings,
                           completed=(status==0 and ts[-1]==end))

def solve(amplitude,mix=Mixture(),eta=1.,mode='coupled',settings=Settings(),path=None):
    p=Parameters();coupled=mode=='coupled'
    if not coupled and path is None: raise ValueError('One-way requires saved mechanics')
    def state(t,y):
        return (y[0],y[1],T0*y[2]) if coupled else (*path(t),T0*y[0])
    def fun(t,y):
        R,U,T=state(t,y)
        if R<=0 or T<=0 or not np.all(np.isfinite(y)): raise InvalidState('Nonpositive/nonfinite trial state')
        try: dT=mix.temperature_rate(R,U,T,eta,trial=True)
        except ValueError as exc: raise InvalidState(str(exc)) from exc
        work=mix.pressure(R,T)*4*math.pi*R*R*U/mix.Eref
        loss=mix.heat_loss(R,T,eta)/mix.Eref
        ledger=[work,loss,max(-work,0.)]
        return np.array([U,acceleration(t,R,U,T,amplitude,mix,eta),dT/T0,*ledger] if coupled else [dT/T0,*ledger])
    guards=[]
    if coupled:
        guards.extend([('technical_radius_1pm',lambda t,y:y[0]-1e-12),
                       ('technical_KM_denominator',lambda t,y:1-y[1]/p.c+4*p.mu/(p.rho*p.c*y[0])-1e-6)])
    guards.append(('coefficient_lower_boundary',lambda t,y:state(t,y)[2]-max(mix.Tmin,1e-6)))
    if math.isfinite(mix.Tmax): guards.append(('coefficient_upper_boundary',lambda t,y:mix.Tmax-state(t,y)[2]))
    y0=[p.R0,0.,1.,0.,0.,0.] if coupled else [1.,0.,0.,0.]
    atol=[1e-16,1e-8,2e-10,1e-8,1e-8,1e-8] if coupled else [2e-10,1e-8,1e-8,1e-8]
    end=5/p.f if coupled else path.end
    sol=integrate(fun,y0,end,settings,guards,atol)
    sol.fun,sol.state,sol.mix,sol.eta,sol.mode,sol.amplitude,sol.path=fun,state,mix,eta,mode,amplitude,path
    return sol
