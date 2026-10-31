"""Accepted-state integration, event localization and threshold residence."""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
from scipy.integrate import DOP853, Radau, OdeSolution
from scipy.optimize import brentq


class InvalidState(ValueError):
    pass


class DomainDOP853(DOP853):
    def _step_impl(self):
        for _ in range(40):
            try:
                return super()._step_impl()
            except InvalidState:
                self.h_abs *= 0.25
        return False, 'Trial states could not remain in the positive domain'


@dataclass(frozen=True)
class Settings:
    method: str = 'DOP853'
    rtol: float = 2e-9
    steps_per_period: int = 600


def integrate(fun, y0, end, period, atol, settings=Settings(), guards=()):
    cls = DomainDOP853 if settings.method=='DOP853' else Radau
    options={}
    if settings.method=='Radau':
        # A fixed-scale centered Jacobian avoids adaptive finite-difference overflow
        # on the identically zero source ledger and passive integral columns.
        scales=np.maximum(np.abs(y0),1.)
        scales[0]=abs(y0[0])
        def jac(t,y):
            delta=np.sqrt(np.finfo(float).eps)*np.maximum(np.abs(y),scales)
            columns=[]
            for j,h in enumerate(delta):
                step=np.zeros(len(y));step[j]=h
                columns.append((fun(t,y+step)-fun(t,y-step))/(2*h))
            return np.array(columns).T
        options['jac']=jac
    solver = cls(fun,0.,np.asarray(y0,dtype=float),end,rtol=settings.rtol,
                 atol=atol,max_step=period/settings.steps_per_period,**options)
    ts,ys,polys = [0.],[np.asarray(y0,dtype=float)],[]
    status,message=0,'Reached requested end time'
    try:
        while solver.status=='running':
            msg=solver.step()
            if solver.status=='failed':
                status,message=-1,str(msg);break
            dense=solver.dense_output()
            stops=[]
            for guard,name in guards:
                if guard(dense(solver.t))<=0:
                    root=brentq(lambda t:guard(dense(t)),ts[-1],solver.t,xtol=1e-24,
                                rtol=4*np.finfo(float).eps)
                    stops.append((root,name))
            last,why=min(stops) if stops else (solver.t,None)
            ts.append(last);ys.append(dense(last));polys.append(dense)
            if why:
                status,message=1,why;break
    except (ValueError,ArithmeticError,FloatingPointError) as exc:
        status,message=-2,f'{type(exc).__name__}: {exc}; last accepted state retained'
    if not polys:
        raise RuntimeError(message)
    return SimpleNamespace(t=np.array(ts),y=np.array(ys).T,sol=OdeSolution(ts,polys),
         status=status,message=message,completed=(status==0 and ts[-1]==end),nfev=solver.nfev,
         settings=settings)


def roots_from_grid(fun,grid):
    grid=np.unique(grid);v=np.asarray(fun(grid))
    ids=np.flatnonzero(v[:-1]*v[1:]<0)
    roots=[brentq(fun,grid[i],grid[i+1],xtol=1e-24,rtol=4*np.finfo(float).eps) for i in ids]
    roots.extend(grid[v==0])
    return np.unique(roots)


def thermal_extrema(Tfun,dTfun,grid):
    candidates=np.unique(np.r_[grid[0],grid[-1],roots_from_grid(dTfun,grid)])
    values=Tfun(candidates); j=np.argmax(values)
    return float(values[j]),float(candidates[j]),candidates


def hot_intervals(Tfun,grid,threshold):
    roots=roots_from_grid(lambda t:Tfun(t)-threshold,grid)
    edges=np.unique(np.r_[grid[0],roots,grid[-1]])
    mid=(edges[1:]+edges[:-1])/2
    return [(float(a),float(b)) for a,b,hot in zip(edges[:-1],edges[1:],Tfun(mid)>threshold) if hot]
