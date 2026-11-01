"""Coupled KM, NASA9 fixed gas, moving conservative water and dynamic Au.

Scaled solver time is microseconds. Physical APIs/archives are SI. Tangent state
is a complex-step directional derivative, never a subtraction of close Tmax.
Only trial stages can bracket coefficient/contact events. Accepted tails stop.
"""
from dataclasses import dataclass, asdict
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
from .upstream import thermo,Parameters,T0,RU
from .transport import MaterialGrid,center_radius
from .particle import Particle

TS=1e-6

@dataclass(frozen=True)
class Config:
    amplitude_bar: float = 1.
    gas: str = 'air'
    vapor_factor: float = 1.
    cells: int = 64
    padding_m: float = 25e-6
    k_water: float = .6
    eta: float = 1.
    cycles: float = 5.
    rtol: float = 2e-9
    atol: float = 2e-11
    max_step_s: float = 1e-7
    method: str = 'BDF'
    clamp_water: bool = False
    formal_material: bool = False
    incident_atol: float | None = None
    grid_faces_m3: tuple | None = None

class Model:
    def __init__(self,cfg=Config(),particle=Particle()):
        if cfg.eta<0 or cfg.cycles<=0 or cfg.rtol<=0 or cfg.atol<=0 or cfg.max_step_s<=0:
            raise ValueError('Nonnegative gas conductance and positive numerical scales required')
        self.cfg,self.particle=cfg,particle
        self.p=Parameters();self.mix=thermo.Mixture(gas=cfg.gas,vapor_factor=cfg.vapor_factor)
        self.grid=MaterialGrid(cfg.cells,cfg.padding_m,cfg.k_water,faces_m3=cfg.grid_faces_m3)
        self.cross=particle.optical()
        self.abs_per_fluence=float(optics_energy(particle,self.cross))
        self.au=3+cfg.cells;self.ledger=self.au+1;self.dim=self.ledger+6

    def initial(self):
        y=np.zeros(self.dim);y[0]=y[2]=1.
        return y

    def state(self,t,y,path=None):
        R,U=(y[0]*self.p.R0,y[1]*self.p.c) if path is None else path(t)
        Tg=y[2]*T0;tl=y[3:self.au]*T0;ta=y[self.au]*T0
        return R,U,Tg,tl,ta

    def capacity(self,T):
        # Same accepted NASA9 polynomial; complex extension is derivative-only.
        if np.iscomplexobj(T):
            return self.mix.nt*RU*(thermo.cp_R(self.mix._a(float(T.real)),T)-1)
        return self.mix.capacity(float(T),trial=True)

    def powers(self,t,y,absorbed_W=0.,incident_W=0.,path=None):
        R,U,Tg,tl,ta=self.state(t,y,path)
        H=4*np.pi*self.cfg.eta*self.mix.kg*R
        ki,_,_=self.grid.conductances(R)
        H_eff=H if self.cfg.clamp_water else H*ki/(H+ki)
        qg=H_eff*((0. if self.cfg.clamp_water else tl[0])-(Tg-T0))
        rnp=center_radius(R,self.particle.radius_m,self.particle.gap_m,self.particle.motion)
        w=self.grid.surface_weights(R,rnp,self.particle.radius_m)
        qa=self.particle.count*4*np.pi*self.particle.radius_m**2*self.particle.G_W_m2K*(ta-w@tl) if self.particle.count else 0.*R
        pg=self.mix.nt*RU*Tg/(4*np.pi*R**3/3)
        work=pg*4*np.pi*R**2*U
        water,bath=self.grid.water_rate(R,tl,qg,qa,w)
        return qg,qa,bath,work,water,w,rnp,pg

    def rhs(self,t_us,y,absorbed_W=0.,incident_W=0.,path=None):
        t=t_us*TS;p=self.p
        R,U,Tg,tl,ta=self.state(t,y,path)
        if np.real(R)<=0 or np.real(Tg)<=0: raise ValueError('Nonpositive trial state')
        qg,qa,qb,work,water,w,rnp,pg=self.powers(t,y,absorbed_W,incident_W,path)
        dT=(-work+qg)/self.capacity(Tg)
        dpg=pg*(dT/Tg-3*U/R)
        omega=2*np.pi*p.f;amplitude=self.cfg.amplitude_bar*1e5
        dp=pg-(p.p0-amplitude*np.sin(omega*t))-2*p.sigma/R-4*p.mu*U/R
        B=dpg+2*p.sigma*U/R**2+amplitude*omega*np.cos(omega*t)+4*p.mu*U**2/R**2
        D=R*(1-U/p.c)+4*p.mu/(p.rho*p.c)
        accel=((1+U/p.c)*dp/p.rho+R*B/(p.rho*p.c)-1.5*(1-U/(3*p.c))*U**2)/D
        out=np.zeros_like(y)
        if path is None: out[:2]=[U/p.R0,accel/p.c]
        out[2]=dT/T0
        if not self.cfg.clamp_water:
            out[3:self.au]=water/T0
            out[self.au]=(absorbed_W-qa)/(self.particle.count*self.particle.capacity*T0) if self.particle.count else 0.
        out[self.ledger:]=np.array([work,qg,qa,qb,absorbed_W,incident_W])/self.mix.Eref
        return out*TS

    def guards(self,path=None,tangent=False):
        def make(name,fun,direction=-1):
            def event(t,y):
                return float(fun(t*TS,y[:self.dim]))
            event.terminal=True;event.direction=direction;event.name=name
            return event
        s=lambda t,y:self.state(t,y,path)
        ev=[make('gas_coefficient_lower_boundary',lambda t,y:s(t,y)[2]-max(self.mix.Tmin,1e-6)),
            make('technical_radius',lambda t,y:s(t,y)[0]-1e-12),
            make('technical_KM_denominator',lambda t,y:1-s(t,y)[1]/self.p.c+4*self.p.mu/(self.p.rho*self.p.c*s(t,y)[0])-1e-6)]
        if np.isfinite(self.mix.Tmax):
            ev.append(make('gas_coefficient_upper_boundary',lambda t,y:self.mix.Tmax-s(t,y)[2]))
        if self.particle.count:
            def gap(t,y):
                R=s(t,y)[0]
                return center_radius(R,self.particle.radius_m,self.particle.gap_m,self.particle.motion)-R-self.particle.radius_m
            ev.append(make('particle_contact',gap))
            ev.append(make('particle_outside_material_domain',lambda t,y:self.grid.radii(s(t,y)[0])[1][-1]-center_radius(s(t,y)[0],self.particle.radius_m,self.particle.gap_m,self.particle.motion)-self.particle.radius_m))
        if not self.cfg.formal_material and not tangent:
            ev.extend([make('water_reference_373K',lambda t,y:373.15-(T0+np.max(s(t,y)[3]))),
                       make('Au_moderate_rise_30K',lambda t,y:30-s(t,y)[4])])
        return ev

    def solve(self,fluence=0.,center_s=0.,duration_s=10e-9,path=None,tangent=False,end_s=None):
        if fluence<0 or duration_s<=0: raise ValueError('Invalid pulse')
        end=self.cfg.cycles/self.p.f if end_s is None else end_s
        onset=center_s-duration_s/2;offset=center_s+duration_s/2
        if fluence and (onset<0 or offset>end): raise ValueError('Pulse must lie in requested interval')
        y0=self.initial(); y0=np.r_[y0,np.zeros_like(y0)] if tangent else y0
        edges=sorted(set([0.,end]+[x for x in (onset,offset) if 0<x<end and (fluence or tangent)]))
        segments=[];status='completed';nfev=0
        Eunit=self.abs_per_fluence
        Iunit=np.pi*self.particle.beam_radius_m**2
        guards=self.guards(path,tangent)
        for a,b in zip(edges[:-1],edges[1:]):
            on=onset < (a+b)/2 < offset
            P=Eunit/duration_s if on else 0.; incident=Iunit/duration_s if on else 0.
            def fun(t,y):
                if not tangent:return self.rhs(t,y,P*fluence,incident*fluence,path)
                base=y[:self.dim];dy=y[self.dim:]
                f=self.rhs(t,base,0.,0.,path)
                # All min/max choices use real geometry, preserving the local branch.
                direction=self.rhs(t,base.astype(complex)+1j*1e-20*dy,0.,0.,path).imag/1e-20
                source=self.rhs(t,base,P,incident,path)-f
                return np.r_[f,direction+source]
            atol=np.full(len(y0),self.cfg.atol)
            # The independently integrated incident beam energy can exceed the
            # tiny thermal perturbation by many orders of magnitude. It has no
            # feedback into the physics; give its scaled ledger its own units.
            if self.cfg.incident_atol is not None:
                atol[self.ledger+5]=self.cfg.incident_atol
            # Water/gold perturbation states are zero based, with enough absolute precision.
            if tangent:atol[self.dim:]=min(self.cfg.atol,1e-15)
            def jac(t,y):
                columns=[]
                for j in range(len(y)):
                    if not tangent:
                        step=np.zeros(len(y),complex);step[j]=1j*1e-24
                        columns.append(fun(t,y.astype(complex)+step).imag/1e-24)
                    else:
                        h=2e-6*max(abs(y[j]),1e-5)
                        step=np.zeros(len(y));step[j]=h
                        columns.append((fun(t,y+step)-fun(t,y-step))/(2*h))
                return np.array(columns).T
            try:
                sol=solve_ivp(fun,(a/TS,b/TS),y0,method=self.cfg.method,dense_output=True,
                    rtol=self.cfg.rtol,atol=atol,max_step=self.cfg.max_step_s/TS,events=guards,jac=jac)
            except (ValueError,FloatingPointError) as exc:
                if not segments: raise
                status='solver_exception: '+str(exc);break
            segments.append(sol);nfev+=sol.nfev;y0=sol.y[:,-1]
            if sol.status!=0:
                status=next((g.name for g,ts in zip(guards,sol.t_events) if len(ts)),sol.message);break
        if not segments: raise RuntimeError('No accepted interval')
        end_actual=segments[-1].t[-1]*TS
        def dense(t):
            values=np.asarray(t);flat=values.ravel();out=np.empty((len(y0),len(flat)))
            if np.any(flat<0) or np.any(flat>end_actual+1e-19): raise ValueError('No extrapolation')
            ends=np.array([x.t[-1]*TS for x in segments])
            ids=np.minimum(np.searchsorted(ends,flat,side='left'),len(segments)-1)
            for j,seg in enumerate(segments):
                mask=ids==j
                if mask.any():out[:,mask]=seg.sol(flat[mask]/TS)
            return out[:,0] if values.ndim==0 else out.reshape((len(y0),)+values.shape)
        return SimpleNamespace(model=self,segments=segments,sol=dense,t=np.unique(np.concatenate([s.t*TS for s in segments])),
            completed=status=='completed' and abs(end_actual-end)<1e-18,status=status,end=end_actual,
            requested_end=end,fluence=fluence,center=center_s,duration=duration_s,path=path,tangent=tangent,nfev=nfev)

    def sensitivity(self,base,center_s,duration_s=10e-9):
        """Linear variational solve along a retained no-laser dense trajectory.

        This avoids re-integrating its large baseline during a weak-source solve.
        The directional complex step differentiates geometry, both thermal fluxes,
        NASA9 capacity and the full mechanical pressure derivative together.
        """
        onset=center_s-duration_s/2;offset=center_s+duration_s/2
        end=base.end
        edges=sorted(set([0.,end]+[x for x in (onset,offset) if 0<x<end]))
        segments=[];y0=np.zeros(self.dim)
        for a,b in zip(edges[:-1],edges[1:]):
            on=onset<(a+b)/2<offset
            P=self.abs_per_fluence/duration_s if on else 0.
            inc=np.pi*self.particle.beam_radius_m**2/duration_s if on else 0.
            def fun(t,dy):
                yy=base.sol(t*TS)[:self.dim]
                tangent=self.rhs(t,yy.astype(complex)+1j*1e-20*dy,path=base.path).imag/1e-20
                src=np.zeros(self.dim)
                if self.particle.count:src[self.au]=P/(self.particle.count*self.particle.capacity*T0)*TS
                src[self.ledger+4:self.ledger+6]=np.array([P,inc])/self.mix.Eref*TS
                return tangent+src
            def jac(t,dy):
                yy=base.sol(t*TS)[:self.dim].astype(complex)
                columns=[]
                for j in range(self.dim):
                    step=np.zeros(self.dim,complex);step[j]=1j*1e-24
                    columns.append(self.rhs(t,yy+step,path=base.path).imag/1e-24)
                return np.array(columns).T
            atols=np.full(self.dim,1e-14);atols[:3]=1e-15
            atols[self.au]=1e-11;atols[self.ledger:]=1e-12;atols[-1]=1e-6
            atols*=min(1.,self.cfg.rtol/2e-8)
            s=solve_ivp(fun,(a/TS,b/TS),y0,method=self.cfg.method,jac=jac,dense_output=True,
                  rtol=self.cfg.rtol,atol=atols,max_step=self.cfg.max_step_s/TS)
            if not s.success:raise RuntimeError('Variational solve failed at t_us='+str(s.t[-1])+', max_state='+str(np.max(abs(s.y[:,-1])))+': '+s.message)
            segments.append(s);y0=s.y[:,-1]
        def dense(t):
            values=np.asarray(t);flat=values.ravel()
            if np.any(flat<0) or np.any(flat>end+1e-19):raise ValueError('No extrapolation')
            out=np.empty((self.dim,len(flat)));ends=np.array([x.t[-1]*TS for x in segments])
            ids=np.minimum(np.searchsorted(ends,flat),len(segments)-1)
            for j,seg in enumerate(segments):
                mask=ids==j
                if mask.any():out[:,mask]=seg.sol(flat[mask]/TS)
            out=out[:,0] if values.ndim==0 else out.reshape((self.dim,)+values.shape)
            return np.concatenate([base.sol(t)[:self.dim],out],axis=0)
        return SimpleNamespace(model=self,segments=segments,sol=dense,
            t=np.unique(np.r_[base.t,np.concatenate([s.t*TS for s in segments])]),
            completed=base.completed,status=base.status,end=base.end,requested_end=base.requested_end,
            fluence=0.,center=center_s,duration=duration_s,path=base.path,tangent=True,
            nfev=sum(s.nfev for s in segments),baseline=base)

def optics_energy(particle,cross):
    from .upstream import optics
    energy=optics.pulse_energy(1.,cross.absorption_m2)
    if not np.isclose(energy,optics.absorbed_power(1e8,cross.absorption_m2)*1e-8,rtol=1e-14,atol=0):
        raise ArithmeticError('MR5 optical power and pulse-energy API disagree')
    return particle.count*particle.absorption_scale*energy
