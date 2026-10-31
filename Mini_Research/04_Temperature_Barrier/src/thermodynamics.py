"""Frozen composition, ideal EOS, NASA9 sensible energy with continuous joins.

NASA fit evaluation outside its domain is allowed ONLY for internal solver trial
stages bracketing a terminal event. Public evaluations reject out-of-range T.
No accepted trajectory past a data boundary is reported.
"""
from dataclasses import dataclass
import json
import math
import numpy as np
from .upstream import ROOT, Parameters, RU

T0 = 293.15
DATA = json.loads((ROOT/'results/species_coefficients.json').read_text())

def volume(R):
    return 4*np.pi*np.asarray(R)**3/3

def cp_R(a,T):
    return a[0]/T**2+a[1]/T+a[2]+T*(a[3]+T*(a[4]+T*(a[5]+T*a[6])))

def primitive_cv_R(a,T):
    return -a[0]/T+a[1]*math.log(T)+(a[2]-1)*T+a[3]*T**2/2+a[4]*T**3/3+a[5]*T**4/4+a[6]*T**5/5

def species_cp(species,T):
    if species=='Ar_ground': return 2.5*RU
    data=DATA['species'][species]
    if not data['intervals'][0]['Tmin_K']<=T<=data['intervals'][-1]['Tmax_K']:
        raise ValueError('Temperature outside species coefficient range')
    interval=next((x for x in data['intervals'] if T<=x['Tmax_K']),data['intervals'][-1])
    return RU*cp_R(interval['coefficients'][:7],T)

@dataclass(frozen=True)
class Mixture:
    gas: str = 'air'
    cv_mode: str = 'variable'
    vapor_factor: float = 1.
    convention: str = 'mixture'

    def __post_init__(self):
        if self.gas not in ('air','argon') or self.cv_mode not in ('variable','frozen','mr3'):
            raise ValueError('Unknown gas/capacity')
        if self.convention not in ('mixture','legacy') or self.vapor_factor<0:
            raise ValueError('Invalid vapor convention')
        p=Parameters(); peq=p.p0+2*p.sigma/p.R0
        if self.convention=='legacy' and self.vapor_factor!=0:
            raise ValueError('Separate constant vapor pressure cannot coexist with H2O inventory')
        nt=(p.pg0 if self.convention=='legacy' else peq)*float(volume(p.R0))/(RU*T0)
        xv=0. if self.convention=='legacy' else self.vapor_factor*p.pv/peq
        if not 0<=xv<1: raise ValueError('Vapor fraction outside [0,1)')
        fractions=({'N2':.79*(1-xv),'O2':.21*(1-xv)} if self.gas=='air'
                   else {'Ar_ground':1-xv})
        if xv: fractions['H2O']=xv
        kg=.0264 if self.gas=='air' else .0177
        mass={'N2':.0280134,'O2':.0319988,'Ar_ground':.039948,'H2O':.01801528}
        low=max(DATA['species'][s]['intervals'][0]['Tmin_K'] for s in fractions if s!='Ar_ground') if any(s!='Ar_ground' for s in fractions) else 0.
        high=min(DATA['species'][s]['intervals'][-1]['Tmax_K'] for s in fractions if s!='Ar_ground') if any(s!='Ar_ground' for s in fractions) else math.inf
        breaks=sorted({1000.,6000.,high}) if math.isfinite(high) else [math.inf]
        breaks=[x for x in breaks if x<=high]
        intervals=[];left=low
        for right in breaks:
            mid=(left+right)/2 if math.isfinite(right) else T0
            a=np.zeros(7)
            for s,x in fractions.items():
                if s=='Ar_ground': a[2]+=2.5*x
                else:
                    d=next(d for d in DATA['species'][s]['intervals'] if d['Tmin_K']<=mid<=d['Tmax_K'])
                    a+=x*np.array(d['coefficients'][:7])
            intervals.append((left,right,a));left=right
        cv0=sum(x*(species_cp(s,T0)-RU) for s,x in fractions.items())
        if self.cv_mode=='mr3':
            if self.convention!='legacy': raise ValueError('MR3 limit requires legacy vapor convention')
            cv0=RU/(.4 if self.gas=='air' else 2/3)
        if self.cv_mode!='variable': low,high=0.,math.inf
        for key,val in dict(nt=nt,fractions=fractions,kg=kg,mass_mol=sum(fractions[s]*mass[s] for s in fractions),
                            Tmin=low,Tmax=high,intervals=intervals,cv0=cv0).items():
            object.__setattr__(self,key,val)

    @property
    def p_constant(self): return Parameters().pv if self.convention=='legacy' else 0.
    @property
    def Eref(self): return self.nt*RU*T0
    def _a(self,T):
        return next((a for lo,hi,a in self.intervals if T<=hi),self.intervals[-1][2])
    def cv(self,T,trial=False):
        if T<=0 or (not trial and not self.Tmin-1e-7<=T<=self.Tmax+1e-7):
            raise ValueError('Temperature outside closure domain')
        cv=self.cv0 if self.cv_mode!='variable' else RU*(cp_R(self._a(T),T)-1)
        if cv<=0 or not math.isfinite(cv): raise ValueError('Nonpositive heat capacity')
        return cv
    def capacity(self,T,trial=False): return self.nt*self.cv(T,trial)
    def cv_array(self,T):
        """Same polynomial arithmetic for vectorized diagnostics inside the domain."""
        T=np.asarray(T,dtype=float)
        if np.any(T<=0) or np.any(T<self.Tmin-1e-7) or np.any(T>self.Tmax+1e-7):
            raise ValueError('Temperature outside closure domain')
        if self.cv_mode!='variable': return np.full_like(T,self.cv0)
        value=np.empty_like(T);left=-math.inf
        for lo,hi,a in self.intervals:
            mask=(T>left)&(T<=hi+1e-7);value[mask]=RU*(cp_R(a,T[mask])-1);left=hi
        return value
    def du(self,T):
        """Sensible U(T)-U(T0), J; integrate Cv exactly through all joins."""
        if self.cv_mode!='variable': return self.nt*self.cv0*(T-T0)
        self.cv(T)
        a,b=sorted((T0,T));value=0.
        for lo,hi,coef in self.intervals:
            aa,bb=max(a,lo),min(b,hi)
            if bb>aa: value+=RU*(primitive_cv_R(coef,bb)-primitive_cv_R(coef,aa))
        return self.nt*math.copysign(value,T-T0)
    def du_array(self,T):
        T=np.asarray(T,dtype=float);self.cv_array(T)
        if self.cv_mode!='variable': return self.nt*self.cv0*(T-T0)
        low=np.minimum(T0,T);high=np.maximum(T0,T);value=np.zeros_like(T)
        def primitive(a,x):
            return -a[0]/x+a[1]*np.log(x)+(a[2]-1)*x+a[3]*x**2/2+a[4]*x**3/3+a[5]*x**4/4+a[6]*x**5/5
        for lo,hi,a in self.intervals:
            aa=np.maximum(low,lo);bb=np.minimum(high,hi);mask=bb>aa
            value[mask]+=RU*(primitive(a,bb[mask])-primitive(a,aa[mask]))
        return self.nt*value*np.sign(T-T0)
    def pressure(self,R,T): return self.nt*RU*T/float(volume(R))
    def partial_pressures(self,R,T):
        return {s:x*self.pressure(R,T) for s,x in self.fractions.items()}
    def pressure_rate(self,R,U,T,dT):
        pg=self.pressure(R,T)
        return pg*(dT/T-3*U/R)
    def heat_loss(self,R,T,eta):
        if eta<0: raise ValueError('eta must be nonnegative')
        return 4*math.pi*eta*self.kg*R*(T-T0)
    def temperature_rate(self,R,U,T,eta,trial=False):
        return (-self.pressure(R,T)*4*math.pi*R*R*U-self.heat_loss(R,T,eta))/self.capacity(T,trial)
