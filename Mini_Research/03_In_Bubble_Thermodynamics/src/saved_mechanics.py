"""Audit and interpolate published data, without integrating old mechanics."""
import hashlib
import json
import subprocess
import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from .bubble_thermo import ROOT, REPO, MR2, MECHANICAL_COMMIT, Parameters


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_inputs():
    paths = subprocess.check_output(['git', 'ls-files', '-z', 'Mini_Research/01_Bubble_Dynamics',
        'Mini_Research/02_Compressible_Bubble_Dynamics', 'Plasma-plan.tex', 'Library'], cwd=REPO).decode().split('\0')
    inventory = []
    for rel in filter(None, paths):
        file = REPO / rel
        expected = subprocess.check_output(['git', 'show', f'{MECHANICAL_COMMIT}:{rel}'], cwd=REPO)
        actual = sha256(file)
        if hashlib.sha256(expected).hexdigest() != actual:
            raise RuntimeError(f'Input differs from accepted MR2 snapshot: {rel}')
        inventory.append(dict(path=rel, sha256=actual, bytes=file.stat().st_size))
    (ROOT/'results/input_inventory.json').write_text(json.dumps(dict(
        accepted_snapshot=MECHANICAL_COMMIT, all_match=True, files=inventory), indent=2)+'\n')
    admission = pd.read_csv(MR2/'results/admission.csv')
    summary = pd.read_csv(MR2/'results/amplitude_summary.csv')
    rows = []
    for a in admission.to_dict('records'):
        key = f"{a['model']}_{a['Pa_bar']:.3f}bar"
        meta = json.loads((MR2/f'results/trajectories/{key}.json').read_text())
        s = meta['summary']
        with np.load(MR2/f'results/trajectories/{key}_accepted.npz') as z:
            t, R, U = [z[k] for k in ['t_s','R_m','U_m_s']]
        tab = summary[(summary.model == a['model']) & np.isclose(summary.Pa_bar,a['Pa_bar'])].iloc[0]
        assert bool(a['completed']) == s['completed'] == bool(tab.completed)
        assert np.all(np.diff(t)>0) and np.all(R>0) and np.isfinite(U).all()
        assert t[0] == 0 and t[-1] == s['last_time_s']
        assert not s['completed'] or np.isclose(t[-1],5/Parameters().f,rtol=1e-14,atol=0)
        assert s['extrema_scope'] == ('global_on_requested_interval' if s['completed'] else 'observed_interval_censored')
        rows.append({**a, 'extrema_scope':s['extrema_scope'], 'last_time_s':t[-1],
                     'solver_status':s['solver_status'], 'solver_message':s['solver_message'],
                     'accepted_nodes':len(t)})
    pd.DataFrame(rows).to_csv(ROOT/'results/input_admission.csv',index=False)
    return rows


class SavedMechanics:
    """C1 Hermite reconstruction: R at nodes, stored U as slope, U=dR/dt.

    This is a NEW interpolant, not the original integrator dense polynomial.
    Only near-roundoff time duplicates are coalesced; event times take priority.
    No missing tail is extrapolated.
    """
    def __init__(self, model, amplitude, stride=1, grid='export'):
        self.model, self.amplitude = model, amplitude
        self.key = f'{model}_{amplitude:.3f}bar'
        base = MR2 / 'results/trajectories' / self.key
        self.meta = json.loads((base.parent / (self.key+'.json')).read_text())
        self.summary = self.meta['summary']
        self.events = self.meta['events']
        # Round-trip decimal parsing is essential for the sub-picosecond RP core:
        # the default fast parser can perturb times by several 1e-18 s.
        df = pd.read_csv(base.parent/(self.key+'.csv.gz'),float_precision='round_trip')
        required = ['t_s','R_m','U_m_s','p_infty_Pa','pg_Pa','wall_Mach']
        assert all(c in df for c in required)
        self.export = df
        if grid == 'accepted':
            with np.load(base.parent/(self.key+'_accepted.npz')) as z:
                t, R, U = [z[k].copy() for k in required[:3]]
        else:
            t, R, U = [df[c].to_numpy()[::stride] for c in required[:3]]
        # Include authoritative radius extrema and endpoints at all resolutions.
        e = [x for x in self.events if x['kind'] in ('minimum','maximum')]
        t = np.r_[t, [x['time_s'] for x in e], df.t_s.iloc[-1]]
        R = np.r_[R, [x['R_m'] for x in e], df.R_m.iloc[-1]]
        U = np.r_[U, [x['U_m_s'] for x in e], df.U_m_s.iloc[-1]]
        order = np.argsort(t, kind='stable'); t,R,U = t[order],R[order],U[order]
        # Keep last (authoritative event) among exactly equal times.
        keep = np.r_[np.diff(t)>0, True]; t,R,U = t[keep],R[keep],U[keep]
        event_times = set(x['time_s'] for x in e)
        ids = [0]
        for i in range(1,len(t)):
            if t[i]-t[ids[-1]] > 8*np.spacing(max(t[i],1/Parameters().f)):
                ids.append(i)
            elif t[i] in event_times or i == len(t)-1:
                ids[-1] = i
        self.t,self.R,self.U = t[ids],R[ids],U[ids]
        self.radius = CubicHermiteSpline(self.t,self.R,self.U,extrapolate=False)
        self.velocity = self.radius.derivative()
        if np.any(self.radius((self.t[1:]+self.t[:-1])/2)<=0):
            raise RuntimeError('Nonpositive reconstructed radius')
        self.completed = self.summary['completed']
        self.end = self.t[-1]

    def __call__(self,t):
        return self.radius(t),self.velocity(t)

    def radius_extrema_times(self):
        roots = self.velocity.roots(extrapolate=False)
        return np.unique(np.r_[0.,self.end,roots[(roots>0)&(roots<self.end)]])
