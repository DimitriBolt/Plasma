"""SI Keller--Miksis equation and an unmodified imported RP reference."""
from dataclasses import dataclass
import importlib.util
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import DOP853, Radau, OdeSolution
from scipy.optimize import brentq, minimize_scalar
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
RP_PATH = ROOT.parent / '01_Bubble_Dynamics/src/bubble_dynamics.py'
spec = importlib.util.spec_from_file_location('accepted_rp', RP_PATH)
rp = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = rp
spec.loader.exec_module(rp)
Parameters = rp.Parameters
BASE_COMMIT = 'd59bc0d0119b514152f763e2537a0477c7c63663'


def gas_pressure(R, p):
    return p.pg0 * (p.R0 / np.asarray(R)) ** (3 * p.kappa)


def pressure_jump(t, R, U, Pa, p):
    return (gas_pressure(R, p) + p.pv - rp.external_pressure(t, Pa, p)
            - 2 * p.sigma / R - 4 * p.mu * U / R)


def denominator(R, U, p):
    return R * (1 - U / p.c) + 4 * p.mu / (p.rho * p.c)


def pressure_derivative(t, R, U, A, Pa, p):
    return (-3 * p.kappa * gas_pressure(R, p) * U / R
            + 2 * p.sigma * U / R**2
            + Pa * (2 * np.pi * p.f) * np.cos(2 * np.pi * p.f * t)
            - 4 * p.mu * A / R + 4 * p.mu * U**2 / R**2)


def km_rhs(t, y, Pa, p):
    R, U = y
    if not np.all(np.isfinite(y)) or R <= 0:
        raise rp.InvalidTrialRadius('Nonfinite state or nonpositive trial radius')
    D = denominator(R, U, p)
    if D == 0:
        raise ArithmeticError('Keller--Miksis denominator is exactly zero')
    B = (-3 * p.kappa * gas_pressure(R, p) * U / R
         + 2 * p.sigma * U / R**2
         + Pa * (2 * np.pi * p.f) * np.cos(2 * np.pi * p.f * t)
         + 4 * p.mu * U**2 / R**2)
    N = ((1 + U / p.c) * pressure_jump(t, R, U, Pa, p) / p.rho
         + R * B / (p.rho * p.c) - 1.5 * (1 - U / (3 * p.c)) * U**2)
    return np.array([U, N / D])


def implicit_residual(t, R, U, A, Pa, p):
    left = (1 - U / p.c) * R * A + 1.5 * (1 - U / (3 * p.c)) * U**2
    right = ((1 + U / p.c) * pressure_jump(t, R, U, Pa, p) / p.rho
             + R * pressure_derivative(t, R, U, A, Pa, p) / (p.rho * p.c))
    return left - right


@dataclass(frozen=True)
class Settings:
    method: str = 'DOP853'
    rtol: float = 1e-10
    atol_R: float = 1e-16
    atol_U: float = 1e-8
    steps_per_period: int = 1000
    cycles: float = 5.
    # Numerical guard, NOT a Mach applicability threshold. No clipping is used.
    denominator_guard: float = 1e-6
    km_radius_guard_m: float = 1e-12


def solve(Pa, model='KM', p=Parameters(), settings=Settings(), y0=None):
    """Low-level adaptive driver retains the last accepted state on failure.

    RP imports both its RHS and positive-stage DOP853 class. Its equation and
    6 nm technical floor remain unchanged. KM has no Mach-based termination.
    Dense polynomials from every accepted step are retained in memory.
    """
    rp.external_pressure(0, Pa, p)
    if model not in ('RP', 'KM') or settings.method not in ('DOP853', 'Radau'):
        raise ValueError('Unknown model or method')
    if settings.cycles <= 0 or settings.steps_per_period <= 0:
        raise ValueError('Require positive duration and step resolution')
    initial = np.array([p.R0, 0.] if y0 is None else y0, dtype=float)
    floor = p.R_floor if model == 'RP' else settings.km_radius_guard_m
    if initial.shape != (2,) or initial[0] <= floor or not np.all(np.isfinite(initial)):
        raise ValueError('Initial state outside numerical domain')
    if model == 'KM' and denominator(*initial, p) / initial[0] <= settings.denominator_guard:
        raise ValueError('Initial state on or beyond denominator guard')
    fun = lambda t, y: (rp.rhs if model == 'RP' else km_rhs)(t, y, Pa, p)
    cls = rp.PositiveDOP853 if settings.method == 'DOP853' else Radau
    end = settings.cycles * p.T
    solver = cls(fun, 0., initial, end, rtol=settings.rtol,
                 atol=[settings.atol_R, settings.atol_U],
                 max_step=p.T / settings.steps_per_period)
    ts, ys, pieces = [0.], [initial], []
    status, message = 0, 'Reached requested end time'
    try:
        while solver.status == 'running':
            msg = solver.step()
            if solver.status == 'failed':
                status, message = -1, str(msg)
                break
            dense = solver.dense_output()
            # Check accepted endpoint, locate a numerical-domain crossing.
            # The RHS denominator, pressure and velocity are never modified.
            guards = [(lambda t: dense(t)[0] - floor, 'radius_guard_reached')]
            if model == 'KM':
                guards.append((lambda t: denominator(*dense(t), p) / dense(t)[0]
                               - settings.denominator_guard,
                               'denominator_guard_reached'))
            stops = []
            for guard, name in guards:
                if guard(solver.t) <= 0:
                    stop = brentq(guard, ts[-1], solver.t, xtol=1e-24,
                                  rtol=4 * np.finfo(float).eps)
                    stops.append((stop, name))
            last, reason = min(stops) if stops else (solver.t, None)
            pieces.append(dense)
            ts.append(last)
            ys.append(dense(last))
            if reason:
                status, message = 1, reason
                break
    except (ValueError, ArithmeticError, FloatingPointError) as exc:
        status, message = -2, f'{type(exc).__name__}: {exc}; last accepted point retained'
    if not pieces:
        raise RuntimeError('No accepted step: ' + message)
    sol = SimpleNamespace(t=np.array(ts), y=np.array(ys).T,
                          sol=OdeSolution(ts, pieces), status=status,
                          success=status >= 0, message=message, nfev=solver.nfev,
                          njev=solver.njev, nlu=solver.nlu, model=model, Pa=Pa,
                          parameters=p, settings=settings)
    sol.events = find_events(sol, fun) if Pa else []
    return sol


def find_events(sol, fun):
    """Search accepted-step halves, then refine roots on unchanged dense output.

    Radius extrema use U=0; velocity extrema use acceleration=0. This independent
    search does not inherit solve_ivp's ~1e-15 s event tolerance.
    """
    t = np.unique(np.r_[sol.t, (sol.t[:-1] + sol.t[1:]) / 2])
    y = sol.sol(t)
    a = np.array([fun(x, z)[1] for x, z in zip(t, y.T)])
    events = []
    for values, component in ((y[1], 'radius'), (a, 'velocity')):
        ids = np.where(values[:-1] * values[1:] < 0)[0]
        for j in ids:
            f = ((lambda x: sol.sol(x)[1]) if component == 'radius'
                 else (lambda x: fun(x, sol.sol(x))[1]))
            root = brentq(f, t[j], t[j + 1], xtol=1e-24, rtol=4*np.finfo(float).eps)
            if root <= sol.parameters.T * 1e-9:
                continue
            R, U = sol.sol(root)
            kind = (('maximum' if values[j] > 0 else 'minimum') if component == 'radius'
                    else ('velocity_maximum' if values[j] > 0 else 'velocity_minimum'))
            events.append(dict(kind=kind, time_s=float(root), R_m=float(R), U_m_s=float(U)))
    return sorted(events, key=lambda e: e['time_s'])


def sample_times(sol):
    """Uniform grid + adaptive subdivisions + event and local log-scale grids."""
    t = sol.t
    mesh = t[:-1, None] + np.diff(t)[:, None] * np.arange(4)[None, :] / 4
    extra = []
    for e in sol.events:
        te = e['time_s']
        extra.append(te)
        if e['kind'] == 'minimum':
            j = min(max(np.searchsorted(t, te) - 1, 0), len(t) - 2)
            extra.extend(np.linspace(t[j], t[j+1], 65))
    # Explicit successive scales around the deepest local collapse, without
    # assuming that picoseconds or nanoseconds are physically required.
    mins = [e for e in sol.events if e['kind'] == 'minimum']
    if mins:
        tc = min(mins, key=lambda e: e['R_m'])['time_s']
        offsets = np.geomspace(max(np.spacing(tc)*4, 1e-17), 2e-6, 900)
        extra.extend(tc + np.r_[-offsets, 0., offsets])
    out = np.unique(np.r_[np.linspace(t[0], t[-1], 5001), mesh.ravel(), t, extra])
    return out[(out >= t[0]) & (out <= t[-1])]


def summarize(sol):
    p, Pa = sol.parameters, sol.Pa
    maxima = [e for e in sol.events if e['kind'] == 'maximum']
    minima = [e for e in sol.events if e['kind'] == 'minimum']
    radius_events = maxima + minima
    rtimes = np.r_[0., sol.t[-1], [e['time_s'] for e in radius_events]]
    utimes = np.r_[0., sol.t[-1], [e['time_s'] for e in sol.events if 'velocity_' in e['kind']]]
    R = sol.sol(rtimes)[0]
    U = sol.sol(utimes)[1]
    Rmax, Rmin = float(R.max()), float(R.min())
    plus, inward = max(0., float(U.max())), max(0., -float(U.min()))
    deepest = min(minima, key=lambda e: e['R_m']) if minima else None
    # First strong collapse: operationally the first minimum below R0/2.
    # This is an event selection rule, not a physical validity cutoff.
    strong = next((e for e in minima if e['R_m'] < p.R0 / 2), None)
    def pair(e):
        prev = [x for x in maxima if e and x['time_s'] < e['time_s']]
        return prev[-1] if prev else None
    prev, firstprev = pair(deepest), pair(strong)
    sample = sample_times(sol)
    rr, uu = sol.sol(sample)
    d = denominator(rr, uu, p)
    # Bound min D/R through velocity events and R events as well as dense export;
    # additionally minimize over the local bracket of the sampled minimum.
    q = d / rr
    j = int(np.argmin(q))
    left, right = sample[max(0,j-1)], sample[min(len(sample)-1,j+1)]
    if right > left:
        z = minimize_scalar(lambda s: denominator(*sol.sol(left+(right-left)*s),p)
                            /sol.sol(left+(right-left)*s)[0],bounds=(0,1),method='bounded',
                            options={'xatol':1e-12})
        qmin = min(float(q[j]), float(z.fun))
    else:
        qmin = float(q[j])
    pL = gas_pressure(rr, p) + p.pv - 2*p.sigma/rr - 4*p.mu*uu/rr
    chi = np.abs(pL - p.p0) / (p.rho * p.c**2)
    completed = sol.status == 0 and sol.t[-1] == sol.settings.cycles * p.T
    ratio = Rmax / Rmin
    # Ratios of global extrema can span different cycles. Also export paired
    # extrema for the characteristic collapse and every individual cycle.
    paired_ratio = prev['R_m']/deepest['R_m'] if prev and deepest else np.nan
    result = dict(model=sol.model, Pa_bar=Pa/1e5, Rmax_m=Rmax, Rmin_m=Rmin,
                  Umax_positive_m_s=plus, Umax_inward_m_s=inward,
                  max_abs_U_m_s=max(plus,inward), max_wall_Mach=max(plus,inward)/p.c,
                  collapse_time_s=deepest['time_s'] if deepest else np.nan,
                  preceding_maximum_time_s=prev['time_s'] if prev else np.nan,
                  collapse_duration_s=deepest['time_s']-prev['time_s'] if prev else np.nan,
                  characteristic_minimum_cycle=int(deepest['time_s']/p.T)+1 if deepest else 0,
                  Rmax_over_Rmin=ratio, Vmax_over_Vmin=ratio**3,
                  paired_Rmax_over_Rmin=paired_ratio, paired_Vmax_over_Vmin=paired_ratio**3,
                  first_strong_time_s=strong['time_s'] if strong else np.nan,
                  first_strong_Rmin_m=strong['R_m'] if strong else np.nan,
                  first_strong_duration_s=strong['time_s']-firstprev['time_s'] if firstprev else np.nan,
                  local_maxima_count=len(maxima), local_minima_count=len(minima),
                  min_external_pressure_Pa=p.p0-Pa, max_gas_pressure_Pa=float(gas_pressure(Rmin,p)),
                  max_liquid_pressure_change_over_bulk_modulus=float(chi.max()),
                  min_D_over_R=qmin, min_D_m=float(d.min()),
                  min_export_step_s=float(np.diff(sample).min()),
                  last_time_s=float(sol.t[-1]), last_R_m=float(sol.y[0,-1]),
                  last_U_m_s=float(sol.y[1,-1]), completed=completed,
                  extrema_scope='global_on_requested_interval' if completed else 'observed_interval_censored',
                  solver_status=sol.status, solver_message=sol.message,
                  accepted_steps=len(sol.t)-1, nfev=sol.nfev)
    return result
