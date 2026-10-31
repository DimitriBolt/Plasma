"""Dimensional Rayleigh-Plesset model. All input and solver states use SI units."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import DOP853, solve_ivp
from scipy.optimize import brentq


class InvalidTrialRadius(ValueError):
    """An internal trial stage has left the model domain."""


class PositiveDOP853(DOP853):
    """DOP853: reject invalid trial radii, shrink the step, and retry. No clipping."""
    def _step_impl(self):
        for _ in range(40):
            try:
                return super()._step_impl()
            except InvalidTrialRadius:
                self.h_abs *= 0.25
        return False, "Could not find a positive-radius trial step"


@dataclass(frozen=True)
class Parameters:
    R0: float = 6e-6
    f: float = 30e3
    p0: float = 101325.0
    rho: float = 998.2
    mu: float = 0.001002
    sigma: float = 0.07274
    pv: float = 2339.0
    kappa: float = 1.4
    c: float = 1482.0
    floor_fraction: float = 0.001

    def __post_init__(self):
        positive = (self.R0, self.f, self.p0, self.rho, self.mu,
                    self.sigma, self.kappa, self.c)
        if not all(np.isfinite(x) and x > 0 for x in positive):
            raise ValueError("Physical scales must be finite and positive")
        if not np.isfinite(self.pv) or self.pv < 0 or self.pg0 <= 0:
            raise ValueError("Require pv >= 0 and equilibrium gas pressure > 0")
        if not 0 < self.floor_fraction < 1:
            raise ValueError("Require 0 < floor_fraction < 1")

    @property
    def T(self):
        return 1 / self.f

    @property
    def pg0(self):
        return self.p0 + 2 * self.sigma / self.R0 - self.pv

    @property
    def R_floor(self):
        return self.R0 * self.floor_fraction

    @property
    def omega0(self):
        return np.sqrt((3 * self.kappa * self.pg0 - 2 * self.sigma / self.R0)
                       / (self.rho * self.R0**2))


def external_pressure(t, Pa, p):
    if not np.isfinite(Pa) or Pa < 0:
        raise ValueError("Acoustic amplitude must be finite and nonnegative")
    return p.p0 - Pa * np.sin(2 * np.pi * p.f * np.asarray(t))


def rhs(t, y, Pa, p):
    """Return [U, dU/dt] in m/s and m/s^2, without clipping the radius."""
    R, U = y
    if not np.all(np.isfinite(y)) or R <= 0:
        raise InvalidTrialRadius("Nonfinite state or nonpositive trial radius")
    pressure = (p.pg0 * (p.R0 / R)**(3 * p.kappa) + p.pv
                - external_pressure(t, Pa, p) - 2 * p.sigma / R
                - 4 * p.mu * U / R)
    return np.array([U, pressure / (p.rho * R) - 1.5 * U**2 / R])


def integrate(Pa, p=Parameters(), *, cycles=5, method="DOP853", rtol=1e-10,
              atol=(1e-16, 1e-8), steps_per_period=1000, y0=None):
    """Adaptive SI integration with dense output and four event functions.

    The floor is a numerical termination boundary, never a physical minimum.
    Invalid trial states raise instead of silently modifying the equations.
    """
    external_pressure(0, Pa, p)
    if cycles <= 0 or steps_per_period <= 0:
        raise ValueError("Positive duration and step resolution required")
    initial = np.array([p.R0, 0.] if y0 is None else y0, dtype=float)
    if initial.shape != (2,) or not np.all(np.isfinite(initial)) or initial[0] <= p.R_floor:
        raise ValueError("Initial state must lie strictly above the numerical floor")

    def floor(t, y):
        return y[0] - p.R_floor
    floor.terminal, floor.direction = True, -1

    def maximum(t, y):
        return y[1]
    maximum.direction = -1

    def minimum(t, y):
        return y[1]
    minimum.direction = 1

    def velocity_extremum(t, y):
        return rhs(t, y, Pa, p)[1]

    sol = solve_ivp(lambda t, y: rhs(t, y, Pa, p), (0., cycles * p.T), initial,
                     method=PositiveDOP853 if method == "DOP853" else method, rtol=rtol, atol=np.asarray(atol),
                     max_step=p.T / steps_per_period, dense_output=True,
                     events=(floor, maximum, minimum, velocity_extremum))
    # Relocalize events: SciPy's default absolute event-time tolerance (~1e-15 s)
    # is too coarse for the formal sub-picosecond rebound in the 1.3 bar run.
    # Accepted-step dense output is unchanged; brentq uses SI xtol=1e-24 s.
    for index, fun in ((1, maximum), (2, minimum), (3, velocity_extremum)):
        roots = []
        for te in sol.t_events[index]:
            if te <= p.T * 1e-9:
                continue
            j = min(max(np.searchsorted(sol.t, te, side="right")-1, 0), len(sol.t)-2)
            for k in (j, j-1, j+1):
                if not 0 <= k < len(sol.t)-1:
                    continue
                a, b = sol.t[k:k+2]
                fa, fb = fun(a, sol.sol(a)), fun(b, sol.sol(b))
                if fa * fb <= 0:
                    roots.append(brentq(lambda t:fun(t,sol.sol(t)),a,b,
                                        xtol=1e-24,rtol=4*np.finfo(float).eps))
                    break
            else:
                raise RuntimeError("Could not bracket an extremum on accepted steps")
        sol.t_events[index] = np.unique(roots)
        sol.y_events[index] = (sol.sol(sol.t_events[index]).T if roots else np.empty((0, 2)))
    return sol


def event_times(sol, index, p):
    """Exclude the initial equilibrium point from local-extremum counts."""
    return sol.t_events[index][sol.t_events[index] > p.T * 1e-9]


def sample_times(sol, p):
    """Export uniform points plus four subdivisions of every adaptive interval."""
    mesh = sol.t[:-1, None] + np.diff(sol.t)[:, None] * np.arange(4)[None, :] / 4
    return np.unique(np.concatenate((np.linspace(0, sol.t[-1], 5001),
                                     mesh.ravel(), sol.t, *sol.t_events)))


def summarize(sol, Pa, p):
    completed = bool(sol.success and sol.status == 0
                     and np.isclose(sol.t[-1], 5 * p.T, rtol=1e-12, atol=0))
    # Equilibrium roundoff is not counted as a sequence of physical oscillations.
    maxima = event_times(sol, 1, p) if Pa else np.array([])
    minima = event_times(sol, 2, p) if Pa else np.array([])
    velocity_roots = event_times(sol, 3, p) if Pa else np.array([])
    radius_times = np.unique(np.r_[0., sol.t[-1], maxima, minima])
    speed_times = np.unique(np.r_[0., sol.t[-1], velocity_roots])
    radii = sol.sol(radius_times)[0]
    velocities = sol.sol(speed_times)[1]
    Rmax, Rmin = float(radii.max()), float(radii.min())
    Uplus = float(max(0., velocities.max()))
    Uminus = float(max(0., -velocities.min()))
    tc = dtc = dte = tmax = np.nan
    if len(minima):
        tc = float(minima[np.argmin(sol.sol(minima)[0])])
        preceding = maxima[maxima < tc]
        if len(preceding):
            tmax = float(preceding[-1])
            dtc = tc - tmax
            before = minima[minima < tmax]
            dte = tmax - (float(before[-1]) if len(before) else 0.)
    mach = max(Uplus, Uminus) / p.c
    warnings = []
    if Pa > p.p0:
        warnings.append("negative_far_field_absolute_pressure")
    if mach >= 0.1:
        warnings.append("wall_Mach_ge_0.1")
    if mach >= 1:
        warnings.append("wall_Mach_ge_1")
    if sol.status == 1:
        warnings.append("radius_floor_reached_metrics_censored")
    if not completed:
        warnings.append("incomplete_solution")
    return dict(Pa_Pa=Pa, Pa_bar=Pa / 1e5, Rmax_m=Rmax, Rmin_m=Rmin,
                Rmax_over_R0=Rmax / p.R0, Rmin_over_R0=Rmin / p.R0,
                Umax_positive_m_s=Uplus, Umax_inward_m_s=Uminus,
                max_abs_U_m_s=max(Uplus, Uminus), max_wall_Mach=mach,
                collapse_time_s=tc, preceding_maximum_time_s=tmax,
                collapse_duration_s=dtc, preceding_expansion_duration_s=dte,
                local_maxima_count=len(maxima), local_minima_count=len(minima),
                min_external_pressure_Pa=p.p0-Pa, completed=completed,
                solver_status=sol.status, nfev=sol.nfev,
                accepted_steps=len(sol.t)-1, solver_message=sol.message,
                warnings=";".join(warnings) or "none")


def linear_solution(t, Pa, p):
    """Independent closed-form damped linear response with zero initial data."""
    w, w0 = 2 * np.pi * p.f, p.omega0
    beta = 2 * p.mu / (p.rho * p.R0**2)
    wd = np.sqrt(w0**2 - beta**2)
    forcing = Pa / (p.rho * p.R0**2)
    denom = (w0**2 - w**2)**2 + (2 * beta * w)**2
    A = forcing * (w0**2 - w**2) / denom
    B = -forcing * 2 * beta * w / denom
    C, D = -B, (-A*w-beta*B)/wd
    t = np.asarray(t)
    decay = np.exp(-beta*t)
    h = C*np.cos(wd*t)+D*np.sin(wd*t)
    eps = A*np.sin(w*t)+B*np.cos(w*t)+decay*h
    deps = A*w*np.cos(w*t)-B*w*np.sin(w*t)+decay*(
        -beta*h-C*wd*np.sin(wd*t)+D*wd*np.cos(wd*t))
    return np.array([p.R0*(1+eps), p.R0*deps])
