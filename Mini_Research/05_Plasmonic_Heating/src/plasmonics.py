"""SI optical/thermal interface for MR6; no bubble or gas-temperature model.

Optics: isolated homogeneous nonmagnetic sphere, real lossless host, local
room-temperature bulk Au dielectric data. Time convention exp(-i omega t),
n_Au = n + i kappa, outgoing h_n^(1). No extrapolation and no clipping.
Thermal: linear conduction, lumped equilibrated Au, concentric water cells,
constant properties, finite interface conductance, bath at outer_radius_m.
Use only with pulse_duration_s >> electron-phonon relaxation (~1.7 ps) and
Au diffusion time, moderate temperature rises, no phase change/flow.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import re
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import spherical_jn, spherical_yn
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]


def _positive(**values):
    for name, value in values.items():
        if not np.all(np.isfinite(value)) or np.any(np.asarray(value) <= 0):
            raise ValueError(f'{name} must be finite and positive in SI')


def _nonnegative(**values):
    for name, value in values.items():
        if not np.all(np.isfinite(value)) or np.any(np.asarray(value) < 0):
            raise ValueError(f'{name} must be finite and nonnegative in SI')


class GoldData:
    """Read archived Johnson.yml; interpolate n and kappa separately vs lambda0.

    wavelength_m: vacuum wavelength in m, scalar or array. 'pchip' (default)
    preserves local monotonicity; 'linear' is a sensitivity control. Returns
    dimensionless complex index or relative epsilon=(n+i*kappa)^2. Neither
    method extrapolates. Film data do not include nanoscale damping changes.
    """
    def __init__(self, path=ROOT / 'data/Johnson.yml'):
        rows = []
        for line in Path(path).read_text().splitlines():
            if re.fullmatch(r'\s*[0-9.]+\s+[0-9.]+\s+[0-9.]+\s*', line):
                rows.append([float(v) for v in line.split()])
        self.table = np.array(rows)
        if len(rows) != 49 or np.any(np.diff(self.table[:, 0]) <= 0):
            raise ValueError('Expected the 49 ordered Johnson-Christy Au rows')
        self.wavelength_m = self.table[:, 0] * 1e-6
        self.interpolators = [PchipInterpolator(self.wavelength_m, self.table[:, j],
                                               extrapolate=False) for j in (1, 2)]

    def index(self, wavelength_m, method='pchip'):
        w = np.asarray(wavelength_m, dtype=float)
        _positive(wavelength_m=w)
        if np.any(w < self.wavelength_m[0]) or np.any(w > self.wavelength_m[-1]):
            raise ValueError('Wavelength outside measured 187.9--1937 nm interval')
        if method == 'pchip':
            nr, ki = [f(w) for f in self.interpolators]
        elif method == 'linear':
            nr, ki = [np.interp(w, self.wavelength_m, self.table[:, j]) for j in (1, 2)]
        else:
            raise ValueError('Interpolation must be pchip or linear')
        return nr + 1j * ki

    def epsilon(self, wavelength_m, method='pchip'):
        return self.index(wavelength_m, method)**2


@dataclass(frozen=True)
class CrossSections:
    """All cross sections in m^2; x is the dimensionless host size parameter."""
    absorption_m2: object
    scattering_m2: object
    extinction_m2: object
    size_parameter: object


def dipole_cross_sections(radius_m, wavelength_m, n_medium, epsilon_au):
    """Radiation-corrected quasistatic electric dipole (no dynamic depolarization).

    alpha0=4*pi*a^3*(eps-eps_m)/(eps+2 eps_m) [m^3],
    p=epsilon0*eps_m*alpha*E; alpha=alpha0/(1-i*k_m^3*alpha0/(6*pi)).
    Extinction=k_m*Im(alpha), scattering=k_m^4*|alpha|^2/(6*pi).
    Absorption uses the algebraically equivalent material-loss expression,
    k_m*Im(alpha0)/|1-i*k_m^3*alpha0/(6*pi)|^2, avoiding cancellation.
    Quantitative use requires both x << 1 and |m*x| << 1, not just passivity.
    """
    _positive(radius_m=radius_m, wavelength_m=wavelength_m, n_medium=n_medium)
    eps = np.asarray(epsilon_au, dtype=complex)
    if np.any(eps.imag < 0) or not np.all(np.isfinite(eps)):
        raise ValueError('Passive finite epsilon with Im(epsilon)>=0 required')
    km = 2*np.pi*n_medium/np.asarray(wavelength_m)
    alpha0 = 4*np.pi*radius_m**3*(eps-n_medium**2)/(eps+2*n_medium**2)
    denom = 1-1j*km**3*alpha0/(6*np.pi)
    alpha = alpha0/denom
    sca = km**4*np.abs(alpha)**2/(6*np.pi)
    ext = km*alpha.imag
    ab = km*alpha0.imag/np.abs(denom)**2
    return CrossSections(ab, sca, ext, km*radius_m)


def mie_cross_sections(radius_m, wavelength_m, n_medium, index_au, extra_orders=0):
    """Direct spherical-Bessel Mie series, cross sections in m^2.

    radius_m scalar m; wavelength_m vacuum m and index_au dimensionless
    scalar/array; real n_medium>0. index_au uses n+i*kappa with kappa>=0.
    Default order ceil(max(x)+4*max(x)^(1/3)+2). extra_orders >=0 is a
    convergence control. Audited design x<=2.36, |m*x|<6; no claim of
    general large-x or highly absorbing-particle numerical robustness.
    Absorption is the channel sum of Re(a+b)-|a|^2-|b|^2. Negative values
    beyond roundoff raise rather than being silently set to zero.
    """
    _positive(radius_m=radius_m, wavelength_m=wavelength_m, n_medium=n_medium)
    if np.ndim(radius_m) or np.ndim(n_medium) or isinstance(extra_orders, bool) or int(extra_orders) != extra_orders or extra_orders < 0:
        raise ValueError('Scalar radius/host and integer nonnegative extra_orders required')
    ind = np.asarray(index_au, dtype=complex)
    if np.any(ind.real < 0) or np.any(ind.imag < 0) or not np.all(np.isfinite(ind)):
        raise ValueError('Finite passive n+i*kappa required')
    km = 2*np.pi*n_medium/np.asarray(wavelength_m)
    x, m = np.broadcast_arrays(km*radius_m, ind/n_medium)
    z = m*x
    nmax = int(np.ceil(np.max(x)+4*np.max(x)**(1/3)+2))+extra_orders
    ext, sca, ab = [np.zeros(x.shape) for _ in range(3)]
    for order in range(1, nmax+1):
        jx = spherical_jn(order, x)
        px = x*jx
        dpx = jx+x*spherical_jn(order, x, derivative=True)
        yz = spherical_yn(order, x)
        hx = px+1j*x*yz
        dhx = dpx+1j*(yz+x*spherical_yn(order, x, derivative=True))
        jz = spherical_jn(order, z)
        pz = z*jz
        dpz = jz+z*spherical_jn(order, z, derivative=True)
        an = (m*pz*dpx-px*dpz)/(m*pz*dhx-hx*dpz)
        bn = (pz*dpx-m*px*dpz)/(pz*dhx-m*hx*dpz)
        e = (2*order+1)*(an+bn).real
        s = (2*order+1)*(np.abs(an)**2+np.abs(bn)**2)
        ext += e
        sca += s
        ab += e-s
    factor = 2*np.pi/km**2
    ext, sca, ab = ext*factor, sca*factor, ab*factor
    if np.any(ab < -1e-12*np.maximum(ext, np.finfo(float).tiny)):
        raise ArithmeticError('Mie passivity failed; inspect convention or numerics')
    return CrossSections(ab, sca, ext, x)


def optical_cross_sections(radius_m, wavelength_m, n_medium=1.33, model='mie',
                           interpolation='pchip', data=None):
    """Convenience optical API: a[m], lambda0[m], real host index -> CrossSections."""
    data = GoldData() if data is None else data
    index = data.index(wavelength_m, interpolation)
    if model == 'mie':
        return mie_cross_sections(radius_m, wavelength_m, n_medium, index)
    if model == 'dipole':
        return dipole_cross_sections(radius_m, wavelength_m, n_medium, index**2)
    raise ValueError('model must be mie or dipole')


def absorbed_power(intensity_W_m2, absorption_m2):
    """Local (not incident-before-attenuation) intensity * cross section -> W."""
    _nonnegative(intensity_W_m2=intensity_W_m2, absorption_m2=absorption_m2)
    return np.asarray(intensity_W_m2)*np.asarray(absorption_m2)


def pulse_energy(fluence_J_m2, absorption_m2):
    """Local time-integrated intensity [J/m^2] * sigma_abs [m^2] -> J."""
    return absorbed_power(fluence_J_m2, absorption_m2)


def ensemble_energy(number_density_m3, radius_m, thickness_m, beam_area_m2,
                    fluence_J_m2, cross):
    """Depth-integrated collimated single-pass layer, losses at entry excluded.

    Uniform independent spheres; scattered photons leave the modeled beam.
    Returns J in absorbed/scattered/transmitted channels, dimensionless
    fractions/optical depths, gold fraction, spacing [m], and q(0)/I [m^-1].
    Scattered energy is a first-event channel, NOT necessarily escaped energy
    in a real optically thick scattering sample; no diffuse reabsorption.
    Inputs finite, density/fluence >=0, all geometry positive; cross in m^2.
    """
    _positive(radius_m=radius_m, thickness_m=thickness_m, beam_area_m2=beam_area_m2)
    _nonnegative(number_density_m3=number_density_m3, fluence_J_m2=fluence_J_m2,
                 absorption_m2=cross.absorption_m2, scattering_m2=cross.scattering_m2,
                 extinction_m2=cross.extinction_m2)
    ab, sc, ex = [np.asarray(v) for v in (cross.absorption_m2, cross.scattering_m2, cross.extinction_m2)]
    if not np.allclose(ab+sc, ex, rtol=1e-10, atol=0):
        raise ValueError('Cross sections must conserve extinction=absorption+scattering')
    n = np.asarray(number_density_m3)
    tau = n*ex*thickness_m
    removed = -np.expm1(-tau)
    eta_a = np.divide(ab, ex, out=np.zeros_like(ab, dtype=float), where=ex>0)*removed
    eta_s = np.divide(sc, ex, out=np.zeros_like(sc, dtype=float), where=ex>0)*removed
    incident = np.asarray(fluence_J_m2)*beam_area_m2
    with np.errstate(divide='ignore'):
        spacing = np.where(n>0, n**(-1/3), np.inf)
    return dict(incident_J=incident, absorbed_J=incident*eta_a,
                scattered_J=incident*eta_s, transmitted_J=incident*np.exp(-tau),
                absorption_fraction=eta_a, scattering_fraction=eta_s,
                transmission_fraction=np.exp(-tau), tau_ext=tau,
                tau_sca=n*sc*thickness_m, volume_fraction=n*4*np.pi*radius_m**3/3,
                spacing_m=spacing, q_per_intensity_m1=n*ab)


@dataclass(frozen=True)
class ThermalProperties:
    """Baffou & Rigneault (2011), Table I; near-room-temperature constants in SI."""
    rho_gold_kg_m3: float = 19320.
    cp_gold_J_kgK: float = 129.
    k_gold_W_mK: float = 317.
    rho_water_kg_m3: float = 1000.
    cp_water_J_kgK: float = 4187.
    k_water_W_mK: float = 0.60
    interface_W_m2K: float = 1e8
    ambient_K: float = 293.15
    electron_phonon_s: float = 1.7e-12

    def __post_init__(self):
        _positive(**asdict(self))


def thermal_scales(radius_m, props=ThermalProperties()):
    """Capacity [J/K], diffusion/RC times [s], resistance [K/W], Kapitza ratio."""
    _positive(radius_m=radius_m)
    a, p = radius_m, props
    cap = 4*np.pi*a**3/3*p.rho_gold_kg_m3*p.cp_gold_J_kgK
    rw = 1/(4*np.pi*p.k_water_W_mK*a)
    ri = 1/(4*np.pi*a*a*p.interface_W_m2K)
    return dict(capacity_J_K=cap,
                gold_diffusion_s=a*a*p.rho_gold_kg_m3*p.cp_gold_J_kgK/p.k_gold_W_mK,
                water_diffusion_s=a*a*p.rho_water_kg_m3*p.cp_water_J_kgK/p.k_water_W_mK,
                interface_relaxation_s=cap*ri, water_RC_s=cap*rw,
                water_resistance_K_W=rw, interface_resistance_K_W=ri,
                kapitza_ratio=p.k_water_W_mK/(p.interface_W_m2K*a))


def steady_temperature(radius_m, power_W, r_m=None, props=ThermalProperties()):
    """Infinite-water steady rise [K]. r_m>=a -> water; separate Au mean node.

    No convection, phase change, nearby boundary/particle, or bulk accumulation;
    bath fixed at infinity. Interface jump included in Au, excluded from water.
    Uniform Au approximation omits small internal O(k_water/k_gold) gradients.
    """
    _positive(radius_m=radius_m)
    _nonnegative(power_W=power_W)
    if r_m is None:
        r_m = radius_m
    r = np.asarray(r_m)
    if not np.all(np.isfinite(r)) or np.any(r < radius_m):
        raise ValueError('Water profile requires finite r>=a')
    s = thermal_scales(radius_m, props)
    return dict(water_rise_K=power_W/(4*np.pi*props.k_water_W_mK*r),
                gold_rise_K=power_W*(s['water_resistance_K_W']+s['interface_resistance_K_W']))


class RadialHeatModel:
    """Energy-conservative finite-volume water + one uniform gold node.

    Geometric radial cells, exact spherical steady resistances; symmetric
    capacity-weighted eigenproblem supplies exact time propagation for a
    top-hat pulse (within spatial discretization). Domain bath is explicit.
    Radiation, pressure waves, electron temperature, phase changes absent.
    """
    def __init__(self, radius_m, props=ThermalProperties(), cells=180, outer_radius_m=3e-6):
        _positive(radius_m=radius_m, outer_radius_m=outer_radius_m)
        if cells < 10 or int(cells) != cells or outer_radius_m <= radius_m:
            raise ValueError('At least 10 cells and outer radius > a required')
        self.a, self.props, self.outer = radius_m, props, outer_radius_m
        edges = np.geomspace(radius_m, outer_radius_m, cells+1)
        self.r = np.sqrt(edges[:-1]*edges[1:])
        cgold = thermal_scales(radius_m, props)['capacity_J_K']
        cwater = props.rho_water_kg_m3*props.cp_water_J_kgK*4*np.pi/3*np.diff(edges**3)
        self.capacity = np.r_[cgold, cwater]
        self.sqrtc = np.sqrt(self.capacity)
        p = props
        first_res = 1/(4*np.pi*radius_m**2*p.interface_W_m2K)+(1/radius_m-1/self.r[0])/(4*np.pi*p.k_water_W_mK)
        links = np.r_[1/first_res, 4*np.pi*p.k_water_W_mK/(1/self.r[:-1]-1/self.r[1:])]
        bath_g = 4*np.pi*p.k_water_W_mK/(1/self.r[-1]-1/outer_radius_m)
        k = np.zeros((cells+1, cells+1))
        for i, g in enumerate(links):
            k[i,i] += g; k[i+1,i+1] += g
            k[i,i+1] = k[i+1,i] = -g
        k[-1,-1] += bath_g
        self.conductance = k
        symmetric = k/self.sqrtc[:,None]/self.sqrtc[None,:]
        self.rates, self.modes = eigh_tridiagonal(np.diag(symmetric), np.diag(symmetric,1))
        if np.min(self.rates) <= 0:
            raise ArithmeticError('Heat operator must be strictly dissipative')
        self.weights = self.modes[0]/self.sqrtc[0]/self.rates

    def temperatures(self, time_s, power_W, duration_s):
        """Return rises [K], shape (gold+cells, times). Laser on for 0<=t<=duration."""
        _nonnegative(time_s=time_s, power_W=power_W)
        _positive(duration_s=duration_s)
        t = np.atleast_1d(time_s)
        on = np.minimum(t, duration_s)
        cool = np.maximum(t-duration_s, 0)
        modal = -np.expm1(-self.rates[:,None]*on)*np.exp(-self.rates[:,None]*cool)
        return (self.modes@(self.weights[:,None]*modal))*power_W/self.sqrtc[:,None]

    def pulse(self, fluence_J_m2, absorption_m2, duration_s, time_s):
        """Top-hat pulse -> temperatures, storage, net liquid transfer and balance.

        E_liquid is cumulative NET Au-to-water heat = E_input-E_Au,
        including water energy subsequently reaching the domain bath.
        It is not heat contained in a selected local volume or a bubble.
        """
        _nonnegative(fluence_J_m2=fluence_J_m2, absorption_m2=absorption_m2, time_s=time_s)
        _positive(duration_s=duration_s)
        t = np.atleast_1d(time_s)
        total = float(pulse_energy(fluence_J_m2, absorption_m2))
        temps = self.temperatures(t, total/duration_s, duration_s)
        ein = total*np.minimum(t/duration_s, 1)
        egold = self.capacity[0]*temps[0]
        ewater = self.capacity[1:]@temps[1:]
        # Independent integrated bath flux is available through the exact energy
        # antiderivative of the modal solution, not a residual-defined loss.
        on = np.minimum(t, duration_s); off = np.maximum(t-duration_s, 0)
        rates = self.rates[:,None]
        integral = (on-(-np.expm1(-rates*on))/rates
                    +(-np.expm1(-rates*on))*(-np.expm1(-rates*off))/rates)
        integ_last = ((self.modes[-1]*self.weights)[:,None]*integral).sum(axis=0)*total/duration_s/self.sqrtc[-1]
        bath = self.conductance[-1].sum()*integ_last
        return dict(time_s=t, gold_rise_K=temps[0], water_rise_K=temps[1:],
                    energy_input_J=ein, gold_energy_J=egold, water_energy_J=ewater,
                    bath_energy_J=bath, liquid_net_energy_J=ein-egold,
                    balance_residual_J=ein-egold-ewater-bath)

    def peak_water(self, power_W, duration_s, end_s):
        """Refine maximum interface-water rise and its time [s] after top-hat pulse."""
        # Reconstruct liquid temperature immediately outside Au from first-cell
        # temperature and the same series interface/half-cell resistance.
        rcond = (1/self.a-1/self.r[0])/(4*np.pi*self.props.k_water_W_mK)
        rint = 1/(4*np.pi*self.a**2*self.props.interface_W_m2K)
        surface_modes = (self.modes[1]/self.sqrtc[1]
                         +(self.modes[0]/self.sqrtc[0]-self.modes[1]/self.sqrtc[1])*rcond/(rcond+rint))*self.weights
        def surface(t):
            return power_W*np.sum(surface_modes*(-np.expm1(-self.rates*min(t,duration_s)))
                                  *np.exp(-self.rates*max(t-duration_s,0)))
        grid = np.unique(np.r_[np.linspace(0,duration_s,60), duration_s+np.geomspace(duration_s*1e-7,end_s-duration_s,200)])
        values = np.array([surface(t) for t in grid])
        j = np.argmax(values)
        if j == 0 or j == len(grid)-1:
            return float(values[j]), float(grid[j])
        # Optimize normalized time to avoid scipy's default absolute x tolerance.
        fit = minimize_scalar(lambda u:-surface(u*duration_s),
                              bounds=(grid[j-1]/duration_s,grid[j+1]/duration_s),
                              method='bounded', options={'xatol':1e-10})
        return -float(fit.fun), float(fit.x*duration_s)


def ideal_energy_transfer_upper_bound(absorbed_energy_J, selected_volume_m3,
                                     props=ThermalProperties()):
    """Ideal energy-transfer upper-bound scenario; liquid sensible-equivalent K.

    All specified absorbed energy placed in the chosen liquid volume, with no
    losses. Does not predict deposition geometry, liquid phase, gas temperature,
    bubble coupling efficiency or achieved temperature. Volume is an assumption.
    """
    _nonnegative(absorbed_energy_J=absorbed_energy_J)
    _positive(selected_volume_m3=selected_volume_m3)
    return dict(energy_J=absorbed_energy_J, selected_volume_m3=selected_volume_m3,
                sensible_equivalent_rise_K=absorbed_energy_J/(props.rho_water_kg_m3*props.cp_water_J_kgK*selected_volume_m3))
