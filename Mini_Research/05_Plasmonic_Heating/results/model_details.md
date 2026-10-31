# Models, SI interface and interpretation

## 1. Optical conventions

Vacuum wavelength is `wavelength_m`; `n_medium` is real, positive and lossless.
The metal index is `n + i*kappa` with `exp(-i*omega*t)`. Epsilon is relative.
`alpha` is a volume in m^3 with `p = epsilon0*epsilon_medium*alpha*E`.
Mie uses outgoing spherical Hankel functions of the first kind. The independent
miepython comparison conjugates the index to its documented negative-loss sign.

The default PCHIP interpolates measured n and kappa separately, then squares their
complex sum. The original 49-row table is kept unchanged. No extrapolation,
negative-absorption clipping, size-correction fit or temperature correction is
hidden in the calculation. Radiation reaction corrects the quasistatic dipole;
dynamic depolarization and multipoles remain absent from that approximation.

## 2. Callable interface

Import `plasmonics` from the study's `src/` directory, or use
`from src.plasmonics import ...` when this study is the working directory.
The module only depends on NumPy/SciPy and its own archived data.

| API | Inputs | Outputs and assumptions |
|---|---|---|
| `GoldData.index/epsilon` | vacuum wavelength [m], interpolation name | dimensionless n+i*kappa / relative epsilon; rejects outside 187.9–1937 nm |
| `optical_cross_sections` | a [m], vacuum wavelength [m], host index, mie/dipole | `CrossSections` with absorption/scattering/extinction [m^2], size parameter |
| `mie_cross_sections` | a [m], wavelength [m], host index, absolute Au index, extra orders | direct Mie series; audited x<=2.36, not arbitrary large x |
| `dipole_cross_sections` | a [m], wavelength [m], host index, relative Au epsilon | passive radiation-corrected electric dipole; x and internal phase must be small |
| `absorbed_power` | local I [W/m^2], sigma_abs [m^2] | absorbed W; no extra attenuation applied |
| `pulse_energy` | local fluence [J/m^2], sigma_abs [m^2] | absorbed J; not total beam energy |
| `ensemble_energy` | N [m^-3], a,L [m], beam area [m^2], incident-inside-host F [J/m^2], cross sections | J and fractions in three first-pass channels; optical depths, gold fraction, spacing [m], mu_abs [m^-1] |
| `thermal_scales` | a [m], `ThermalProperties` | capacity [J/K], resistances [K/W], diffusion/interface times [s], Kapitza ratio |
| `steady_temperature` | a [m], absorbed power [W], r>=a [m], properties | water profile rise [K], lumped Au rise [K], infinite bath |
| `RadialHeatModel.pulse` | fluence [J/m^2], cross section [m^2], top-hat duration [s], time [s] | temperature rises [K], input/Au/water/bath/net-liquid energies and balance [J] |
| `RadialHeatModel.peak_water` | absorbed power [W], duration/end [s] | refined water-interface maximum [K] and its time [s] |
| `ideal_energy_transfer_upper_bound` | explicitly selected absorbed energy [J], chosen liquid volume [m^3] | energy, volume, water sensible-equivalent rise [K]; no deposition mechanism or gas prediction |

Example for a 20 nm sphere:

```python
from src.plasmonics import optical_cross_sections, RadialHeatModel
c = optical_cross_sections(20e-9, 532e-9, n_medium=1.33)
history = RadialHeatModel(20e-9).pulse(1.0, c.absorption_m2,
                                     10e-9, [10e-9, 20e-9, 110e-9])
```

All temperature results are **rises**, not absolute temperatures. Add the explicit
ambient temperature only when that interpretation is valid. API documentation
does not turn a formal short-pulse/high-temperature evaluation into an admissible
thermal model. Negative and nonfinite physical inputs are rejected.

## 3. Thermal discretization and balances

Water shells have exact spherical volumes. Centers are geometric means of shell
edges. Between node centers, thermal resistance is
`(1/r_i - 1/r_j)/(4*pi*k_water)`. The first link adds interface resistance
`1/(4*pi*a^2*G)` in series. The outer bath is at a stated finite radius.

With capacities C and positive conductance matrix K, `C theta_dot + K theta = P`.
Diagonalizing `C^(-1/2) K C^(-1/2)` gives a positive dissipative spectrum. Modal
top-hat propagation has no time-step integration error. The bath energy is an
analytic integral of the independently propagated outer flux, not defined as a
residual. An independent BDF solver and interface-flux quadrature check it.

Gold storage plus water storage plus bath energy equals absorbed energy up to
roundoff. Net liquid transfer equals absorbed minus gold storage; water storage
alone excludes energy already reaching the bath. None of these integrals implies
energy residing inside an arbitrarily selected local sphere. Before pulse end,
the denominator for a final-pulse-energy fraction differs from energy absorbed
so far; the saved fields keep both meanings distinct.

The reconstructed immediate liquid-surface temperature uses the same half-cell
and interface series resistances; its spatial error is checked against 720 cells.
Gold's peak is at top-hat pulse end because its collocated step-response derivative
is a positive sum of decaying modes. The water peak can occur later. At fixed
energy, the gold response is the average of that decaying impulse kernel, so the
allowable fluence under a fixed gold-rise cap is nondecreasing with pulse duration.
This proves the 100 ns ensemble branch bounds shorter pulses under the same
fixed concentration, geometry and overlap screens.

## 4. Finite design and objective selection

Optical domain: radii 5–100 nm inclusive in 1 nm steps, 400–1000 nm inclusive in
1 nm steps, hosts 1/1.33/1.5. Finite thermal candidates use water and radii
20–100 nm, durations 1/10/100 ns. For each, fluence is varied continuously up to
the minimum of the 30 K temperature allowance, 1000 J/m^2 and the duration times
1e10 W/m^2. Ties in one-particle energy prefer larger absorption cross section.

Ensembles search the 61-value log concentration lattice, three thicknesses and
three beam radii. Both fraction and energy are optimized separately at each
radius/geometry. The result table retains both winners. These are global
**discrete design** optima, not continuous or universal optima. At fixed N,
volumetric heating follows sigma_abs; at fixed gold volume fraction it follows
sigma_abs / gold-particle volume. Fraction optimization is different again.

Primary material/range motivation is Johnson–Christy (film measurements), Jain
(20/40/80 nm diameter spheres; theoretical), Baffou–Rigneault (short-pulse
diffusion/electron relaxation), and Baffou–Quidant (CW/pulse/interface context).
The larger spheres, beam radii, concentrations and engineering limits are explicit
design choices used to probe finite extremes. They are not claimed as one
published experimental protocol.

## 5. Limits needing new physics

Bulk-film epsilon misses size-dependent damping, surface/ligand effects,
polydispersity, nonsphericity and temperature-dependent dielectric response.
Real water dispersion/absorption, entry losses and exact beam propagation are
neglected. Dispersion screens do not exclude rare close pairs in a random sample.
Single-pass scattering does not model diffuse reabsorption or return paths.

Constant-property one-temperature Au and diffusive water are moderate-rise models.
The selected time screen is duration >=10 times max(electron–phonon relaxation,
Au internal diffusion). Shorter pulses need electron/lattice dynamics; the formal
one-temperature histories are diagnostic continuations. The code does not remove
them. Gold melting, reshaping, evaporation, fragmentation, water nucleation,
superheating, latent heats, changing contact and motion require new models.
373.15 K is not imposed as a universal nanobubble threshold. The bulk gold point
1337.33 K is not imposed as a universal nanoparticle melting threshold.

The interface supplies liquid energy only. No bubble radius change, coupled
laser–bubble system, new gas-temperature history or fusion calculation appears.
