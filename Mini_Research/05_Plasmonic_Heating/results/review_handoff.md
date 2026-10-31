# Mini-Research 5: independent scientific review handoff

Please perform an independent mathematical, numerical, physical and documentation
audit of Mini-Research 5 before any scientific acceptance. No MR6 work has begun.
Publication does not constitute acceptance or closure by RoadMap.

The dispatch message supplies the full published commit SHA and immutable GitHub,
PDF and package links. In a clone, `git rev-parse HEAD` identifies the containing
snapshot. The starting published main was
`803079ef68c0d3b02514dc88e211892079beeaff`. All 1008 tracked
upstream inputs were verified byte-identical; the existing unrelated user files
were preserved. Only `Mini_Research/05_Plasmonic_Heating/` is added.

Report: `report/05_Plasmonic_Heating.pdf`, nine pages: seven main-text pages,
one verification appendix page and one references page. Author: Dimitri Bolt.
Academic supervisor: Prof. Dr. Ildar R. Gabitov. Scientific text and figures are
English, and the eight figures are JPG only. Code, raw table, tests, pinned
dependencies, parameter tables, CSV/NPZ results, LaTeX/PDF, results notebook,
FILE_LIST.txt and SHA256SUMS form the complete package.

## Mathematical models

- Experimental Johnson–Christy room-temperature film n and k, 49 raw rows,
  archived without alteration. PCHIP of n and k separately against vacuum
  wavelength, then epsilon=(n+i*k)^2. No extrapolation outside 187.9–1937 nm.
- Volume-normalized dipole: p=epsilon0*epsilon_medium*alpha*E,
  alpha0=4*pi*a^3*(epsilon-epsilon_medium)/(epsilon+2*epsilon_medium),
  alpha=alpha0/(1-i*k_medium^3*alpha0/(6*pi)). Absorption uses its nonnegative
  material-loss expression, algebraically equal to extinction minus scattering.
- Direct spherical-Bessel Mie series, positive-imaginary-loss convention and
  outgoing first-kind Hankel functions; independent miepython 3.0.2 comparison
  explicitly conjugates its opposite sign convention. No negative clipping.
- P_abs=I_local*sigma_abs and E_abs=F_local*sigma_abs. The finite single-pass
  layer integrates attenuation: eta_abs=(sigma_abs/sigma_ext)*(1-exp(-N*sigma_ext*L)).
  Absorption, first-event scattering and transmission sum to the incident beam
  energy. Entry-window loss, solvent absorption and diffuse reabsorption excluded.
- Lumped equilibrated Au plus radial finite-volume water conduction, finite
  interface conductance, constant coefficients and an explicit outer bath.
  Exact-in-time top-hat modal propagation; independent interface-flux integral,
  independent BDF solver and spatial/domain controls. Au, water and bath storage
  remain distinct from cumulative net energy transferred to liquid.

## Key parameters and finite domain

| Parameter | Value | Unit / interpretation |
|---|---:|---|
| Optical radius grid | 5–100, step 1 | nm; 96 radii |
| Optical vacuum wavelength grid | 400–1000, step 1 | nm; 601 wavelengths |
| Real optical host index | 1 / 1.33 / 1.5 | only 1.33 combined with water heating |
| Optical grid points | 173088 | all retained |
| Moderate thermal radii | 20–100 | nm |
| Moderate pulse durations | 1 / 10 / 100 | ns, top-hat |
| Maximum fluence | 1000 | J/m^2 = 100 mJ/cm^2 |
| Maximum peak intensity | 1e10 | W/m^2 |
| Selected gold-rise screen | 30 | K above 293.15 K, not a material threshold |
| Au rho / cp / k | 19320 / 129 / 317 | kg/m^3; J/(kg K); W/(m K) |
| Water rho / cp / k | 1000 / 4187 / 0.60 | same SI units |
| Interface G | 50 / 100 / 200 | MW/(m^2 K); moderate upper allowances use 50 |
| Layer N lattice | 1e14–1e19, 61 values | m^-3, logarithmic |
| Layer thickness | 0.1 / 1 / 10 | mm |
| Equivalent flat-beam radius | 25 / 100 / 500 | micrometre |
| Screening conditions | phi<=1e-3, d>=10a, tau_sca<=0.1 | dispersion/first-pass assumptions |
| Geometry/overlap screens | L<=pi*n*b^2/lambda; d>=2a+4*sqrt(Dw*200ns) | approximate collimation; short-time thermal separation |
| Water domain | 180 cells; outer radius 3 | micrometre; compared with 90/360/720 cells and 6 micrometre bath |

Laser, geometry and concentration bounds are finite design choices, not one
experiment or universal limits. Material constants and relaxation scales are
literature-grounded. The conservative radii reduce small-particle film-to-sphere
uncertainty without claiming it vanishes. Parameter provenance is explicit.

## Resonances, scattering and distinct optical objectives

For water, at radii 5/10/20/40/60/80/100 nm, Mie absorption maxima occur at
522/522/526/541/517/524/506 nm. Corresponding extinction maxima are
522/522/527/550/591/658/561 nm. Multipolar competition makes absorption peaks
nonmonotone; absorption and extinction resonances are not interchangeable.
At the absorption maxima, scattering fractions are 0.10/0.79/6.34/39.70/47.30/
54.61/51.95 percent. Radiative dipole peak errors are -1.3/-5.0/-18.0/-41.4/
-31.6/-71.5/-82.9 percent. It stays passive but is quantitatively insufficient
for common larger spheres.

The full optical-grid maximum single cross section is 5.30938e-14 m^2 at
100 nm, 499 nm, n=1. The water maximum is 5.23979e-14 m^2 at 100 nm, 506 nm.
The maximum absorption per gold volume is 1.59417e8 m^-1 at 15 nm, 537 nm,
n=1.5; the water counterpart is 1.06197e8 m^-1 at 21 nm, 526 nm.
At fixed number density volumetric absorption follows sigma_abs; at fixed gold
volume fraction it follows sigma_abs / particle volume. These objectives and
finite-beam absorbed fraction are explicitly separated.

## Thermal results and strongest moderate candidates

Reference: a=20 nm, 532 nm, F=1 J/m^2, 10 ns, G=100 MW/(m^2 K).
Absorbed energy 3.43373 fJ; gold peak rise 2.31891 K; liquid-surface peak
1.64060 K near pulse end. Net transfer to liquid is 94.36% by pulse end,
99.57% after another 10 ns, 99.98% after another 100 ns.
Au capacity 8.35171e-17 J/K; Au internal diffusion 3.145 ps, water diffusion
2.791 ns, interface relaxation 0.166 ns. Gold and liquid maxima are evaluated
separately; short pulses can give delayed water maxima. The mathematical
instantaneous limit is not interpreted as a physical one-temperature fs model.
The analytic CW control and interface jump are checked separately.

The largest one-particle moderate energy in the stated design is at a=100 nm,
506 nm, 100 ns, F=60.7279 J/m^2, I_peak=6.07279e8 W/m^2:
E_abs=3.182015 pJ, Au rise=30 K, net liquid energy=3.137235 pJ (98.5927%)
by pulse end plus 100 ns. This lies at the upper radius and duration boundaries.
The finite thermal limit scales with the selected temperature screen while
laser limits are inactive; no universal maximum follows.

The separate maximum-energy ensemble is a=20 nm, 613 nm, N=1.46780e17 m^-3,
L=10 mm, b=500 micrometre, 100 ns, F=1000 J/m^2. Incident energy is 785.398
microjoule; absorbed 235.051 microjoule (29.9276%), first-event scattered
62.020 microjoule (7.8966%), transmitted 488.328 microjoule (62.1758%).
Tau_ext=0.475204, tau_sca=0.099209, phi=4.91864e-6, spacing=1.89574
micrometre. Entrance Au rise=25.586 K and liquid-surface peak=15.396 K for
G=50 MW/(m^2 K). Mean layer sensible rise is only 0.007148 K. Net liquid
energy after another 100 ns is 235.015 microjoule.

The distinct fraction maximum is 90.2007% at a=20 nm, 497 nm,
N=1.21153e17 m^-3, L=10 mm; the b=100 micrometre representative uses
F=137.358 J/m^2 and absorbs 3.89237 microjoule. Tau_ext=2.74904 and
phi=4.05987e-6. This is a single-pass model; the scattered fraction can be
redirected in a real sample. Fraction is independent of admitted beam area.

## Formal extremes, upper-bound interface and limits

The expanded thermal scan retains four radii, six durations from 0.1 ps to
100 ns and 37 fluences from 1e-3 to 1e3 J/m^2 (888 cases). At 20 nm,
526 nm, 0.1 ps and F=1000 J/m^2, cold optics gives 3.55514 pJ and an
adiabatic constant-cp sensible-equivalent Au rise about 42568 K; the formal
one-temperature continuation gives about 42556 K. At 100 nm the largest
formal pulse absorption is 52.3979 pJ. These are not attained temperatures or
reliable hot-material absorption predictions. Electron/lattice disequilibrium,
melting/reshaping, evaporation, changing dielectric response, latent heat,
liquid nucleation and moving boundaries are absent. The 1337.33 K bulk gold
point and 373.15 K atmospheric water diagnostic are not nanoparticle phase caps.
NIST's in-range Au vapor-pressure fit is a diagnostic only, not an ablation model.

`src/plasmonics.py` documents all reusable SI interfaces, with no bubble imports.
Its ideal energy-transfer upper-bound scenario assigns selected absorbed energy
to a chosen liquid volume with unit efficiency by assumption. For a 6 micrometre
radius volume the single-particle moderate energy gives about 0.000840 K water
sensible-equivalent; the macroscopic layer energy gives about 62046 K.
The latter is a concentration/energy bookkeeping scenario, NOT a real liquid or
gas temperature, not the computed diffusion profile and not bubble efficiency.
No mechanism brings the macroscopic absorption into that small volume here.

Timing comparison is read-only: inherited f=30 kHz, period=33.33 microsecond;
one saved complete MR4 1.15 bar saturated-air coupled case has 17.85 ns thermal
half-height width. A 100 ns later liquid-transfer total need not be available
during that peak. No acoustic coupling or new gas trajectory is computed.

## Verification and completeness

- 18 automatic tests pass, including positive absorption, optical balance,
  zero-loss/index-matched controls, dipole limit, interpolation bounds,
  finite-beam conservation, exact thermal limits and admissible maximum selection.
- Extinction balance maximum relative residual: 7.902e-16.
- Independent miepython absorption difference: 8.457e-11.
- Eight-extra-order Mie convergence: 8.127e-13.
- Linear/PCHIP spectral difference up to 6.46%; peak value difference up to 4.37%.
  A 58 nm competing-peak switch is explicitly retained, not hidden as solver error.
- Dipole-limit relative difference at 0.1 nm: 7.125e-6, mathematical test only.
- 180 vs 720 cell moderate Au-peak difference: 8.031e-04;
  liquid-surface peak: 2.240e-04;
  liquid transferred-energy difference: 2.106e-07.
- Thermal energy residual: 1.929e-12;
  independent BDF gold comparison: 2.547e-10;
  interface-flux quadrature: 1.554e-15.
- Finite-beam partition error: 4.441e-16;
  independent resolved-depth quadrature difference: 1.110e-16.

Fresh locked-environment reproduction passed: 18 tests; all 19 CSV tables,
all arrays in both NPZ files, result JSON values, eight JPG figures and extracted
PDF text reproduced exactly in a separate archive and new environment. All nine
PDF pages were visually inspected; no clipping, overlap or unreadable equations
were found. The PDF was rebuilt successfully. Full package integrity and staged
archive checks passed before dispatch; the final counts are supplied by the
dispatch and saved manifests. See `clean_reproduction.json` and `validation.md`.

## Unresolved scientific questions for independent review

1. Adequacy of uncorrected bulk-film optics for each selected radius, ligand and
   temperature; meaningful material uncertainty beyond interpolation sensitivity.
2. Experimental calibration of absorption/thermal transfer, exact beam profile,
   solvent/entry losses and particle size/morphology distribution.
3. Whether independence, tau_sca and thermal-overlap screens justify this
   approximation for a real dispersion; diffuse reabsorption and rare aggregates.
4. Temperature-dependent interface conductance/transport, the 30 K design screen,
   electron/lattice times and the uniform-Au approximation at the boundary cases.
5. Actual reshaping, melting, vapor-nucleation and destruction thresholds; replacing
   formal continuations with enthalpy/phase-change and two-temperature physics.
6. Localization, capture, timing and energy transfer to the acoustic interface,
   reserved for MR6; no coupling efficiency is inferred in MR5.
7. Global versus discrete design optima, competing spectral peaks, parameter
   bounds and the difference between energy, fraction and temperature objectives.

Please independently audit the derivations and units, optical conventions and
passivity, Mie and thermal convergence, depth/beam/heat balances, all maximum
claims, physical applicability labels, source interpretation and report/data
consistency. Identify required corrections before independent acceptance.
