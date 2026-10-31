# Saved data and units

Every numeric API uses SI. Plot labels may use SI prefixes (nm, ns, pJ, μJ).
All temperatures are rises except explicitly named ambient/material reference
temperatures. Real host indices and optical depths are dimensionless.

- `optical_grid.npz`: axes `n_medium` (3), `radius_m` (96), `wavelength_m` (601).
  Arrays `sigma_abs_m2`, `sigma_sca_m2`, `sigma_ext_m2`,
  `sigma_abs_dipole_m2` have shape `(3,96,601)`; units m^2. No interpolation beyond
  the archived data interval and no filtering of the finite design.
- `cross_sections.csv`: samples of the importable optical function. Extra
  laboratory samples at 1064 nm are within the measured data interval but outside
  the primary 400–1000 nm optimization grid.
- `resonances.csv`: one row per radius/host; maximum absorption and extinction
  vacuum wavelengths [m]; sections [m^2]; albedo and dipole relative error at the
  **Mie absorption** peak. `sigma_abs_per_gold_volume_m1` is absorption divided
  by particle volume, useful at fixed gold volume fraction.
- `thermal_response.csv`: each radius, G [W/m^2/K], duration [s]. Gold peak
  response is K/J; liquid transfer fractions refer to final absorbed pulse
  energy at pulse end or 10/100 ns later. Heat capacity [J/K], diffusion/RC times
  [s], resistances [K/W], Kapitza ratio dimensionless. G=50/100/200 MW/m^2/K.
- `pulse_reference.npz`: `time_s`, `radius_water_m`; gold rises [time], water
  rises [water cell,time]; input/storage/bath/net-liquid/balance fields in J.
  `energy_input_J` is absorbed **so far**; the final pulse energy is sigma_abs
  times reference fluence. Bath flux is integrated independently.
- `cw_heating.csv`: a and vacuum wavelength [m], local intensity [W/m^2], power
  [W], water-surface and Au steady rises [K]. CW energy requires an exposure time.
- `bounded_single_particle.csv`: best fluence/wavelength for each radius and
  duration subject to laser and 30 K gold-rise constraints at G=50 MW/m^2/K.
  `conservative_geometry` admits 20–100 nm. Energies in J, times in s. Ties at
  the thermal cap prefer larger cross section.
- `bounded_ensemble.csv`: separately optimized `fraction` and `energy` rows
  for each admitted radius/thickness/beam-radius combination. N [m^-3], geometry
  [m], duration [s], fluence [J/m^2], energies [J], spacing [m], and dimensionless
  fractions/optical depths/gold volume fraction. All rows satisfy the stated
  screens. `summary.json` selects the global winners of the declared lattice.
- `concentration.csv`: 1 J/m^2, 532 nm, L=1 mm, b=100 μm for 20/40/80 nm radii;
  full N scan includes invalid independent/single-pass regimes on purpose.
  Scattering denotes the first interaction channel, not guaranteed escape.
- `formal_regimes.csv`: all 4 radii × 6 durations × 37 fluences (888 rows).
  One-temperature and 30 K screens are Boolean diagnostics, not data filters.
  Gold and water peak rises [K] are **formal** linear continuations. Bulk melting
  and atmospheric boiling diagnostics do not assert actual damage/nucleation.
  The water peak time is independent of amplitude in this linear model.
- `ideal_upper_bound.csv`: explicitly chosen source energy [J], selected liquid
  radius [m] and volume [m^3], sensible-equivalent rise [K]. Not a gas prediction.
- `independent_mie.csv`, `mie_convergence.csv`, `dipole_limit.csv`: relative
  differences (dimensionless), sizes/wavelengths in m. Absorption is never clipped.
- `thermal_convergence.csv`: relative differences against 720 water cells;
  peak times [s], scaled energy residuals. 0.1 ns runs may fail physical time
  screens despite converging numerically.
- `interpolation_sensitivity.csv`: linear vs PCHIP differences and peak shifts
  [m]. This is input-model sensitivity, not a solver error bound.
- `optical_refinement.csv`: local fine-radius and optimized-wavelength controls
  around each host's coarse absorption maximum; not a new global continuous scan.
- `beam_energy_checks.csv`: depth quadrature and exact partition fractions.
  Very large test concentrations check algebra only, not optical independence.
- `literature_comparison.csv`: reported **calculated/estimated** reference values
  with units in the reference name; no literature curve is labeled experimental.
- `parameters.csv/json`, `environment.json`, `validation_summary.json`,
  `input_inventory.json`, `clean_reproduction.json`: full explicit design,
  software versions, numerical checks, byte preservation and clean-run evidence.

`FILE_LIST.txt` includes itself and SHA256SUMS. The checksum manifest hashes every
listed file except itself. It intentionally excludes virtual environments,
runtime caches and intermediate LaTeX files. Metadata lives outside scientific
report text. Reproduction can change PDF timestamps and test durations while
preserving numeric outputs.
