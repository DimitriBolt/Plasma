# Data schema and interpretation

All computations use SI. CSV values retain float64 reporting precision; report
rounding does not claim comparable physical accuracy. Every old tracked input
is fingerprinted in `accepted_tree_sha256.json` (1491 files).

- `mechanics.csv`: 18 fresh one-cycle integrations. `Pa_bar` is amplitude in bar;
  radii are m, speeds m/s, times s, pressures Pa; `chi` and Mach are dimensionless.
  `max_M_local` and `max_density_ratio` apply only to Gilmore. Identical MR2
  gas law, reference state and forcing are used for both models. No temperature.
- `mechanics_{model}_{amplitude}.npz`: strict-control dense/event-aware traces,
  arrays `t_s`, `R_m`, `U_m_s`. No interpolation beyond the solved cycle.
- `convergence.csv`: discrepancies from the strict DOP853 control at the same
  amplitude/model, including independent Radau. Empirical, not rigorous bounds.
- `inherited_physical_diagnostics.csv`: exact values read from MR7, with explicit
  complete/censored status. Half-height widths are reconstructed thermal widths,
  not ion, emission or nuclear confinement lifetimes.
- `energy_audit.csv`: actual saved endpoint storage, work, absorption, incident
  energy and signed integrated powers, plus independently reassembled residuals.
  Energies are J. Water/Au values are relative to the reference state. Acoustic
  background is included; they are not all laser-tagged energies. 1.30 bar has a
  shorter censored endpoint.
- `paired_energy_audit.csv`: same-prefix differences. Caloric delta-U is newly
  integrated from the saved tiny delta-T without subtracting two large energies.
  Work and heat increments come directly from the archived nonlinear ledgers.
  The stored total paired residual refers to original independent state ledgers;
  full paired water states are not reconstructed or claimed to be available.
- `ideal_caloric_allowance.csv`: invert ∫Cv(T)dT=Eabs at fixed volume and fixed
  inventory, with all absorption assigned directly to gas. This defined operation
  excludes transport and further work. It is not a dynamic/universal bound.
  Fluences 1, 10, 100, 1000 J/m² are algebraic scenarios; higher fluences do not
  acquire material admissibility from this calculation.
- `particle_kinematics.csv`: imposed constant-gap velocity minus passive-water
  velocity, and particle Reynolds number. No trapping force is derived. The
  Stokes relaxation value is a low-Re reference, not a strong-collapse drag law.
- `scales.json`: reference material/time/energy scales and their definitions.
- `formal_extrema.json`: selected MR3 extrema, keeping observed/censored versus
  complete scope. JSON NaN denotes an inherited undefined quantity, never zero.
- `scenarios.csv`: full parameter/status/applicability matrix for representative
  scenarios, including coefficient domains and missing physical information.
- `source_inventory.json`, `figure_provenance.json`: hashes and source paths.
- `verification_values.json`: machine-readable numerical maxima used in reports.

Gas wet-air inventory is about 1.86% H2O; molecular heat capacities are NASA9
fixed-species data. H2O domain 200–6000 K, N2/O2 200–20000 K; analytic ground-state
Ar has no polynomial ceiling but does not model excitation/ionization.
No chemistry or coefficient extrapolation is performed. Reference atomic or
phase energies do not determine composition, ionization fraction or temperature.
